/*
 * Processing the PesTrie index to serve queries.
 *
 * Initiated 2011.8
 * Improved 2012.7: improve the alias query to O(lgn) practical performance.
 * Improved 2012.10: fix bugs and code refactoring.
 */

#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include <ctime>
#include <set>
#include "query.h"
#include "profile_helper.h"
#include "kvec.hh"
#include "histogram.hh"
#include "bitmap.h"

using namespace std;

// Point
class Point
{
public:
  int y1;

  Point() {}
  Point( const Point& pt ) { y1 = pt.y1; }
  Point( int y ): y1(y) {}
  
  Point& operator=( const Point& other )
  {
    y1 = other.y1;
    return *this;
  }

  static void init_buffer( int sz ) { 
    Point::pt_buf = new Point[sz]; 
    Point::pt_buf_cur = 0; 
  }
  static Point* get_point(int y1) { 
    Point* pt = Point::pt_buf + (Point::pt_buf_cur++); 
    pt->y1 = y1; 
    return pt; 
  }
  static void release_buffer() { 
    delete[] Point::pt_buf; 
  }

private:
  static Point* pt_buf;
  static int pt_buf_cur;
};

Point* Point::pt_buf = NULL;
int Point::pt_buf_cur = 0;

// Vertical Line
class VLine : public Point
{
public:
  int y2;

  VLine() { }
  VLine( int y1, int y2 ) { this->y1 = y1; this->y2 = y2; }
  VLine( const VLine* other ) { this->y1 = other->y1; this->y2 = other->y2; }
  VLine& operator=( const VLine& other )
  {
    y1 = other.y1;
    y2 = other.y2;
    return *this;
  }

  static void init_buffer( int sz ) { 
    VLine::vline_buf = new VLine[sz]; 
    VLine::vline_buf_cur = 0; 
  }
  static VLine* get_vline(int y1, int y2) { 
    VLine* vl = VLine::vline_buf + (VLine::vline_buf_cur++); 
    vl->y1 = y1; vl->y2 = y2; 
    return vl; 
  }
  static void release_buffer() { 
    delete[] VLine::vline_buf; 
  }
  
private:
  static VLine* vline_buf;
  static int vline_buf_cur;
};

VLine* VLine::vline_buf = NULL;
int VLine::vline_buf_cur = 0;


class QHeader
{
public:
  kvec_t(Point*) points;
  kvec_t(VLine*) rects;
  bitmap pointsto;

  QHeader()
  {
    kv_init( Point*, points);
    kv_init( VLine*, rects );
    pointsto = BITMAP_ALLOC(NULL);
  }  
};

/*
 * Cut off the rectangles into strip pieces.
 */
static QHeader** unitRoots;

static int index_type;
static int max_store_prev;

// The number of pointers and objects
static int n, m, n_trees, n_es;
static int vertex_num, n_rects, n_verticals, n_horizontals, n_points;
// Mapping from pointer ID to tree ID
static int *tree;
// Mapping from pointer ID to PesTrie ID
static int *preV;
//
static int *obj_rank;

// Mapping from PesTrie node ID to tree ID
static int *prev_to_tree;
// Used for testing if a pre-order belongs to a root
static int *rootPrevs;      

//
static kvec_t(int) *es2pointers;
static kvec_t(int) *es2objs;

//Statistics
static int cnt_same_tree = 0;


static QHeader* get_root(int x)
{
  struct QHeader *p = unitRoots[x];
  if ( p == NULL ) {
    p = new QHeader;
    unitRoots[x] = p;
  }
  return p;
}


static void
insert_horizs( int x1, int x2, Point* pt )
{
  QHeader *p;
  
  for ( int x = x1; x <= x2; ++x ) {
    p = get_root(x);
    p->points.push_back(pt);
  }
}


static void
insert_verticals( int x, VLine* pl )
{
  QHeader *p;

  p = get_root(x);
  p->rects.push_back(pl);
}


static void
insert_rects( int x1, int x2, VLine* pv )
{
  QHeader *p;
  
  for ( int x = x1; x <= x2; ++x ) {
    p = get_root(x);
    p->rects.push_back(pv);
  }
}


static void
insert_points( int x, int y )
{
  Point* pt = Point::get_point(y);
  get_root(x)->points.push_back(pt);
}

static bool compare_points( Point* r1, Point* r2 )
{
  return r1->y1 < r2->y1;
}

static void 
build_query_header( int x, QHeader *p )
{
  // First sort all these objects
  kvec_t(Point*) &points = p->points;
  kvec_t(VLine*) &rects = p->rects;
  bitmap pointsto = p->pointsto;

  sort( points.begin(), points.end(), compare_points );
  sort( rects.begin(), rects.end(), compare_points );

#define VISIT_OBJ(v)				\
  do{						\
    if ( rootPrevs[v] != 0 ) {			\
      int tr = prev_to_tree[v];			\
      bitmap_set_bit(pointsto, tr);		\
    }						\
  }while (0)

  int size = points.size();
  for ( int i = 0; i < size; ++i ) {
    Point* pt = points[i];
    int v = pt->y1;
    VISIT_OBJ(v);
  }
  
  size = rects.size();
  for ( int i = 0; i < size; ++i ) {
    VLine *r = rects[i];
    int v = r->y1;
    VISIT_OBJ(v);
  }

  // last, x -> tree[x]
  bitmap_set_bit( pointsto, tree[x] );
}

// Concatenate continuous segments into large segments
// Be careful some segments are shared
/*
static void
merge_segments( VECTOR(VLine*) &rects )
{
  int size = rects.size();
  VLine* p_last = rects[0];
  int i_last = 0;
  int k = 1;

  for ( int i = 1; i < size; ++i ) {
    VLine* pl = rects[i];
    
    if ( pl->y1 > p_last->y2 + 1 ) {
      // They are not continuous
      rects[k] = pl;
      p_last = pl;
      i_last = k++;
    }
    else {
      // merge
      if ( p_last->counts > 1 ) {
	// We copy on write
	VLine* p_new = new VLine(p_last);
	p_last->counts--;
	rects[i_last] = p_new;
	p_last = p_new;
      }
      p_last->y2 = pl->y2;
      
      // Update
      pl->counts--;
      if ( pl->counts == 0 ) delete pl;
    }
  }
  
  //if ( k < size ) printf( "%d\n", size - k );
  rects.reset_end(k);
}
*/

static bool 
read_index()
{
  char magic_code[8];
  FILE *fp;

  fp = fopen( input_file, "rb" );
  if ( fp == NULL )
    return false;

  fread( magic_code, sizeof(char), 4, fp );
  magic_code[4] = 0;
  index_type = UNDEFINED_MATRIX;

  if ( strcmp( magic_code, PESTRIE_PT_1 ) == 0)
    index_type = PT_MATRIX;
  else if ( strcmp( magic_code, PESTRIE_SE_1 ) == 0 )
    index_type = SE_MATRIX;

  if ( index_type == UNDEFINED_MATRIX ) {
    fprintf( stderr, "This is an INVALID PesTrie index file.\n" );
    return false;
  }
  
  fread( &n, sizeof(int), 1, fp );
  fread( &m, sizeof(int), 1, fp );
  fread( &vertex_num, sizeof(int), 1, fp );
  fread( &n_rects, sizeof(int), 1, fp );
  fread( &n_verticals, sizeof(int), 1, fp );
  fread( &n_horizontals, sizeof(int), 1, fp );
  fread( &n_points, sizeof(int), 1, fp );

  Point::init_buffer(n_points*2 + n_verticals + n_horizontals);
  VLine::init_buffer(n_rects*2 + n_verticals + n_horizontals);

  // Now initialize the data structures
  ++vertex_num;
  tree = new int[n+m];
  preV = new int[n+m];
  obj_rank = new int[m+1];
  es2pointers = new kvec_t(int)[vertex_num];
  es2objs = new kvec_t(int)[vertex_num];

  // Read in the pre-order descriptors for both pointers and objects
  fread( preV, sizeof(int), n+m, fp );

  // Now we rebuild the object permutation
  n_trees = 0;
  rootPrevs = new int[vertex_num];
  prev_to_tree = new int[vertex_num];
  memset ( rootPrevs, 0, sizeof(int) * vertex_num );
  
  for ( int i = 0; i < m; ++i ) { 
    int v = preV[n+i];
    if ( v != -1 ) {
      if ( rootPrevs[v] == 0 ) {
	obj_rank[n_trees++] = v;
	rootPrevs[v] = 1;
      }
      // Object i belongs to equivalent set v
      es2objs[v].push_back(i);
    }
  }
  
  // Objects with smaller preV must be earlier roots
  sort( obj_rank, obj_rank + n_trees );
  for ( int i = 0; i < n_trees; ++i ) {
    int v = obj_rank[i];
    prev_to_tree[v] = i;
  }

  // Sentinels
  obj_rank[n_trees] = vertex_num;
  if ( index_type == SE_MATRIX ) max_store_prev = obj_rank[m/2];

  // We re-discover the tree codes for pointers through the binary search
  n_es = 0;
  for ( int i = 0; i < n; ++i ) {
    int preI = preV[i];
    if ( preI > n_es ) n_es = preI;

    if ( preI != -1 ) {
      int s = 0, e = n_trees;
      int mid;
      /*
       * Invariant: s is the larget root <= pre-label for i
       * Since multiple roots may have the same label, s is the latest one.
       */
      while ( e - s > 1 ) {
	mid = (s + e) / 2;
	if ( obj_rank[mid] <= preI )
	  s = mid;
	else
	  e = mid;
      }
      
      // Assign the tree code
      tree[i] = s;
      es2pointers[preI].push_back(i);
    }
    else
      tree[i] = -1;
  }

  // #of ES groups = largest label + 1
  ++n_es;

  // Now we read in the index figures
  unitRoots = new QHeader*[vertex_num];
  memset( unitRoots, 0, sizeof(void*) * vertex_num );

  int max_buf = 0;
  if ( n_rects*4 > max_buf ) max_buf = n_rects*4;
  if ( n_verticals*3 > max_buf ) max_buf = n_verticals*3;
  if ( n_horizontals*3 > max_buf ) max_buf = n_horizontals*3;
  if ( n_points*2 > max_buf ) max_buf = n_points*2;

  int *labels = new int[max_buf];
 
  // Rectangles
  fread( labels, sizeof(int), 4*n_rects, fp );
  for ( int i = 0, k = 0; i < n_rects; ++i, k+=4 ) {
    VLine* pv1 = VLine::get_vline(labels[k+1], labels[k+3]);
    VLine* pv2 = VLine::get_vline(labels[k], labels[k+2]);

    // Original
    insert_rects( pv2->y1, pv2->y2, pv1 );
    // Inverse
    insert_rects( pv1->y1, pv1->y2, pv2 );
  }
  
  // Vertical lines
  fread( labels, sizeof(int), 3*n_verticals, fp );
  for ( int i = 0, k = 0; i < n_verticals; ++i, k+=3 ) {
    int x = labels[k+1];
    Point* pt = Point::get_point(x);
    VLine* pv = VLine::get_vline(labels[k], labels[k+2]);

    // Original
    insert_verticals( x, pv );    
    // It's inverse
    insert_horizs(pv->y1, pv->y2, pt);
  }

  // Horizontal lines
  fread( labels, sizeof(int), 3*n_horizontals, fp );
  for ( int i = 0, k = 0; i < n_horizontals; ++i, k+=3 ) {
    int x = labels[k+1];
    Point* pt = Point::get_point(x);
    VLine* pv = VLine::get_vline( labels[k], labels[k+2] );

    // Original
    insert_horizs(pv->y1, pv->y2, pt);
    // Inverse
    insert_verticals( x, pv );
  }
  
  // Points
  // special
  while ( n_points > 0 ) {
    // Read the X label and #of shared points
    int count;
    fread( labels, sizeof(int), 1, fp );    
    fread( &count, sizeof(int), 1, fp );
    fread( labels + 1, sizeof(int), count, fp );
    n_points -= count;
    
    for ( int k = 0; k < count; ++k ) {
      int x = labels[0];
      int y = labels[k+1];
      insert_points( x, y );
      insert_points( y, x );
    }
  }
  
  fclose( fp );

  // Sort and build up the querying structure
  for ( int i = 0; i < vertex_num; ++i )
    if ( unitRoots[i] != NULL )
      build_query_header( i, unitRoots[i] );

  fprintf( stderr, "\n-------Input: %s-------\n", input_file );
  show_res_use( "Index loading" );
  
  delete[] labels;
  delete[] rootPrevs;
  
  return true;
}

static void 
profile_pestrie()
{
  fprintf( stderr, "Trees = %d, Nodes = %d\n", n_trees, vertex_num - 1 );
}

static int
iterate_equivalent_set( kvec_t(int) *es_set )
{
  int ans = 0;
  int size = es_set->size();

  for ( int i = 0; i < size; ++i ) {
    int q = es_set->at(i);
    if ( q >= 0 )
      ++ans;
  }

  return ans;
}

static bool
IsAlias( int x, int y )
{
  int tr1, tr2;
  struct QHeader *p;
  int yy;
  int mid, s, e;

  tr1 = tree[x];
  if ( tr1 == -1 ) return false;
  tr2 = tree[y];
  if ( tr2 == -1 ) return false;
  if ( tr1 == tr2 ) {
    ++cnt_same_tree;
    return true;
  }

  x = preV[x];
  p = unitRoots[x];
  if ( p == NULL ) return false;
  y = preV[y];

  // Now we use binary search for static point location query
  kvec_t(Point*) &points = p->points;
  kvec_t(VLine*) &rects = p->rects;
  
  // Points
  s = 0; e = points.size();
  while ( e > s ) {
    mid = (s+e) / 2;
    Point* pt = points[mid];
    int y1 = pt->y1;
    if ( y1 >= y ) {
      if ( y1 == y ) return true;
      e = mid;
    }
    else
      s = mid + 1;
  }
  
  // Rectangles
  s = 0; e = rects.size();
  while ( e > s ) {
    mid = (s+e) / 2;
    VLine* tt = rects[mid];
    if ( tt->y2 >= y ) {
      if ( tt->y1 <= y ) {
	// Found the closest one
	return true;
      }
      e = mid;
    }
    else
      s = mid + 1;
  }
    
  return false;
}

// List query in real use should be passed in a handler.
// That handler decide what to do with the query answer.
// This is to avoid the compiler smartly deleting the code.
static int 
ListPointsTo( int x )
{
  QHeader *p;
  int ans = 0;

  x = preV[x];
  if ( x != -1 ) {
    p = unitRoots[x];
    if ( p != NULL ) {
      bitmap pointsto = p->pointsto;
      unsigned o;
      bitmap_iterator bi;
      EXECUTE_IF_SET_IN_BITMAP( pointsto, 0, o, bi ) {
	kvec_t(int) *objs = &es2objs[o];
	ans += iterate_equivalent_set( objs );
      }
    }
  }

  return ans;
}

// Find all the pointers y that *x and *y is an alias pair
static int 
ListAliases( int x, kvec_t(int) *es2baseptrs )
{
  QHeader *p;
  int size;

  int tr = tree[x];
  if ( tr == -1 ) return 0;
  
  int ans = 0;
  x = preV[x];

#define VISIT(es)							\
    do {								\
      kvec_t(int) *ptrs = ( es2baseptrs == NULL ? &es2pointers[es] : &es2baseptrs[es] ); \
      ans += iterate_equivalent_set( ptrs );				\
    } while(0)

  // We first extract the ES groups that belong to the same subtree
  {
    int upper = obj_rank[tr+1];
    for ( int i = obj_rank[tr]; i < upper; ++i ) {
      VISIT(i);
    }
  }
  
  // traverse index figures
  p = unitRoots[x];
  if ( p != NULL ) {
    kvec_t(Point*) &points = p->points;
    kvec_t(VLine*) &rects = p->rects;

    // Points
    size = points.size();
    for ( int i = 0; i < size; ++i ) {
      Point* pt = points[i];
      VISIT(pt->y1);
    }

    // Rects
    size = rects.size();
    for ( int i = 0; i < size; ++i ) {
      VLine *r = rects[i];
      int lower = r->y1;
      int upper = r->y2;
      
      for ( ; lower <= upper; ++lower ) {
	VISIT(lower);
      }
    }
  }
  
  return ans;
}

static int 
ListPointedTo( int o )
{
  return ListAliases( o + n, NULL );
}

static int 
ListModRefVars( int x )
{
  return ListPointsTo( x );
}

static int
ListConflicts( int x )
{
  int ans = 0;

  if ( preV[x] < max_store_prev )
    ans = ListAliases( x, NULL );
  
  return ans;
}

void
execute_query_plan()
{
  int x, y;

  FILE *fp = fopen( query_plan, "r" );
  if ( fp == NULL ) {
    fprintf( stderr, "Cannot open the query plan file. Simulation exits.\n" );
    return;
  }

  // Read
  kvec_t(int) queries;
  kv_init( int, queries );
  kvec_t(int) *es2baseptrs = new kvec_t(int)[vertex_num];
    
  while ( fscanf( fp, "%d", &x ) != EOF ) {
    queries.push_back(x);
    int es = preV[x];
    if ( es >= vertex_num ) printf( "oops\n" );
    if ( es != -1 )
      es2baseptrs[es].push_back(x);
  }
  
  fclose( fp );
  n_query = kv_size( queries );  
  //fprintf( stderr, "Query plan loaded : %d entries.\n", n_query );
  show_res_use( NULL );

  // Execute
  for ( int i = 0; i < n_query; ++i ) {
    x = queries[i];
    
    switch ( query_type ) {
    case IS_ALIAS:
      // No difference in PesTrie cass
      for ( int j = i + 1; j < n_query; ++j ) {
	y = queries[j];
	bool ans = IsAlias( x, y );
	if ( print_answers )
	  printf( "(%d, %d) : %s\n", x, y, ans == true ? "true" : "false" );
      }
      break;

    case LIST_POINTS_TO:
      {
	int ans = ListPointsTo( x );
	if ( print_answers )
	  printf( "%d : %d\n", x, ans );
      }
      break;
      
    case LIST_ALIASES:
      {
	int ans = ListAliases( x, es2baseptrs );
	if ( print_answers )
	  printf( "%d : %d\n", x, ans );
      }
      break;
    }
  }

  delete[] es2baseptrs;
}

// We generate random queries
void 
simulate_queries()
{
  int x, y;
  int ans = 0;

  if ( (index_type == PT_MATRIX && query_type >= LIST_ACC_VARS) ||
       (index_type == SE_MATRIX && query_type < LIST_ACC_VARS) ) {
    fprintf( stderr, "The query commands are supported by input index file.\n" );
    return;
  }

  if ( index_type == SE_MATRIX ) n_query = n;
  //printf ( "max_store_prev = %d\n", max_store_prev );

  for ( int i = 0; i < n_query; ++i ) {
    switch (query_type) {
    case IS_ALIAS:
      x = rand() % n; y = rand() % n;
      IsAlias( x, y );
      break;
      
    case LIST_POINTS_TO:
      x = rand() % n;
      ListPointsTo( x );
      break;
      
    case LIST_POINTED_TO:
      x = rand() % m;
      ListPointedTo( x );
      break;
      
    case LIST_ALIASES:
      x = rand() % n;
      ListAliases( x, NULL );
      break;
  
    case LIST_ACC_VARS:
      ans += ListModRefVars(i);
      break;
      
    case LIST_CONFLICTS:
      ans += ListConflicts(i);
      //printf ("%d\n", preV[i] );
      break;
    }
  }

  //fprintf( stderr, "Total Conflicts = %d\n", ans );
}

int main( int argc, char** argv )
{
  if ( parse_options( argc, argv ) == false )
    return -1;

  srand( time(NULL) );
  bitmap_obstack_initialize(NULL);

  if ( read_index() == false )
    return -1;

  if ( do_profile ) profile_pestrie();

  if ( query_plan != NULL )
    execute_query_plan();
  else
    simulate_queries();

  // if ( query_type == IS_ALIAS ) {
  //   long total = (long)n_query * (n_query -1);
  //   fprintf ( stderr, "Direct answers = %d, total = %lld\n", cnt_same_tree, total/2 );
  // }

  char buf[128];
  sprintf( buf, "%s querying", query_strs[query_type] );
  show_res_use( buf );

  Point::release_buffer();
  VLine::release_buffer();

  return 0;
}
