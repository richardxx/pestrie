/*
 * Processing the PesTrie index to serve queries.
 *
 * Initiated 2011.8
 * Improved 2012.7: improve the alias query to O(lgn) practical performance.
 * Improved 2012.10: fix bugs and code refactoring.
 * Improved 2014.02: fix bugs and improve performance significantly.
 */

#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include <ctime>
#include <set>
#include "query.h"
#include "pes-querier.hh"
#include "profile_helper.h"
#include "histogram.hh"
#include "bitmap.h"

using namespace std;

/*
 * Cut off the rectangles into strip pieces.
 */
static QHeader **unitRoots, *segRoot;

static int index_type;
static int max_store_prev;

// The number of pointers and objects
static int n, m, n_trees, n_es;
static int vertex_num, n_figures;
// Mapping from pointer ID to tree ID
static int *tree;
// Mapping from pointer ID to PesTrie ID
static int *preV;
// Mapping tree ID to pre-order stamp
static int *root_prevs;
// Mapping the pre-order of root node to tree ID
static int *root_tree;

//
static VECTOR(int) *es2pointers;
static VECTOR(int) *es2objs;

//Statistics
static int cnt_same_tree = 0;


static void
prepare_pestrie_info( FILE* fp )
{
  fread( &n, sizeof(int), 1, fp );
  fread( &m, sizeof(int), 1, fp );
  fread( &vertex_num, sizeof(int), 1, fp );
  fread( &n_figures, sizeof(int), 1, fp );

  tree = new int[n+m];
  preV = new int[n+m];
  root_prevs = new int[m+1];
  root_tree = new int[vertex_num+1];
  es2pointers = new VECTOR(int)[vertex_num];
  es2objs = new VECTOR(int)[vertex_num];

  // Read in the pre-order descriptors for both pointers and objects
  fread( preV, sizeof(int), n+m, fp );

  // We label the time-stamps that could be roots
  n_trees = 0;
  memset ( root_tree, -1, sizeof(int) * (vertex_num+1) );
  
  for ( int i = 0; i < m; ++i ) { 
    int v = preV[n+i];
    if ( v != -1 ) {
      if ( root_tree[v] == -1 ) {
	root_prevs[n_trees++] = v;
	root_tree[v] = 1;
      }
    }
  }
  
  // Objects with smaller preV must be earlier roots
  sort( root_prevs, root_prevs + n_trees );
  for ( int i = 0; i < n_trees; ++i ) {
    int v = root_prevs[i];
    root_tree[v] = i;
  }

  // Collect the equivalent objects
  for ( int i = 0; i < m; ++i ) { 
    int v = preV[n+i];
    if ( v != -1 ) {
      int tr = root_tree[v];
      es2objs[tr].push_back(i);
      tree[i+n] = tr;
    }
  }

  if ( index_type == SE_MATRIX ) 
    max_store_prev = root_prevs[m/2];

  // We re-discover the tree codes for pointers through the binary search
  // Sentinels
  root_prevs[n_trees] = vertex_num + 1;
  
  for ( int i = 0; i < n; ++i ) {
    int preI = preV[i];

    if ( preI != -1 ) {
      int s = 0, e = n_trees;
      int mid;
      
      while ( e - s > 1 ) {
	mid = (s + e) / 2;
	if ( root_prevs[mid] > preI )
	  e = mid;
	else
	  s = mid;
      }
      
      // Assign the tree code
      tree[i] = s;
      es2pointers[preI].push_back(i);
    }
    else
      tree[i] = -1;
  }
}

static QHeader*
build_seg_tree( int l, int r )
{
  int x = (l+r) / 2;
  QHeader* p = new QHeader; 
  
  p->x1 = l;
  p->x2 = r;
  
  if ( l != r ) {
    if ( l <= x ) {
      p->left = build_seg_tree( l, x );
      p->left->parent = p;
    }
    
    if ( x < r ) {
      p->right = build_seg_tree( x + 1, r ); 
      p->right->parent = p;
    }
  }

  unitRoots[x] = p;
  return p;
}

/*
 * We update the parent links to skip the empty nodes.
 */ 
static void
optimize_seg_tree( QHeader* p ) 
{
  QHeader* q = p->parent;
  
  if ( q != NULL &&
       q->n_of_rects() == 0 ) {
    p->parent = q->parent;
  }
  
  optimize_seg_tree( p->left );
  optimize_seg_tree( p->right );
}


static void
insert_point( int x, int y )
{
  VLine* p = new VLine(y, y);
  unitRoots[x]->add_vertis(p);
  
  p = new VLine(x, x);
  unitRoots[y]->add_vertis(p);
}

/*
 * We add log(n) references to each figure
 */
static void 
insert_rect( int x1, int x2, VLine* pr, QHeader* p )
{
  if ( x1 <= p->x1 && x2 >= p->x2 ) {
    // We only consider the full coverage
    p->add_rect(pr);
    return;
  }

  int x = (p->x1 + p->x2) / 2;
  p = unitRoots[x];
  
  if ( x1 < x ) insert_rect( x1, x2, pr, p->left );
  if ( x2 > x ) insert_rect( x1, x2, pr, p->right );
}

static bool comp_rect( VLine* r1, VLine* r2 )
{
  return r1->y1 < r2->y1;
}

static void
process_figures( FILE* fp )
{
  // First we initialize the segment tree
  unitRoots = new QHeader*[vertex_num];
  build_seg_tree( 0, vertex_num );
  segRoot = unitRoots[vertex_num/2];
  segRoot->merged = true;

  int buf_size;
  int *labels = new int[vertex_num * 3];
  VECTOR(Rectangle*) all_rects; 

  for ( int x1 = 0; x1 < vertex_num; ++x1 ) {
    fread( &buf_size, sizeof(int), 1, fp );
    if ( buf_size == 0 ) continue;
    fread( labels, sizeof(int), buf_size, fp );
    fprintf( stderr, "%d\n", buf_size );

    int i = 0;
    while ( i < buf_size ) {
      int y1 = labels[i++];
      int x2, y2;

      if ( (y1&SIG_RECT) == 0 ) {
	// This is a point, directly insert it
	insert_point( x1, y1);
	continue;
      }

      if ( (y1&SIG_VERTICAL) == SIG_VERTICAL ) {
	y1 ^= SIG_VERTICAL;
	y2 = labels[i++];
	x2 = x1;
      }
      else if ( (y1&SIG_HORIZONTAL) == SIG_HORIZONTAL ) {
	y1 ^= SIG_HORIZONTAL;
	x2 = labels[i++];
	y2 = y1;
      }
      else {
	y1 ^= SIG_RECT;
	x2 = labels[i++];
	y2 = labels[i++];
      }

      // First, the reversed rect
      VLine* p = new VLine(x1, x2);
      if ( y1 == y2 ) {
	unitRoots[y1]->add_vertis(p);
      }
      else {
	insert_rect( y1, y2, p, segRoot );
      }

      // Second, cache it
      Rectangle* r = new Rectangle(x1, x2, y1, y2);
      if ( x1 == x2 ) {
	unitRoots[x1]->add_vertis(r);
      }
      else {
	all_rects.push_back(r);
      }
    }
  }

  // We insert the cached rectangle
  sort( all_rects.begin(), all_rects.end(), comp_rect );
  int size = all_rects.size();

  for ( int i = 0; i < size; ++i ) {
    Rectangle* r = all_rects[i];
    insert_rect( r->x1, r->x2, r, segRoot );
  }

  optimize_seg_tree( segRoot );
  delete[] labels;
}

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
  
  // Loading and decoding the persistence file
  prepare_pestrie_info(fp);
  process_figures(fp);

  fprintf( stderr, "\n-------Input: %s-------\n", input_file );
  show_res_use( "Index loading" );

  fclose( fp );  
  return true;
}

static void 
profile_pestrie()
{
  int non_empty_nodes = 0;
  int n_pts = 0, n_rects = 0;

  for ( int i = 0; i < vertex_num; ++i ) {
    // Test 1:
    VECTOR(int) ptrs = es2pointers[i];
    if ( ptrs.size() > 0 )
      ++non_empty_nodes;

    // Test 2:
    struct QHeader* p = unitRoots[i];
    n_rects += p->rects.size();
    n_pts += p->vertis.size();
  }
  
  n_pts += non_empty_nodes;
  
  fprintf( stderr, "Trees = %d, ES = %d, Non-empty ES = %d\n", 
	   n_trees, vertex_num, non_empty_nodes );
  fprintf( stderr, "Number of rects = %d\n", 
	   n_rects + n_pts ); 
}

static int
iterate_equivalent_set( VECTOR(int) *es_set )
{
  int ans = 0;
  int size = es_set->size();

  for ( int i = 0; i < size; ++i ) {
    //    int q = es_set->at(i);
    //    if ( q > 0 )
      // This is to avoid the compiler smartly deleting the code.
      ans++;
  }

  return ans;
}

// We do merging sort
static void
recursive_merge( QHeader* p )
{
  if ( p->merged == true ) return;
    
  if ( p->parent->merged == false )
    recursive_merge( p->parent);

  VECTOR(VLine*) &list1 = p->parent->rects;
  int sz1 = list1.size();
  int sz2 = p->rects.size();

  if ( sz1 != 0 ) {
    // We have something to merge top down
    if ( sz2 == 0 ) {
      // Fast path, just copy
      p->rects.copy( list1 );
    }
    else {
      VECTOR(VLine*) &list3 = p->rects;
      VECTOR(VLine*) list2(list3);
      
      int i = 0, j = 0;
      VLine *r1 = list1[0], *r2 = list2[0];
      list3.clear();

      while ( r1 != NULL && r2 != NULL ) {
	if ( r2 == NULL || 
	     ( r1 != NULL && r1->y1 < r2->y1 ) ) {
	  list3.push_back(r1);
	  ++i;
	  if ( i >= sz1 ) r1 = NULL;
	  else r1 = list1[i];
	}
	else {
	  list3.push_back(r2);
	  ++j;
	  if ( j >= sz2 ) r2 = NULL;
	  else r2 = list2[j];
	}
      }
    }
  }

  p->merged = true;
}

static bool
binary_search( VECTOR(VLine*) &rects, int y )
{
  int mid, s, e;
  VLine *r, *tt;

  // Rectangles
  s = 0; e = rects.size();
  while ( e > s ) {
    mid = (s+e) / 2;
    tt = rects[mid];
    
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


static bool
IsAlias( int x, int y )
{
  struct QHeader *p;

  int tr_x = tree[x];
  if ( tr_x == -1 ) return false;
  int tr_y = tree[y];
  if ( tr_y == -1 ) return false;
  if ( tr_x == tr_y ) {
    ++cnt_same_tree;
    return true;
  }

  x = preV[x];
  y = preV[y];
  p = unitRoots[x];

  // Search the verticals special to this x axis
  if ( binary_search( p->vertis, y ) ) 
    return true;

  // We traverse the segment tree bottom up
  while ( p != NULL ) {    
    if ( binary_search( p->rects, y ) ) 
      return true;
    p = p->parent;
  }
  
  return false;
}

// List query in real use should be passed in a handler.
// That handler decide what to do with the query answer.
static int 
ListPointsTo( int x )
{
  QHeader *p;
  VECTOR(int) *objs;
  int ans = 0;

  int tr = tree[x];
  if ( tr == -1 ) return 0;
  
  // Don't forget x -> tree[x]
  objs = &es2objs[tr];
  ans += iterate_equivalent_set( objs );
  
  x = preV[x];
  p = unitRoots[x];
  
  return ans;
}

// Find all the pointers y that *x and *y is an alias pair
static int 
ListAliases( int x, VECTOR(int) *es2baseptrs )
{
  QHeader *p;
  VLine *r;
  int size;

  int tr = tree[x];
  if ( tr == -1 ) return 0;
  
  int ans = 0;
  x = preV[x];

#define VISIT_POINT(v)							\
  do {									\
    VECTOR(int) *ptrs = ( es2baseptrs == NULL ? &es2pointers[v] : &es2baseptrs[v] ); \
    ans += iterate_equivalent_set( ptrs );				\
  } while(0)

#define VISIT_RECT(r)							\
  do {									\
    int lower = r->y1;							\
    int upper = r->y2;							\
    do {								\
      VISIT_POINT(lower);						\
      ++lower;								\
    } while ( lower <= upper );						\
  } while(0)

  // We first extract the ES groups that belong to the same subtree
  {
    int upper = root_prevs[tr+1];
    for ( int i = root_prevs[tr]; i < upper; ++i ) {
      VISIT_POINT(i);
    }
  }
  
  // Visit the verticals
  VECTOR(VLine*) &vertis = p->vertis;
  size = vertis.size();
  for ( int i = 0; i < size; ++i ) {
    VLine* r = vertis[i];
    VISIT_RECT(r);
  }

  // traverse the rectangles
  p = unitRoots[x];
  while ( p != NULL ) {
    VECTOR(VLine*) &rects = p->rects;
    size = rects.size();
    for ( int i = 0; i < size; ++i ) {
      r = rects[i];
      VISIT_RECT(r);
    }
    p = p->parent;
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
  VECTOR(int) queries;
  VECTOR(int) *es2baseptrs = new VECTOR(int)[vertex_num];
    
  while ( fscanf( fp, "%d", &x ) != EOF ) {
    queries.push_back(x);
    int es = preV[x];
    if ( es >= vertex_num ) printf( "oops\n" );
    if ( es != -1 )
      es2baseptrs[es].push_back(x);
  }
  
  fclose( fp );
  int n_query = queries.size();
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
	bool aliased = IsAlias( x, y );
	if ( print_answers )
	  printf( "(%d, %d) : %s\n", x, y, aliased == true ? "true" : "false" );
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

void
traverse_result()
{
  int x, y;
  int ans = 0;

  if ( (index_type == PT_MATRIX && query_type >= LIST_ACC_VARS) ||
       (index_type == SE_MATRIX && query_type < LIST_ACC_VARS) ) {
    fprintf( stderr, "The query commands are supported by input index file.\n" );
    return;
  }

  int n_query = ( query_type == LIST_POINTED_TO ? m : n );
  
  // We permute the pointers/objects
  /*
  int *queries = new int[n_query];
  for ( int i = 0; i < n_query; ++i ) queries[i] = i;

  for ( int i = 0; i < n_query - 1; ++i ) {
    int j = queries[i];
    int k = rand() % ( n_query - i ) + i;
    queries[i] = queries[k];
    queries[k] = j;
  }
  */

  for ( int i = 0; i < n_query; ++i ) {
    //x = queries[i];
    x = i;

    switch (query_type) {
    case IS_ALIAS:
      y = n_query - x;
      ans += IsAlias( x, y ) == true ? 1 : 0;
      break;
      
    case LIST_POINTS_TO:
      ans += ListPointsTo(x);
      break;
      
    case LIST_POINTED_TO:
      ans += ListPointedTo(x);
      break;
      
    case LIST_ALIASES:
      ans += ListAliases(x, NULL);
      break;
  
    case LIST_ACC_VARS:
      ans += ListModRefVars(x);
      break;
      
    case LIST_CONFLICTS:
      ans += ListConflicts(x);
      break;
    }
  }

  fprintf( stderr, "\nReference answer = %d\n", ans );
  //delete[] queries;
}

int main( int argc, char** argv )
{
  if ( parse_options( argc, argv ) == false )
    return -1;

  srand( time(NULL) );
  bitmap_obstack_initialize(NULL);

  if ( read_index() == false ) return -1;
  if ( do_profile ) profile_pestrie();

  query_plan != NULL ? execute_query_plan() : traverse_result();
  
  // if ( query_type == IS_ALIAS ) {
  //   long total = (long)n_query * (n_query -1);
  //   fprintf ( stderr, "Direct answers = %d, total = %lld\n", cnt_same_tree, total/2 );
  // }

  char buf[128];
  sprintf( buf, "%s querying", query_strs[query_type] );
  show_res_use( buf );

  return 0;
}
