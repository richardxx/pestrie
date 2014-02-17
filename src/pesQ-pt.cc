/*
 * Processing the PesTrie index to serve queries.
 * by richardxx, 2011.8
 * modifed by richardxx, 2012.7
 */

#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include "options.h"
#include "profile_helper.h"
#include "kvec.hh"
#include "treap.hh"
#include "histogram.hh"


using namespace std;

struct segTreeNode
{
  // Start and end positions of this segment
  int s, e;
  // The middle position of this segment and its tree code
  int mid, tr_mid;
  struct segTreeNode *left, *right;

  // Group 1
  kvec_t(Rectangle*) sortByY1;
  // Group 2
  kvec_t(Rectangle*) sortByX1, sortByX2;
  // Group 3
  kvec_t(Rectangle*) sortByX1Object, sortByX2Object;

  segTreeNode()
  {
    left = NULL;
    right = NULL;
    kv_init( Rectangle*, sortByY1 );
    kv_init( Rectangle*, sortByX1 );  kv_init( Rectangle*, sortByX2 );
    kv_init( Rectangle*, sortByX1Object ); kv_init( Rectangle*, sortByX2Object );
  }  
};

struct segTreeNode *segRoot;
struct Treap_node<Rectangle>** roots;

// The number of pointers and objects
int n, m;
int vertex_num, n_rects, n_verticals, n_horizontals, n_points;
int *tree, *preV, *backPrev, *treePrev;
int *vis;          // used for indicating of the visiting trail

// Containers
int *tree_members;
vector<Rectangle*> internal_container;

//Statistics
int cnt_grp = 0;

struct segTreeNode *build_segtree( int s, int e )
{
  struct segTreeNode *r = NULL;

  if ( s <= e ) {
    r = new segTreeNode;
    r -> s = s;
    r -> e = e;
    r -> mid = ( s + e ) / 2;
    r -> tr_mid = treePrev[ r->mid ];
  }

  return r;
}

void insert_segment_tree( Rectangle* pr, bool is_points_to )
{
  
  struct segTreeNode *p;

  if ( pr->x2 - pr->x1 >= 32 ) {
    p = segRoot;
    while ( p != NULL ) {
      if ( pr->x1 <= p->mid && pr->x2 >= p->mid ) {
	// Place this figure here
	p->sortByY1.push_back( pr );
	p->sortByX1.push_back( pr );
	p->sortByX2.push_back( pr );
	if ( is_points_to ) {
	  p->sortByX1Object.push_back( pr );
	  p->sortByX2Object.push_back( pr );
	}
	break;
      }
      
      if ( pr->x1 < p->mid ) {
	if ( p->left == NULL )
	  p->left = build_segtree( p->s, p->mid-1 );
	p = p -> left;
      }
      else {
	if ( p->right == NULL )
	  p->right = build_segtree( p->mid+1, p->e );
	p = p -> right;
      }
    }
  }
  else {
    for ( int i = pr->x1; i <= pr->x2; ++i )
      roots[i] = insert_treap<Rectangle>( roots[i], pr );
  }
}

bool firstY1Small( Rectangle* r1, Rectangle* r2 )
{
  return r1->y1 < r2->y1;
}

bool firstX1Small( Rectangle* r1, Rectangle* r2 )
{
  return r1->x1 < r2->x1;
}

bool firstX2Large( Rectangle *r1, Rectangle* r2 )
{
  return r1->x2 > r2->x2;
}

void sort_segment_tree( segTreeNode *p )
{
  if ( p == NULL )
    return;

  // We first sort by Y1
  sort( p->sortByY1.begin(), p->sortByY1.end(), firstY1Small );
  // Sort by other rules
  sort( p->sortByX1.begin(), p->sortByX1.end(), firstX1Small );
  sort( p->sortByX2.begin(), p->sortByX2.end(), firstX2Large );
  sort( p->sortByX1Object.begin(), p->sortByX1Object.end(), firstX1Small );
  sort( p->sortByX2Object.begin(), p->sortByX2Object.end(), firstX2Large );

  sort_segment_tree( p->left );
  sort_segment_tree( p->right );
}

void insert_an_index_figure( Rectangle* pr )
{
  int s, e, mid;

  // We first binary search to see if this rectangle has points-to information
  s = n;
  e = n + m;
  while ( e - s > 1 ) {
    mid = (s+e) / 2;
    // pr->y1 is the preorder label for a PES root, if any
    if ( preV[mid] <= pr->y1 )
      s = mid;
    else
      e = mid;
  }

  // Insert this rectangle
  insert_segment_tree( pr, preV[s] == pr->y1 );
  
  // Then, we insert an inverse rectangle, and it cannot represent the points-to information whatever
  Rectangle* pr_prime = new Rectangle(pr->y1, pr->x1, pr->y2, pr->x2);
  insert_segment_tree( pr_prime, false );
}

bool read_index()
{
  char magic_code[8];
  FILE *fp;

  fp = fopen( input_file, "rb" );
  if ( fp == NULL )
    return false;

  fread( magic_code, sizeof(char), 4, fp );
  magic_code[4] = 0;
  if ( strcmp( magic_code, "PTP1" ) != 0 ) {
    fprintf( stderr, "This file is not a valid PesTrie index file.\n" );
    return false;
  }
  
  fread( &n, sizeof(int), 1, fp );
  fread( &m, sizeof(int), 1, fp );
  fread( &n_rects, sizeof(int), 1, fp );
  fread( &n_verticals, sizeof(int), 1, fp );
  fread( &n_horizontals, sizeof(int), 1, fp );
  fread( &n_points, sizeof(int), 1, fp );

  // Now initialize the data structures
  tree = new int[n+m];
  preV = new int[n+m];
  backPrev = new int[n+m];
  treePrev = new int[n];
  vis = new int[n+m];
  tree_members = new int[m];

  // Read in the pre-order descriptors
  fread( preV, sizeof(int), n+m, fp );

  // We re-discover the tree codes through the binary search
  // Also we recover the number of tree verteices in the PesTrie index
  vertex_num = -1;
  memset( vis, 0, sizeof(int) * (n+m) );
  memset( tree_members, 0, sizeof(int) * m );
  for ( int i = 0; i < n; ++i ) {
    int pes_node = preV[i];

    if ( pes_node != -1 ) {
      int s = n, e = n + m;
      int mid;
      
      while ( e - s > 1 ) {
	mid = (s + e) / 2;
	if ( preV[mid] <= pes_node )
	  s = mid;
	else
	  e = mid;
      }
      
      // Assign the tree code
      tree[i] = s - n;
      treePrev[pes_node] = s - n;
      // Count how many tree nodes for each tree
      if ( vis[pes_node] == 0 )	tree_members[s-n]++;
      ++vis[pes_node];
    }
    else
      tree[i] = -1;
    
    if ( pes_node > vertex_num ) 
      vertex_num = preV[i];
  }

  // The tree codes and the back mappings for the objects
  for ( int i = 0; i < m; ++i ) {
    tree[ i + n ] = i;
    backPrev[ preV[i+n] ] = i;
  }

  roots = new Treap_node<Rectangle>*[vertex_num];

  // Now we read in the index figures
  segRoot = build_segtree( 0, vertex_num );
  int labels[] = {0, 0, 0, 0};
 
  // Rectangles
  for ( int i = 0; i < n_rects; ++i ) {
    fread( labels, sizeof(int), 4, fp );
    Rectangle* pr = new Rectangle;
    pr->x1 = labels[0]; pr->y1 = labels[1]; pr->x2 = labels[2]; pr->y2 = labels[3];
    insert_an_index_figure( pr );
  }
  
  // Vertical lines
  for ( int i = 0; i < n_verticals; ++i ) {
    fread( &labels[1], sizeof(int), 3, fp );
    Rectangle* pr = new Rectangle;
    pr->y1 = labels[1]; pr->x2 = labels[2]; pr->y2 = labels[3];
    pr->x1 = pr->x2;
    insert_an_index_figure( pr );
  }

  // Horizontal lines
  for ( int i = 0; i < n_horizontals; ++i ) {
    fread( labels, sizeof(int), 3, fp );
    Rectangle* pr = new Rectangle;
    pr->x1 = labels[0]; pr->y1 = labels[1]; pr->x2 = labels[2];
    pr->y2 = pr->y1;
    insert_an_index_figure( pr );
  }

  // Points
  for ( int i = 0; i < n_points; ++i ) {
    fread( labels, sizeof(int), 2, fp );
    Rectangle* pr = new Rectangle;
    pr->x1 = labels[0]; pr->y1 = labels[1];
    pr->x2 = pr->x1;    pr->y2 = pr->y1;
    insert_an_index_figure( pr );
  }
  
  fclose( fp );

  sort_segment_tree( segRoot );

  delete[] treePrev;

  show_res_use( "Input" );
  return true;
}

void profile_pestrie()
{
  int i, j;
  int ans = 0;

  memset( vis, 0, sizeof(int) * m );
  
  for ( i = 0; i < n; ++i )
    vis[ tree[i] ]++;

  // In permille
  int scales[] = { 1, 5, 10, 50, 100, 200, 300, 400, 500};
  for ( i =0, j = 0; j < 9; ++j ) {
    int firstN = scales[j] * m / 1000;
    for ( ; i < firstN; ++i ) ans += vis[i];
    fprintf( stderr, "Pointed-to matrix (first %d\%-\%) = %.3lf\%\n", scales[j], (double)ans / n * 100 );
  }

  show_res_use( NULL );
}

bool IsAlias( int x, int y )
{
  int tr1, tr2;
  struct segTreeNode *p;
  struct Rectangle *r, *tt;

  tr1 = tree[x];
  if ( tr1 == -1 ) return false;
  tr2 = tree[y];
  if ( tr2 == -1 ) return false;
  if ( tr1 == tr2 ) {
    ++cnt_grp;
    return true;
  }

  x = preV[x]; y = preV[y];
  Rectangle *res = find_treap<Rectangle>( roots[x], y );
  //return (res == NULL ? false : res->y2 >= y);
  if ( res != NULL ) return res->y2 >= y;

  if ( tr1 < tr2 ) { tr1 = tr2; y = x; } 
  p = segRoot;
  while ( p ) {
    // We first check if x and p->mid belongs to the same PES tree
    // This small trick helps quickly walk down the tree
    
    if ( tr1 != p->tr_mid ) {
      p = ( x < p->mid ? p->left : p->right );
    }
    else {
      kvec_t(Rectangle*)& rects = p->sortByY1;
      int mid, s = 0, e = rects.size();
      r = NULL;
      
      if ( e < 5 ) {
	// Directly go over the list
	while ( s < e ) {
	  tt = rects[s];
	  if ( tt -> y2 >= y ) {
	    if ( tt -> y1 <= y ) r = tt;
	    break;
	  }
	  ++s;
	}
      }
      else {
	// Binary search for the enclosing figure
	while ( e - s > 1 ) {
	  mid = (s+e) / 2;
	  tt = rects[mid];
	  if ( tt->y1 <= y ) {
	    if ( tt->y2 >= y ) {
	      // Found the closest one
	      r = tt;
	      break;
	    }
	    s = mid;
	  }
	  else
	    e = mid;
	}
      }
      
      if ( x >= p->mid ) {
	if ( r != NULL ) {
	  if ( x == p->mid ) {
	    //printf ("%d Hit : [%d, %d]\n", x, r->x1, r->x2 );
	    return true;
	  }
	  if ( r->x2 >= x ) {
	    //printf ("%d Hit : [%d, %d]\n", x, r->x1, r->x2 );
	    return true;
	  }
	}
	if ( x == p->mid ) break;
	p = p -> right;
      }
      else {
	if ( r != NULL && r->x1 <= x ) {
	  //printf ("%d Hit : [%d, %d]\n", x, r->x1, r->x2 );
	  return true;
	}
	p = p -> left;
      }
    }
  }
   
  return false;
}

// List query in real use should be passed in a handler.
// That handler decide what to do with the query answer.
// We implement a placeholder statement at the space a qualified result is found.
// This is to avoid the compiler smartly deleting the code.
int ListPointsTo( int x )
{
  int X, tr;
  segTreeNode *p;
  Rectangle *r;
  int ans = 0;

  X = preV[x];
  if ( X == -1 ) return 0;
  p = segRoot;
  tr = tree[x];
  
  printf( "pointer %d (prev = %d) points-to : %d", x, X, tree[x] );
  ++ans;       // The root of this pointer is also a points-to answer

  while ( p ) {
    int size;
    
    if ( tr != p->tr_mid ) {
      p = ( x < p->mid ? p->left : p->right );
    }
    else {
      if ( X >= p-> mid ) {
	kvec_t(Rectangle*) &rects = p->sortByX2Object; 
	size = rects.size(); 
	
	for ( int i = 0; i < size; ++i ) {
	  r = rects[i];
	  if ( r->x2 < X ) break;
	  // r->y1 is an answer
	  printf( " %d", backPrev[r->y1] );
	  ans++;
	}
	
	if ( X == p->mid ) break;
	p = p -> right;
      }
      else {
	kvec_t(Rectangle*) &rects = p->sortByX1Object; 
	size = rects.size(); 
	
	for ( int i = 0; i < size; ++i ) {
	  r = rects[i];
	  if ( r->x1 > X ) break;
	  // r->y1 is an answer
	  printf( " %d", backPrev[r->y1] );
	  ans++;
	}
	
	p = p -> left;
      }
    }
  }

  printf( "\n" );
  return ans;
}

// Find all the pointers y that *x and *y is an alias pair
int ListAliases( int x )
{
  int X;
  segTreeNode *p;
  Rectangle *r;
  int size;
  int ans = 0;

  // We first extract the ES nodes belong to the same subtree
  int tr = tree[x];
  if ( tr == -1 ) return 0;
  for ( int i = tree_members[X]; i > 0; --i )
    ++ans;

  X = preV[x];
  p = segRoot;
  while ( p ) {
    if ( tr != p->tr_mid ) {
      p = ( x < p->mid ? p->left : p->right );
    }
    else if ( X >= p-> mid ) {
      kvec_t(Rectangle*) &rects = p->sortByX2; 
      size = rects.size(); 

      for ( int i = 0; i < size; ++i ) {
	r = rects[i];
	if ( r->x2 < X ) break;
	for ( int j = r->y1; j <= r->y2; ++j ) {
	  // j is an answer	 
	  ++ans;
	}
      }
      
      if ( X == p->mid ) break;
      p = p -> right;
    }
    else {
      kvec_t(Rectangle*) &rects = p->sortByX1; 
      size = rects.size(); 

      for ( int i = 0; i < size; ++i ) {
	r = rects[i];
	if ( r->x1 > X ) break;
	for ( int j = r->y1; j <= r->y2; ++j )
	  // j is an answer
	  ++ans;
      }
      
      p = p -> left;
    }
  }

  return ans;
}

int ListPointedTo( int o )
{
  return ListAliases( o );
}

// A binary research procedure for answering the bulk ListAliases Query
// At each segment tree node, we binary search for the two pointers tha are closeset to the middle line
// then we only use those two pointers to search for the intersected rectangles
// It would be significantly faster than the naive implementation that uses up all the pointers for searching.
static void recursive_enumerate_alias_or( segTreeNode* p, int s, int e, kvec_t(int)& query_points )
{
  Rectangle* r;
  int pi, pos;
  Rectangle **it, **ie;

  if ( p == NULL || s >= e )
    return;

  it = p->sortByX1.begin();
  ie = p->sortByX1.end();

  if ( e - s > 22 ) {
    // We use binary search
    int begin = s, end = e;
    while ( end - begin > 1 ) {
      pi = (begin+end) / 2;
      if (query_points[pi] <= p->mid)
	begin = pi;
      else
	end = pi - 1;
    }
    
    if ( query_points[begin] <= p->mid )
      pi = begin + 1;
    else
      pi = begin;
  }
  else {
    // Linear scan may be faster
    for ( pi = s; pi < e && query_points[pi] <= p->mid; ++pi );
  }

  if ( pi > s ) {
    // At least one of the querying point is to the left of the middle line
    pos = query_points[pi-1];
    while ( it != ie ) {
      r = *it;
      if ( r->x1 > pos ) break;
      internal_container.push_back(r);
      r->selected = true;
      it++;
    }

    recursive_enumerate_alias_or( p->left, s, pi, query_points );

    // Special case
    if ( pos == p->mid ) {
      recursive_enumerate_alias_or( p->right, pi+1, e, query_points );
      return;
    }
  }
  
  if ( pi < e ) {
    // Now we should go over the sorted list to see if we miss something
    it = p->sortByX2.begin();
    ie = p->sortByX2.end();
    pos = query_points[pi];
    while ( it != ie ) {
      r = *it;
      if ( r->x2 < pos ) break;
      if ( r->selected == false ) {
	internal_container.push_back(r);
	r->selected = true;
      }
      it++;
    }
  
    recursive_enumerate_alias_or( p->right, pi, e, query_points );
  }
}

int ListAliasesBulk( kvec_t(int)& query_points )
{
  int ans;

  sort( query_points.begin(), query_points.end() );
  internal_container.clear();
  recursive_enumerate_alias_or( segRoot, 0, query_points.size(), query_points );

  // Decode the answers
  ans = 0;
  memset( vis, 0, sizeof(int) * vertex_num ); 
  
  for ( int i = 0; i < internal_container.size(); ++i ) {
    Rectangle *r = internal_container[i];
    r->selected = false;
    for ( int j = r->y1; j <= r->y2; ++j )
      if ( vis[j] == 0 ) {
	++ans;
	vis[j] = 1;
      }
  }

  return ans;
}

void simulate_queries()
{
  int x, y;
  kvec_t(int) queries;
  FILE *fp;

  kv_init(int, queries);

  if ( query_plan_file != NULL ) {
    fp = fopen( query_plan_file, "r" );
    if ( fp == NULL ) {
      fprintf( stderr, "Cannot open the query plan file. Simulation exits.\n" );
      return;
    }

    // Read
    while ( fscanf( fp, "%d", &x ) != EOF )
      kv_push( int, queries, x );

    fclose( fp );
    n_query = kv_size( queries );
    
    show_res_use( "Query plan loaded" );
  }

  // IsAlias query is treated specially  
  if ( query_type == IS_ALIAS_QUERY ) {
    /*
    n_query = 10;
    ListPointsTo( 21 );
    ListPointsTo( 74 );
    ListPointsTo( 75 );
    */
    //n_query = 10;
    if ( query_plan_file != NULL ) {
      for ( int i = 0; i < n_query; ++i ) {
	x = queries[i];
	for ( int j = i + 1; j < n_query; ++j ) {
	  y = queries[j];
	  bool res = IsAlias( x, y );
	  //printf( "(%d, %d) : %s\n", x, y, res == true ? "true" : "false" );
	}
      }
    }
    else {	  
      for ( int i = 0; i < n_query; ++i ) {
	x = rand() % n;
	y = rand() % n;
	IsAlias( x, y );
      }
    }
  }
  else if ( query_type == BULK_LIST_ALIASES_QUERY ) {
    // Special case again
    if ( query_plan_file == NULL ) {
      for ( int i = 0; i < n_query; ++i ) {
	int x = rand() % n;
	kv_push( int, queries, x );
      }
    }
    
    ListAliasesBulk( queries );
  }
  else {
    int res;

    for ( int i = 0; i < n_query; ++i ) {
      if ( query_plan_file != NULL )
	x = kv_A( queries, i );
      else {
	if ( query_type != LIST_POINTED_TO_QUERY )
	  x = rand() % n;
	else
	  x = n + rand() % m;
      }
      
      switch (query_type ) {
      case LIST_POINTS_TO_QUERY:
	res = ListPointsTo( x );
	break;
	
      case LIST_ALIASES_QUERY:
	res = ListAliases( x );
	break;
	
      case LIST_POINTED_TO_QUERY:
	res = ListPointedTo( x );
	break;
	
      default:
	break;
      }

      //printf( "%d\n", res );
    }
  }

  fprintf ( stderr, "Direct answers = %d, total = %d\n", cnt_grp, n_query * (n_query-1) / 2 );
  show_res_use( "Querying" );
}

int main( int argc, char** argv )
{
  if ( parse_options_querier( argc, argv ) == 0 )
    return -1;

  if ( read_index() == false )
    return -1;

  if ( do_profile ) profile_pestrie();

  simulate_queries();

  return 0;
}
