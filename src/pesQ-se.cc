/*
 * Processing the PesTrie index to serve the side-effect queries.
 * by richardxx, 2012.7
 */

#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include "options.h"
#include "profile_helper.h"
#include "kvec.hh"

using namespace std;

struct Rectangle
{
  int x1, y1, x2, y2;
  bool selected;

public:
  Rectangle()
  {
    selected = false;
  }

  Rectangle( int X1, int Y1, int X2, int Y2 )
  {
    x1 = X1;
    y1 = Y1;
    x2 = X2;
    y2 = Y2;
    selected = false;
  }
};

struct segTreeNode
{
  // Start and end positions of this segment
  int s, e;
  // The middle position and its tree code of this segment
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

// The number of pointers and objects
int n, m;
int vertex_num, n_rects, n_verticals, n_horizontals, n_points;
int *tree, *preV, *treePrev, *modref_flag;
int *vis;          // used for bulk query

// Containers
int *tree_members;
vector<Rectangle*> internal_container;

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
  struct segTreeNode *p = segRoot;

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

  // We first sort by X1 and pick out all the alpha-rectangles
  sort( p->sortByY1.begin(), p->sortByY1.end(), firstY1Small );

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

  // We first binary search to see if this rectangle covers a PES root
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
  
  insert_segment_tree( pr, preV[s] == pr->y1 );
  
  // Then, we insert an inverse rectangle, and it cannot represent the points-to information whatever
  Rectangle* pr_prime = new Rectangle(pr->y1, pr->x1, pr->y2, pr->x2);
  insert_segment_tree( pr_prime, false );
}

// Read the comments in the pes-se-indexer.cc for index file format details
bool read_index()
{
  int half_m;
  char magic_code[8];
  FILE *fp;

  fp = fopen( input_file, "rb" );
  if ( fp == NULL )
    return false;

  fread( magic_code, sizeof(char), 4, fp );
  magic_code[4] = 0;
  if ( strcmp( magic_code, "SEP1" ) != 0 ) {
    fprintf( stderr, "This file is not a valid PesTrie index for side-effect information.\n" );
    return false;
  }
  
  fread( &n, sizeof(int), 1, fp );
  fread( &m, sizeof(int), 1, fp );
  fread( &n_rects, sizeof(int), 1, fp );
  fread( &n_verticals, sizeof(int), 1, fp );
  fread( &n_horizontals, sizeof(int), 1, fp );
  fread( &n_points, sizeof(int), 1, fp );

  // Allocate resource
  tree = new int[n+m];
  preV = new int[n+m];
  treePrev = new int[n+m];
  modref_flag = new int[n];
  vis = new int[n+m];
  tree_members = new int[m];
  
  // Read in the pre-order descriptors
  fread( preV, sizeof(int), n+m, fp );
  
  // We re-discover the tree codes through the binary search
  // Also we recover the number of tree verteices in the PesTrie index
  // And also, we discover the type (mod/ref) of each statement
  vertex_num = -1;
  half_m = m / 2;
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

      // Tree code
      tree[i] = s - n;
      treePrev[pes_node] = s - n;
      // mod/ref flag
      if ( (s-n) < half_m )
	modref_flag[i] = MOD_SET;
      else
	modref_flag[i] = REF_SET;
      // Count how many tree nodes for each tree
      if ( vis[pes_node] == 0 ) tree_members[s-n]++;
      ++vis[pes_node];
    }
    else
      tree[i] = -1;
    
    if ( pes_node > vertex_num ) 
      vertex_num = preV[i];
  }

  // The tree codes for the objects
  for ( int i = 0; i < m; ++i ) tree[ i + n ] = i;

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

  // We sort the indexed rectangles for future search
  sort_segment_tree( segRoot );

  // Delete temporary resource
  delete[] treePrev;

  show_res_use( "Input" );
  return true;
}

// List query in real use should be passed in a handler.
// That handler decide what to do with the query answer.
// We implement a placeholder statement at the space a qualified result is found.
// This is to avoid the compiler smartly deleting the code.
int ListLoadsStores( int x )
{
  int tr, X;
  segTreeNode *p;
  Rectangle *r;
  int size;
  int ans = 0;

  // We first extract the pointers belong to the same subtree
  tr = tree[x];
  if ( tr == -1 ) return 0;
  if ( modref_flag[x] == MOD_SET ) {
    for ( int i = tree_members[tr]; i > 0; --i )
      ++ans;
  }
  
  X = preV[x];
  p = segRoot;
  while ( p ) {
    if ( tr != p->tr_mid ) {
      p = ( x < p->mid ? p->left : p->right );
    }
    else {
      if ( X >= p-> mid ) {
	kvec_t(Rectangle*) &rects = p->sortByX2; 
	size = rects.size(); 
	
	for ( int i = 0; i < size; ++i ) {
	  r = rects[i];
	  if ( r->x2 < X ) break;
	  for ( int j = r->y1; j <= r->y2; ++j )
	    // j is an answer
	    ++ans;
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
  else {
    n_query = n;
    //printf( "%d\n", n_query );
  }

  // IsAlias query is treated specially  
  for ( int i = 0; i < n_query; ++i ) {
    if ( query_plan_file != NULL )
      x = queries[i];
    else {
      x = i;
    }
    
    if ( query_type == LIST_LOADS_QUERY
	 && modref_flag[i] != MOD_SET )
      continue;
    if ( query_type == LIST_STORES_QUERY
	 && modref_flag[i] != REF_SET )
      continue;

    int ans = ListLoadsStores(x);
    //printf( "%d\n", ans );
  }

  show_res_use( "Querying" );
}

int main( int argc, char** argv )
{
  if ( parse_options_querier( argc, argv ) == 0 )
    return -1;

  if ( read_index() == false )
    return -1;

  simulate_queries();

  return 0;
}
