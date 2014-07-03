/*
 * Processing the PesTrie index to serve queries.
 *
 * Initiated 2011.8
 * Improved 2012.7: improve the alias query to O(lgn) practical performance.
 * Improved 2012.10: fix bugs and code refactoring.
 * Improved 2014.02: fix bugs and improve performance.
 * Improved 2014.05: significantly improve the decoding performance and querying memory usage
 */

#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include <ctime>
#include <set>
#include "shapes.hh"
#include "profile_helper.h"
#include "histogram.hh"
//#include "bitmap.h"
#include "query.hh"
#include "query-inl.hh"

using namespace std;

// The segment tree node
// We still use segment tree as the fundamental querying structure
class SegNode
{
public:
  // This node represents the range [l, r]
  int l, r;
  SegNode *left, *right, *parent;
  bool merged;

  VECTOR(VLine*) rects;
  VECTOR(VLine*) vertis;
  //VECTOR(int) pointsto;

  SegNode()
  {
    left = right = NULL;
    parent = NULL;
    merged = false;
  }

  void add_rect( VLine* p ) 
  { 
    rects.push_back(p);
  }

  void add_vertis( VLine* p ) 
  { 
    vertis.push_back(p);
  }

  int n_of_rects() 
  {
    return rects.size();
  }
};

// The segment tree structure for querying system
class SegTree
{
public:
  SegTree(int maxN)
  {
    segUnits = new SegNode[maxN];
    segRoot = build_seg_tree(0, maxN);
  }

  void optimize_seg_tree( SegNode* p );
  void recursive_merge( SegNode* p );

private:
  SegNode* build_seg_tree( int l, int r );

private:
  // The querying segment tree
  SegNode **segUnits, *segRoot;
};


SegNode*
SegTree::build_seg_tree( int l, int r )
{
  int x = (l+r) / 2;
  SegNode* p = segUnits[x]; 
  
  p->l = l;
  p->r = r;

  if ( l < x ) {
    p->left = build_seg_tree( l, x - 1 );
    p->left->parent = p;
  }
  
  if ( x + 1 < r ) {
    p->right = build_seg_tree( x + 1, r ); 
    p->right->parent = p;
  }
  
  return p;
}

/*
 * We update the parent links to skip the empty nodes.
 */ 
void
PesQS::optimize_seg_tree( SegNode* p ) 
{
  SegNode* q = p->parent;
  
  if ( q != NULL &&
       q->n_of_rects() == 0 ) {
    q = q->parent;
    p->parent = q;
  }
  
  if ( q == NULL )
    p->merged = true;

  if ( p->left != NULL )
    optimize_seg_tree( p->left );
  if ( p->right != NULL )
    optimize_seg_tree( p->right );
}

void
SegTree::insert_point( int x, int y )
{
  VLine* p = new VLine(y, y);
  segUnits[x]->add_vertis(p);
}

/*
 * We add at most log(n) references to each figure.
 * [x1, x2]: the X range of insertion rectangle
 */
void
SegTree::insert_rect( int x1, int x2, VLine* pr, SegNode* p )
{
  if ( x1 <= p->x1 && x2 >= p->x2 ) {
    // We only consider the full coverage
    p->add_rect(pr);
    return;
  }

  int x = (p->x1 + p->x2) / 2;
  p = segUnits[x];
  
  if ( x1 < x ) insert_rect( x1, x2, pr, p->left );
  if ( x2 > x ) insert_rect( x1, x2, pr, p->right );
}



class PesQS : public IQuery
{
public:
  // Interface functions
  bool IsAlias( int x, int y);
  int ListPointsTo( int x, const IFilter* filter );
  int ListAliases( int x, const IFilter* filter );
  int ListPointedBy( int o, const IFilter* filter );
  int ListModRefVars( int x, const IFilter* filter );
  int ListConflicts( int x, const IFilter* filter );

public:
  int getPtrEqID(int x) { return preV[x]; }
  int getObjEqID(int x) { return preV[x+n]; }

public:
  int nOfPtrs() { return n; }
  int nOfObjs() { return m; }
  int getIndexType() { return index_type; }

public:
  PesQS(int n_ptrs, int n_objs, int n_vertex, int type, int d_merging)
  {
    n = n_ptrs; m = n_objs; vertex_num = n_vertex;
    index_type = type;
    demand_merging = d_merging;

    tree = new int[n_ptrs+n_objs];
    preV = new int[n_ptrs+n_objs];
    root_prevs = new int[n_objs+1];
    root_tree = new int[n_vertex+1];
    es2ptrs = new VECTOR(int)[n_vertex];
    
    if ( type == PT_MATRIX ) {
      es2objs = new VECTOR(int)[n_vertex];
    }

    qtree = new SegTree(n_vertex);
  }

public:
  void rebuild_eq_groups( FILE* fp );
  
private:
  SegTree* qtree;

  // The maximum preorder timestamp for the store statements
  int max_store_prev;

  // The number of pointers, objects, trees, and nodes
  int n, m, n_trees, vertex_num;
  // Mapping from pointer to tree ID
  int *tree;
  // Mapping from pointer to PesTrie ID
  int *preV;
  // Recording pre-order of the roots
  int *root_prevs;
  // Mapping the pre-order of root node to tree ID (-1 means this node is not a root)
  int *root_tree;

  // Mapping from equivalent set to pointers/objects
  VECTOR(int) *es2ptrs, *es2objs;

  // Points-to or side-effect information
  int index_type;
  
  // Merging the aliasing information bottom up on demand
  bool demand_merging;
};


//Statistics
//static int cnt_same_tree = 0;


void
PesQS::rebuild_eq_groups( FILE* fp )
{
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
    pesqs->max_store_prev = root_prevs[m/2];

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
      es2ptrs[preI].push_back(i);
    }
    else
      tree[i] = -1;
  }
}


static bool 
comp_rect( VLine* r1, VLine* r2 )
{
  return r1->y1 < r2->y1;
}

void
PesQS::load_figures( FILE* fp )
{
  // First we initialize the segment tree
  segRoot = build_seg_tree( 0, vertex_num );

  int n_labels;
  int *labels = new int[vertex_num * 3];
  VECTOR(Rectangle*) all_rects(vertex_num);

  // Loading and processing each figure
  for ( int x1 = 0; x1 < vertex_num; ++x1 ) {
    fread( &n_labels, sizeof(int), 1, fp );
    if ( n_labels == 0 ) continue;
    fread( labels, sizeof(int), n_labels, fp );

    int i = 0;
    while ( i < n_labels ) {
      int y1 = labels[i++];
      int x2, y2;

      if ( (y1&SIG_FIGURE) == SIG_POINT ) {
	// This is a point, directly insert it
	insert_point( x1, y1 );
	insert_point( y1, x1 );
	continue;
      }

      if ( (y1&SIG_FIGURE) == SIG_VERTICAL ) {
	y1 &= ~SIG_VERTICAL;
	y2 = labels[i++];
	x2 = x1;
      }
      else if ( (y1&SIG_FIGURE) == SIG_HORIZONTAL ) {
	y1 &= ~SIG_HORIZONTAL;
	x2 = labels[i++];
	y2 = y1;
      }
      else {
	y1 &= ~SIG_RECT;
	x2 = labels[i++];
	y2 = labels[i++];
      }

      // First, insert the reversed figure directly
      VLine* p = new VLine(x1, x2);
      if ( y1 == y2 ) {
	segUnits[y1]->add_vertis(p);
      }
      else {
	insert_rect( y1, y2, p, segRoot );
      }

      // Second, we cache the original figure
      Rectangle* r = new Rectangle(x1, x2, y1, y2);
      if ( x1 == x2 ) {
	segUnits[x1]->add_vertis(r);
      }
      else {
	all_rects.push_back(r);
      }
    }
  }

  // We sort the cached figures by y1
  sort( all_rects.begin(), all_rects.end(), comp_rect );

  // We insert the cached figures
  int size = all_rects.size();
  for ( int i = 0; i < size; ++i ) {
    Rectangle* r = all_rects[i];
    insert_rect( r->x1, r->x2, r, segRoot );
  }

  // Optimize the parent links
  optimize_seg_tree( segRoot );
  delete[] labels;
}

void 
PesQS::profile_pestrie()
{
  int non_empty_nodes = 0;
  int n_pts = 0, n_rects = 0;

  for ( int i = 0; i < vertex_num; ++i ) {
    // Test 1:
    VECTOR(int) ptrs = es2ptrs[i];
    if ( ptrs.size() > 0 )
      ++non_empty_nodes;

    // Test 2:
    struct SegNode* p = segUnits[i];
    n_rects += p->rects.size();
    n_pts += p->vertis.size();
  }
  
  n_pts += non_empty_nodes;
  
  fprintf( stderr, "Trees = %d, ES = %d, Non-empty ES = %d\n", 
	   n_trees, vertex_num, non_empty_nodes );
  fprintf( stderr, "Number of rects = %d\n", 
	   n_rects + n_pts ); 
}

// We do merging sort
void
PesQS::recursive_merge( SegNode* p )
{
  if ( p->merged == true ||
       p->parent == NULL ) return;
  recursive_merge( p->parent );

  VECTOR(VLine*) &list1 = p->parent->rects;
  int sz1 = list1.size();
  int sz2 = p->rects.size();

  if ( sz1 != 0 ) {
    // We have something to merge top down
    VECTOR(VLine*) &list3 = p->rects;
    if ( sz2 == 0 ) {
      // Fast path, just copy
      list3.add_all( list1 );
    }
    else {
      VECTOR(VLine*) list2;
      list2.add_all( list3 );
      list3.clear();

      int i = 0, j = 0;
      VLine *r1 = list1[0], *r2 = list2[0];
      
      while ( r1 != NULL || r2 != NULL ) {
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

  s = 0; 
  e = rects.size();
  /*
  if ( e < 3 ) {
    // fast path
    while ( s < e ) {
      tt = rects[s];
      if ( tt->y2 >= y &&
	   tt->y1 <= y ) return true;
      ++s;
    }
  }
  else {
  */
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
    //  }

  return false;
}


bool
PesQS::IsAlias( int x, int y )
{
  struct SegNode *p;

  int tr_x = tree[x];
  if ( tr_x == -1 ) return false;
  int tr_y = tree[y];
  if ( tr_y == -1 ) return false;
  if ( tr_x == tr_y ) {
    //++cnt_same_tree;
    return true;
  }

  x = preV[x];
  y = preV[y];
  p = segUnits[x];

  // Search the verticals special to this x axis
  if ( binary_search( p->vertis, y ) ) 
    return true;

  if ( !demand_merging ) {
    // We traverse the segment tree bottom up
    while ( p != NULL ) {    
      if ( binary_search( p->rects, y ) ) 
	return true;
      p = p->parent;
    }
  }
  else {
    recursive_merge(p);
    if ( binary_search( p->rects, y ) )
      return true;
  }

  return false;
}

// List query in real use should be passed in a handler.
// That handler decide what to do with the query answer.
int 
PesQS::ListPointsTo( int x, const IFilter* filter )
{
  SegNode *p;
  VECTOR(int) *objs;
  int ans = 0;

  int tr = tree[x];
  if ( tr == -1 ) return 0;
  
  // Don't forget x -> tree[x]
  objs = &es2objs[tr];
  ans += iterate_equivalent_set( objs );
  
  x = preV[x];
  p = segUnits[x];
  
  return ans;
}

// Find all the pointers y that *x and *y is an alias pair
int 
PesQS::ListAliases( int x, const IFilter* filter )
{
  SegNode *p;
  VLine *r;
  int size;

  int tr = tree[x];
  if ( tr == -1 ) return 0;
  
  int ans = 0;
  x = preV[x];

#define VISIT_POINT(v)							\
  ans += iterate_equivalent_set( ptrs, filter );			
  
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

  p = segUnits[x];

  // Visit the verticals
  VECTOR(VLine*) &vertis = p->vertis;
  size = vertis.size();
  for ( int i = 0; i < size; ++i ) {
    VLine* r = vertis[i];
    VISIT_RECT(r);
  }

  // traverse the rectangles up the tree
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

int 
PesQS::ListPointedBy( int o, const IFilter* filter )
{
  return ListAliases( o + n, filter );
}

int 
PesQS::ListModRefVars( int x, const IFilter* filter )
{
  return ListPointsTo( x, filter );
}

int
PesQS::ListConflicts( int x, const IFilter* filter )
{
  int ans = 0;

  if ( preV[x] < max_store_prev )
    ans = ListAliases( x, filter );
  
  return ans;
}

PesQS*
load_pestrie_index(FILE* fp, int index_type, int d_merging )
{
  // Loading the header info
  fread( &n, sizeof(int), 1, fp );
  fread( &m, sizeof(int), 1, fp );
  fread( &vertex_num, sizeof(int), 1, fp );
  
  // Initialize the querying struture
  PesQS* pesqs = new PesQS( n, m, vertex_num, index_type, d_merging );
  
  // Loading and decoding the persistence file
  pesqs->rebuild_eq_groups(fp);
  pesqs->load_figures(fp, pesqs);
  pesqs->profile_pestrie();
  
  fprintf( stderr, "\n-------Input: %s-------\n", input_file );
  show_res_use( "Index loading" );
  
  return pesqs;
}
