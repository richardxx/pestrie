// Copyright 2014, Hong Kong University of Science and Technology. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

/*
 * Processing the PesTrie index to serve queries.
 *
 * by Xiao Xiao
 * Initial: 2011.8
 * Improved 2012.7: improve the alias query to O(lgn) practical performance.
 * Improved 2012.10: fix bugs and code refactoring.
 * Improved 2014.02: fix bugs and improve performance.
 * Improved 2014.05: significantly improve the decoding performance and querying memory usage
 */

#include <cstdio>
#include <cstring>
#include <algorithm>
#include <ctime>
#include <cassert>
#include "options.hh"
#include "shapes.hh"
#include "query.hh"
#include "query-inl.hh"
#include "profile_helper.h"

using namespace std;

static bool 
comp_strips( VLine* r1, VLine* r2 )
{
  return r1->y1 < r2->y1;
}

static bool 
comp_rect( Rectangle* r1, Rectangle* r2 )
{
  return r1->y1 < r2->y1;
}

static void
push_and_merge( VECTOR(VLine*) &rects, VLine* p )
{
  if ( rects.size() != 0 ) {
    // Return a reference to the last entry
    // Update to vlast will be reflected to the last entry
    VLine* &vlast = rects.back();
    if ( vlast->y2 + 1 == p->y1 ) {
      // Concatenate two rectangles
      int y1 = vlast->y1;
      /* 
	 Warning:
	 Could lead to a memory leak if vlast is the only pointer to the rectangle
	 Can be solved by a manual reference counting but we didn't do it here
      */
      vlast = new VLine(y1, p->y2);
      return;
    }
  }
  
  rects.push_back(p);
}

// The segment tree node
// We still use segment tree as the fundamental querying structure
class SegNode
{
public:
  // This node represents the range [l, r]
  int l, r;
  SegNode *left, *right, *parent;
  bool merged, pt_extracted;
  
  VECTOR(VLine*) rects;
  VECTOR(int) pointsto;
  
  SegNode()
  {
    left = right = NULL;
    parent = NULL;
    merged = false;
    pt_extracted = false;
  }
  
  void add_rect( VLine* p ) 
  { 
    push_and_merge( rects, p );
  }
  
  void add_points_to( int o )
  {
    pointsto.push_back(o);
  }

  int n_of_rects() { return rects.size(); }
};

class SegUnitNode : public SegNode
{
public:
  VECTOR(VLine*) strips;

  SegUnitNode() {}

  void add_strip( VLine* p )
  {
    strips.push_back(p);
  }

  int n_of_strips() { return strips.size(); }
};


// The segment tree structure for querying system
class SegTree
{
public:
  SegTree(int n_range)
  {
    unitNodes = new SegUnitNode*[n_range];
    segRoot = build_seg_tree(0, n_range-1);
    maxN = n_range;
  }

  ~SegTree()
  {
    free(unitNodes);
    free_seg_tree( segRoot );
  }

  SegNode* get_unit_node( int x ) 
  { 
    return unitNodes[x]; 
  }

  void optimize_seg_tree();
  void recursive_merge( SegNode* p );
  void insert_point( int x, VLine* p );
  void insert_rect( int x1, int x2, VLine* pr );
  bool verify();

private:
  SegNode* build_seg_tree( int l, int r );
  void free_seg_tree( SegNode* );
  void __insert_rect( int x1, int x2, VLine* pr, SegNode* p );
  void __opt_seg_tree( SegNode* p );
  
private:
  // The pointers to the unit nodes
  SegUnitNode **unitNodes;
  // The pointer to the root
  SegNode *segRoot;
  // The range of X-aixs
  int maxN;
};


SegNode*
SegTree::build_seg_tree( int l, int r )
{
  SegNode* p;

  if ( l == r ) {
    p = new SegUnitNode;
    unitNodes[l] = (SegUnitNode*)p;
  }
  else {
    p = new SegNode;
    int x = (l+r) / 2;
    if ( l <= x ) {
      p->left = build_seg_tree( l, x );
      p->left->parent = p;
    }
    
    if ( x < r ) {
      p->right = build_seg_tree( x + 1, r ); 
      p->right->parent = p;
    }
  }
  
  p->l = l;
  p->r = r;
  return p;
}

void
SegTree::free_seg_tree( SegNode* p )
{
  if ( p->left != NULL )
    free_seg_tree(p->left);
  
  if ( p->right != NULL )
    free_seg_tree(p->right);

  free(p);
}

/*
 * We update the parent links to skip the empty nodes.
 * The points-to information is also extracted.
 */ 
void
SegTree::__opt_seg_tree( SegNode* p ) 
{
  SegNode* q = p->parent;
  
  if ( q != NULL &&
       q->n_of_rects() == 0 ) {
    q = q->parent;
    p->parent = q;
  }

  if ( p->left != NULL )
    __opt_seg_tree( p->left );
  if ( p->right != NULL )
    __opt_seg_tree( p->right );
}

// Merge sort, the result is given in list1
static void
merge_into( VECTOR(VLine*) &list1, VECTOR(VLine*) &list2 )
{
  int sz1 = list1.size();
  int sz2 = list2.size();

  if ( sz2 == 0 ) return;

  if ( sz1 == 0 ) {
    // Fast path, just copy
    list1.add_all( list2 );
    return;
  }
  
  VECTOR(VLine*) list3(sz1+sz2);
  int i = 0, j = 0;
  VLine *r1 = list1[0], *r2 = list2[0];
  
  while ( i < sz1 || j < sz2 ) {
    if ( j == sz2 ||
	 ( i < sz1 && r1->y1 < r2->y1 ) ) {
      push_and_merge( list3, r1 );
      ++i;
      if ( i < sz1 ) r1 = list1[i];
    }
    else {
      push_and_merge( list3, r2 );
      ++j;
      if ( j < sz2 ) r2 = list2[j];
    }
  }

  list1.swap(list3);
}

void
SegTree::optimize_seg_tree()
{
  // We first merge the figures in the unit nodes
  for ( int i = 0; i < maxN; ++i ) {
    SegUnitNode* p = unitNodes[i];
    if ( p->n_of_strips() != 0 ) {
      merge_into( p->rects, p->strips );
    }
  }

  // We recursively process the figures
  __opt_seg_tree( segRoot );
}

// We do merging sort
void
SegTree::recursive_merge( SegNode* p )
{
  if ( p->merged == true ||
       p->parent == NULL ) return;
  
  // process its parent first
  recursive_merge( p->parent );

  VECTOR(VLine*) &list1 = p->rects;
  VECTOR(VLine*) &list2 = p->parent->rects;
  merge_into( list1, list2 );
  p->merged = true;
}

void
SegTree::insert_point( int x, VLine* p )
{
  unitNodes[x]->add_strip(p);
}

/*
 * We add at most log(n) references to each figure.
 * [x1, x2]: the X range of the rectangle for insersion
 */
void
SegTree::__insert_rect( int x1, int x2, VLine* pr, SegNode* p )
{
  if ( x1 <= p->l && x2 >= p->r ) {
    // We only consider the full coverage
    p->add_rect(pr);
    return;
  }

  int x = (p->l + p->r) / 2;  
  if ( x1 <= x ) __insert_rect( x1, x2, pr, p->left );
  if ( x2 > x ) __insert_rect( x1, x2, pr, p->right );
}

void
SegTree::insert_rect( int x1, int x2, VLine* pr )
{
  __insert_rect( x1, x2, pr, segRoot );
}


bool
SegTree::verify()
{
  for ( int i = 0; i < maxN; ++i ) {
    SegNode* p = unitNodes[i];
    while ( p != NULL && 
	    !p->merged ) {
      VECTOR(VLine*) &rects = p->rects;
      int size = rects.size();
      int lasty = -1;
      for ( int j = 0; j < size; ++j ) {
	VLine* r = rects[j];
	if ( r->y1 < lasty ) {
	  fprintf( stderr, "Error [%d, %d]! X = %d, y1 = %d, lasty = %d, index = %d\n",
		   p->l, p->r, i, r->y1, lasty, j );
	  assert(r->y1 >= lasty);
	}
	lasty = r->y2;
      }

      // Temporarily use merged as a stop indicator
      p -> merged = true;
      p = p->parent;
    }
  }

  for ( int i = 0; i < maxN; ++i ) {
    SegNode* p = unitNodes[i];
    while ( p!= NULL && p->merged ) {
      p->merged = false;
      p = p->parent;
    }
  }

  return true;
}


// The querying interface for Pestrie
class PesQS : public IQuery
{
public:
  // Interface functions
  bool IsAlias( int x, int y);
  int ListPointsTo( int x, IFilter* filter );
  int ListAliases( int x, IFilter* filter );
  int ListPointedBy( int o, IFilter* filter );
  int ListModRefVars( int x, IFilter* filter );
  int ListConflicts( int x, IFilter* filter );

public:
  int getPtrEqID(int x) { return preV[x]; }
  int getObjEqID(int x) { return preV[x+n]; }
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
    root_tree = new int[n_vertex];
    es2ptrs = new VECTOR(int)[n_vertex];
    qtree = new SegTree(n_vertex);
    es2objs = NULL;

    if ( type == PT_MATRIX ) {
      es2objs = new VECTOR(int)[n_vertex];
    }    
  }

  ~PesQS()
  {
    if ( qtree != NULL ) delete qtree;
    if ( tree != NULL ) delete[] tree;
    if ( preV != NULL ) delete[] preV;
    if ( root_prevs != NULL ) delete[] root_prevs;
    if ( root_tree != NULL ) delete[] root_tree;
    if ( es2ptrs != NULL ) delete[] es2ptrs;
    if ( es2objs != NULL ) delete[] es2objs;
  }
  
public:
  void load_figures( FILE* );
  
private:
  void rebuild_mapping_info( FILE* );
  void extract_pointsto(SegNode*);

private:
  SegTree* qtree;
  
  // The maximum preorder timestamp for the store statements
  int max_store_prev;
  // The number of pointers, objects, trees, and nodes
  int n, m, n_trees, vertex_num;
  // Mapping from pointer and object to tree ID
  int *tree;
  // Mapping from pointer and object to pre-order stamp
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


void
PesQS::rebuild_mapping_info( FILE* fp )
{
  // Read in the pre-order descriptors for both pointers and objects
  fread( preV, sizeof(int), n+m, fp );

  // We label the time-stamps that could be roots
  n_trees = 0;
  memset ( root_tree, -1, sizeof(int) * vertex_num );
  
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
      es2objs[v].push_back(i);
      tree[i+n] = tr;
    }
  }

  if ( index_type == SE_MATRIX ) 
    max_store_prev = root_prevs[m/2];

  // We re-discover the tree codes for pointers through the binary search
  // Sentinels
  root_prevs[n_trees] = vertex_num;
  
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

void
PesQS::load_figures( FILE* fp )
{
  int n_points = 0, n_horizs = 0, n_vertis = 0, n_rects = 0;
  int cross_pairs = 0;
  int n_labels;
  int *labels = new int[vertex_num * 3];
  VECTOR(Rectangle*) all_rects(vertex_num);

  // We rebuild the mapping between pointers to Pestrie constructs
  rebuild_mapping_info(fp);

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
	x2 = x1;
	y2 = y1;
	++n_points;
      }
      else if ( (y1&SIG_FIGURE) == SIG_VERTICAL ) {
	y1 &= ~SIG_VERTICAL;
	y2 = labels[i++];
	x2 = x1;
	++n_vertis;
      }
      else if ( (y1&SIG_FIGURE) == SIG_HORIZONTAL ) {
	y1 &= ~SIG_HORIZONTAL;
	x2 = labels[i++];
	y2 = y1;
	++n_horizs;
      }
      else {
	y1 &= ~SIG_RECT;
	x2 = labels[i++];
	y2 = labels[i++];
	++n_rects;
      }

      assert( x1 <= x2 && x2 <= y1 && y1 <= y2 );
      
      // First, insert the reversed figure 
      // (must insert before the original figure)
      VLine* p = new VLine(x1, x2);
      if ( y1 == y2 ) {
	qtree->insert_point(y1, p);
      }
      else {
	qtree->insert_rect(y1, y2, p);
      }

      // Second, we cache the original figure
      if ( x1 == x2 ) {
	p = new VLine(y1, y2);
	qtree->insert_point(x1, p);
      }
      else {
	Rectangle* r = new Rectangle(x1, x2, y1, y2);
	all_rects.push_back(r);
      }
      
      cross_pairs += (x2-x1+1)*(y2-y1+1)*2;
    }
  }

  int size = all_rects.size();

  // We sort the cached figures by y1
  sort( all_rects.begin(), all_rects.end(), comp_rect );

  // We insert the cached figures
  int lasty = -1;
  for ( int i = 0; i < size; ++i ) {
    Rectangle* r = all_rects[i];
    qtree->insert_rect( r->x1, r->x2, r );
  }

  // Optimize the parent links
  qtree->optimize_seg_tree();
  delete[] labels;

#ifdef DEBUG
  qtree->verify();
#endif
  
  // Profile
  int non_empty_nodes = 0;
  int internal_pairs = 0;

  for ( int i = 0; i < vertex_num; ++i ) {
    VECTOR(int) &ptrs = es2ptrs[i];
    if ( ptrs.size() > 0 )
      ++non_empty_nodes;
  }

  for ( int i = 0; i < n_trees; ++i ) {
    int sz = root_prevs[i+1] - root_prevs[i];
    internal_pairs += sz * (sz-1) / 2;
  }
  
  fprintf( stderr, "Trees = %d, ES = %d, Non-empty ES = %d\n", 
	   n_trees, vertex_num, non_empty_nodes );

  fprintf( stderr, 
	   "Points = %d, Verticals = %d, Horizontals = %d, Rectangles = %d\n",
	   n_points, n_vertis, n_horizs, n_rects );

  fprintf( stderr,
	   "Alias pairs = %d\n", internal_pairs + cross_pairs );
}

static bool
binary_search( VECTOR(VLine*) &rects, int y )
{
  int mid, s, e;
  VLine *r, *tt;

  s = 0; 
  e = rects.size();

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

  /*
  // Processing small case
  while ( s < e ) {
    tt = rects[s];
    if ( tt->y2 >= y &&
	 tt->y1 <= y ) return true;
    ++s;
  }
  */

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
  p = qtree->get_unit_node(x);
  if ( p == NULL ) return false;

  if ( !demand_merging ) {
    // We traverse the segment tree bottom up
    do {
      if ( binary_search( p->rects, y ) ) 
	return true;
      p = p->parent;
    } while ( p!= NULL );
  }
  else {
    qtree->recursive_merge(p);
    if ( binary_search( p->rects, y ) )
      return true;
  }

  return false;
}

void
PesQS::extract_pointsto( SegNode* p )
{
  VECTOR(VLine*) &rects = p->rects;
  int size = rects.size();

  for ( int i = 0; i < size; ++i ) {
    VLine* r = rects[i];
    int lower = r->y1;
    int upper = r->y2;
    do {
      if ( root_tree[lower] != -1 )
	p->add_points_to(lower);
      ++lower;
    } while (lower <= upper);
  }

  p->pt_extracted = true;
}

// List query in real use should be passed in a handler.
// That handler decide what to do with the query answer.
int 
PesQS::ListPointsTo( int x, IFilter* filter )
{
  int tr = tree[x];
  if ( tr == -1 ) return 0;
  
  // Don't forget x points-to tree[x]
  int ans = iterate_equivalent_set( es2objs[root_prevs[tr]], filter );
  
  x = preV[x];
  SegNode* p = qtree->get_unit_node(x);

  // traverse the rectangles up the tree
  while ( p != NULL ) {
    if ( p->pt_extracted == false )
      extract_pointsto(p);

    VECTOR(int) &pointsto = p->pointsto;
    int size = pointsto.size();
    for ( int i = 0; i < size; ++i ) {
      int o = pointsto[i];
      ans += iterate_equivalent_set( es2objs[o], filter );	
    }
    p = p->parent;
  }

  return ans;
}

// Find all the pointers y that *x and *y is an alias pair
int 
PesQS::ListAliases( int x, IFilter* filter )
{
  int tr = tree[x];
  if ( tr == -1 ) return 0;
  
  int ans = 0;

  // We first extract the ES groups that belong to the same subtree
  {
    int upper = root_prevs[tr+1];
    for ( int i = root_prevs[tr]; i < upper; ++i ) {
      ans += iterate_equivalent_set( es2ptrs[i], filter );
    }
  }

  x = preV[x];
  SegNode *p = qtree->get_unit_node(x);

  // traverse the rectangles up the tree
  while ( p != NULL ) {
    VECTOR(VLine*) &rects = p->rects;
    int size = rects.size();
    for ( int i = 0; i < size; ++i ) {
      VLine *r = rects[i];
      int lower = r->y1;
      int upper = r->y2;
      do {
	ans += iterate_equivalent_set( es2ptrs[lower], filter );	
	++lower;
      } while ( lower <= upper );			
    }
    p = p->parent;
  }
  
  return ans;
}

int 
PesQS::ListPointedBy( int o, IFilter* filter )
{
  return ListAliases( o + n, filter );
}

int 
PesQS::ListModRefVars( int x, IFilter* filter )
{
  return ListPointsTo( x, filter );
}

int
PesQS::ListConflicts( int x, IFilter* filter )
{
  return ListAliases( x, filter );
}

IQuery*
load_pestrie_index(FILE* fp, int index_type, bool d_merging )
{
  int n, m, vertex_num;

  // Loading the header info
  fread( &n, sizeof(int), 1, fp );
  fread( &m, sizeof(int), 1, fp );
  fread( &vertex_num, sizeof(int), 1, fp );
  
  // Initialize the querying struture
  PesQS* pesqs = new PesQS( n, m, vertex_num, index_type, d_merging );
  fprintf( stderr, "----------Index File Info----------\n" );

  // Loading and decoding the persistence file
  pesqs->load_figures(fp);
  
  return pesqs;
}
