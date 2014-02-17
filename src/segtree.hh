/*
 * The segment tree tailored for the PesTrie algorithm.
 * By richardxx, 2010.6
 * Updated, 2012.9
 */

#ifndef SEGTREE_H
#define SEGTREE_H

#include "treap.hh"

struct SegTreeNode
{
  // The size of segment tree would be larger if we split the rectangles and points
  struct TreapNode *rects;
  struct TreapNode *points;
  //bool left, right;

  SegTreeNode()
  {
    rects = NULL;
    points = NULL;
    //left = right = false;
  }
  
  // Delete everything recursively
  ~SegTreeNode()
  {
    if ( rects != NULL ) delete rects;
    if ( points != NULL ) delete points;
  }
};

// The segment tree and its statistical information
struct SegTree
{
  int maxN;
  SegTreeNode **unitRoots;
  int n_points, n_horizs, n_vertis, n_rects;
  int n_pairs;

  SegTree( int );
  ~SegTree();
};

// To help implement a general visitor
typedef void (*SEG_TREE_VISITOR)( SegTreeNode*, void* );

extern SegTree* build_segtree( int, int);
extern bool query_point( SegTree*, int, int );
extern void insert_segtree_wrapper( SegTree*, const Rectangle& );
extern void visit_segTree( SegTree*, SEG_TREE_VISITOR, void* );

#endif
