/*
 * The segment tree based structure for storing figures.
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

  SegTreeNode()
  {
    rects = NULL;
  }
  
  // Delete everything recursively
  ~SegTreeNode()
  {
    clean_treap(rects);
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

//
SegTree* 
build_segtree( int, int);

bool
query_point( SegTree*, int, int );

void
insert_segtree_wrapper( SegTree*, const Rectangle& );

int
dump_figures( SegTree*, FILE* fp );

void
flush_left_shapes( SegTree* );

#endif
