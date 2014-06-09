/*
 * The segment tree tailored for the PesTrie algorithm.
 * By richardxx, 2010.6
 * Updated, 2012.9
 */

#ifndef SEGTREE_H
#define SEGTREE_H

#include "treap.hh"
#include <cstdio>
using std::FILE;


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
