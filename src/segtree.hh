// Copyright 2014, Hong Kong University of Science and Technology. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

/*
 * The segment tree based structure for storing figures.
 *
 * by Xiao Xiao
 * initial: 2012.7
 * improve: 2012.9
 */

#ifndef SEGTREE_H
#define SEGTREE_H

#include <cstdio>
#include "treap.hh"

struct SegTreeNode
{
  // A balanced tree structure that sorts the figures
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

public:
  bool collect_figures( VECTOR(VLine*)& );
};

// The segment tree and its statistical information
struct SegTree
{
  int maxN;
  SegTreeNode **unitNodes;
  int n_points, n_horizs, n_vertis, n_rects;
  int n_out_points, n_out_horizs, n_out_vertis, n_out_rects;
  int n_pairs;

  SegTree( int );
  ~SegTree();

public:
  bool query_point( int, int );
  void insert_segtree( const Rectangle& );
  void flush_left_shapes();
  int dump_figures( std::FILE* );

private:
  void insert_rectangle( int, int, VLine* );
  SegTreeNode* get_unit_node(int);
  void insert_unit_node(int, VLine*);
};

// Construct a segment tree instance
SegTree* 
build_segtree( int, int );


#endif
