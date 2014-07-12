// Copyright 2014, Hong Kong University of Science and Technology. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

/*
 * A high performance probablistic balanced tree.
 * The implementation is tailored for the PesTrie algorithm.
 *
 * by Xiao Xiao
 * initial: 2012.10
 */

#ifndef TREAP_H
#define TREAP_H

#define INDEX_UTILITY
#include "shapes.hh"
#include "options.hh"

// Finally, the treap node
struct TreapNode
{
  int rkey;                                // A random key used to maintain the heap property
  struct VLine* data;                      // Pointer to the actual shape data
  struct TreapNode *left, *right;

public:
  TreapNode( VLine* );
};


// Standard procedures for operating treap
VLine* find_treap( TreapNode*, int );
TreapNode* insert_treap( TreapNode*, VLine* );
TreapNode* remove_treap( TreapNode*, int );
void inorder_treap( TreapNode*, VECTOR(VLine*)& );
void clean_treap( TreapNode* );


#endif
