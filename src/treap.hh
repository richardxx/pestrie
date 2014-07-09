/*
 * A high performance probablistic balanced tree.
 * The implementation is tailored for the PesTrie algorithm.
 * By richardxx, 2012.10
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
