/*
 * A high performance probablistic balanced tree.
 * The implementation is tailored for the PesTrie algorithm.
 * By richardxx, 2010.6
 */

#ifndef TREAP_H
#define TREAP_H

#include "shapes.hh"
#include <cstdlib>

// Finally, the treap node
struct TreapNode
{
  int rkey;                         // A random key used to maintain the heap property
  Point* data;                      // Only accept the pointer data
  struct TreapNode *left, *right;

  TreapNode( Point *r )
  {
    data = r;
    rkey = std::rand();
    left = right = NULL;
  }
  
  ~TreapNode()
  {
    if ( data != NULL ) delete data;
    if ( left != NULL ) delete left;
    if ( right != NULL ) delete right;
  }
};

typedef void (*TREAP_VISITOR)( TreapNode*, void* );

//
extern Point* find_treap( TreapNode*, int);
extern TreapNode* insert_treap( TreapNode*, Point*);
extern void visit_treap( TreapNode *, TREAP_VISITOR, void* );

#endif
