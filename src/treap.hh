
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
  VLine* data;                      // Only accept the pointer data
  struct TreapNode *left, *right;
  
  TreapNode( VLine *r )
  {
    data = r;
    rkey = std::rand();
    left = right = NULL;
  }  
};


//
VLine* find_treap( TreapNode*, int);
TreapNode* insert_treap( TreapNode*, VLine*);
TreapNode* remove_treap( TreapNode*, int);
void clean_treap( TreapNode* );


// Inorder traversal, for the reason of speeding up the index load
template<class T> void 
visit_treap( TreapNode *p, vector<T>* set )
{
  if ( p == NULL ) return;
  
  if ( p->left != NULL )
    visit_treap<T>( p->left, set );
  
  set->push_back((T)p->data);
  
  if ( p->right != NULL )
    visit_treap<T>( p->right, set );
}


#endif
