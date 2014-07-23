// Copyright 2014, Hong Kong University of Science and Technology. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

/*
 * Implementation of the Treap balanced tree.
 *
 * by Xiao Xiao
 * initial: 2012.9
 */

#include "treap.hh"
#include <cstdlib>

using namespace std;

TreapNode::TreapNode( VLine *r )
{
  data = r;
  rkey = std::rand();
  left = right = NULL;
}

static struct TreapNode* 
rotate_left( struct TreapNode *p )
{
  struct TreapNode *t;
  
  t = p -> right;
  p -> right = t -> left;
  t -> left = p;
  
  return t; 
}

static struct TreapNode* 
rotate_right( struct TreapNode *p )
{
  struct TreapNode *t;

  t = p -> left;
  p -> left = t -> right;
  t -> right = p;

  return t;
}

static TreapNode*
remove_node( TreapNode* p )
{
  // case 1
  if ( p->left == NULL ) {
    TreapNode* right = p->right;
    delete p;
    return right;
  }
  
  // case 2
  if ( p->right == NULL ) {
    TreapNode* left = p->left;
    delete p;
    return left;
  }
  
  // case 3
  if ( p->left->rkey <= p->right->rkey ) {
    // Using left child as the new root
    p = rotate_right(p);
    p->right = remove_node( p->right ); 
  }
  else {
    p = rotate_left(p);
    p->left = remove_node( p->left );
  }
  
  return p;
}

// -------------------------------------------------------

// A Non-recursive procedure to find the nearest figure below y
struct VLine* 
find_treap( struct TreapNode *p, int y )
{
  int diff;
  struct TreapNode *ans = NULL;

  while (p) {
    // We compare with the lower bound
    diff = p->data->y1 - y;
    if ( diff <= 0 ) {
      ans = p;
      if ( diff == 0 ) break;
      p = p -> right;
    }
    else {
      p = p -> left;
    }
  }
  
  return ans ? ans->data : NULL;
}

// We sort the treap nodes by the lower y axis of the contained shapes
struct TreapNode* 
insert_treap( struct TreapNode* p, VLine *r )
{
  if ( p == NULL ) {
    // Now we create a new tree node
    p = new TreapNode(r);
  }
  else {
    if ( r->y1 < p->data->y1 ) {
      p->left = insert_treap( p->left, r );
      if ( p->left->rkey < p->rkey ) p = rotate_right( p );
    }
    else {
      p->right = insert_treap( p->right, r );
      if ( p->right->rkey < p->rkey ) p = rotate_left( p );
    }
  }
  
  return p;
}

// Remove the specified treap node
struct TreapNode*
remove_treap( struct TreapNode* p, int y )
{
  if ( p == NULL ) return NULL;

  int y1 = p->data->y1;
  
  if ( y1 == y ) {
    p = remove_node(p);
  }
  else if ( y1 > y ) {
    p->left = remove_treap( p->left, y );
  }
  else {
    p->right = remove_treap( p->right, y );
  }
  
  return p;
}

// Inorder traversal of treap and collect the figures
void 
inorder_treap( TreapNode *p, VECTOR(VLine*) &collector )
{
  if ( p == NULL ) return;
  
  if ( p->left != NULL )
    inorder_treap( p->left, collector );
  
  collector.push_back(p->data);
  
  if ( p->right != NULL )
    inorder_treap( p->right, collector );
}


void
clean_treap( struct TreapNode* p )
{
  if ( p == NULL ) return;

  if ( p->left != NULL )
    clean_treap( p->left );

  if ( p->right != NULL )
    clean_treap( p->right );
  
  delete p;
}

