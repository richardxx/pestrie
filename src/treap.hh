/*
 * A high performance probablistic balanced tree.
 * The implementation is tailored for the PesTrie algorithm.
 * By richardxx, 2010.6
 */

#ifndef TREAP_H
#define TREAP_H

#include <cstdlib>

struct Point
{
  int x1, y1;
  //bool selected;

  Point() { }
  Point( const Point& pt ): x1(pt.x1), y1(pt.y1) { }
  Point( int x, int y ): x1(x), y1(y) { }

  Point& operator=( const Point& other )
  {
    x1 = other.x1;
    y1 = other.y1;
    return *this;
  }
};

// Vertical Line
struct VLine : public Point
{
  int y2;

  VLine() { }

  VLine& operator=( const VLine& other )
  {
    x1 = other.x1;
    y1 = other.y1;
    y2 = other.y2;
    return *this;
  }

};

// Horizontal Line
struct HLine : public Point
{
  int x2;

  HLine() { }
  HLine& operator=( const HLine& other )
  {
    x1 = other.x1;
    y1 = other.y1;
    x2 = other.x2;
    return *this;
  }
};


/*
 * x1, y1: lower left corner
 * x2, y2: upper right corner
 */
struct Rectangle : public VLine
{
  int x2;
  
  Rectangle& operator=( const Rectangle& other )
  {
    x1 = other.x1;
    y1 = other.y1;
    x2 = other.x2;
    y2 = other.y2;
    return *this;
  }

  Rectangle() { }

  Rectangle( int X1, int Y1, int X2, int Y2 )
  {
    x1 = X1;
    y1 = Y1;
    x2 = X2;
    y2 = Y2;
    //selected = false;
  }
};

// Finally, the treap node
struct TreapNode
{
  int rkey;             // A random key used to maintain the heap property
  Point* data;              // Only accept the pointer data
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
