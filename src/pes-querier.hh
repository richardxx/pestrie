/*
 * The querying data structure.
 * by richardxx, 2014.6
 */

#ifndef PES_QUERIER_H
#define PES_QUERIER_H

#include "kvec.hh"

#define VECTOR(T) kvec_t(T)
//#define VECTOR(T) vector<T>

const int SIG_POINT = 0;
const int SIG_VERTICAL = 0x04000000;
const int SIG_HORIZONTAL = 0x80000000;
const int SIG_RECT = 0xc0000000;


// Point
struct Point
{
public:
  int y1;
  
  Point() { 
    //counts = 1; 
  }
  Point( const Point& pt ): y1(pt.y1) { 
    //counts = 1; 
  }
  Point( int y ): y1(y) { 
    //counts = 1; 
  }
  
  Point& operator=( const Point& other )
  {
    y1 = other.y1;
    return *this;
  }

  virtual bool in_range( int y ) 
  {
    return y == y1;
  }
};

// Vertical Line
struct VLine : public Point
{
public:
  int y2;

  VLine() { }

  VLine( int y1, int y2 ) { this->y1 = y1; this->y2 = y2; }
  
  VLine( const VLine* other ) { this->y1 = other->y1; this->y2 = other->y2; }
  
  VLine& operator=( const VLine& other )
  {
    y1 = other.y1;
    y2 = other.y2;
    return *this;
  }

  bool in_range( int y )
  {
    return y1 <= y && y <= y2;
  }
};


// Full rectangle
struct Rectangle : public VLine
{
public:
  int x1, x2;
  
  Rectangle() { }
  Rectangle( int x1, int x2, int y1, int y2 ) {
    this->x1 = x1;
    this->x2 = x2;
    this->y1 = y1;
    this->y2 = y2;
  }
};


// Structure for query plane
struct QHeader
{
public:
  int x1, x2;
  QHeader *left, *right, *parent;
  bool merged;

  VECTOR(VLine*) rects;
  VECTOR(VLine*) vertis;
  //bitmap pointsto;

  QHeader()
  {
    left = right = NULL;
    parent = NULL;
    merged = false;
    //pointsto = BITMAP_ALLOC(NULL);
  }

  void add_rect( VLine* p ) 
  { 
    rects.push_back(p);
  }
  
  void add_vertis( VLine* p ) 
  { 
    vertis.push_back(p); 
    //bitmap_set_bit(pointsto, p);
  }

  int n_of_rects() 
  {
    return rects.size();
  }
};


#endif
