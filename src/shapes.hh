/*
 * All four shapes (points, vertial lines, horizontal lines, rectangles).
 * By richardxx, 2014.4
 */

#ifndef SHAPE_H
#define SHAPE_H

#include <vector>
using std::vector;

const int SIG_POINT = 0;
const int SIG_VERTICAL = 0x04000000;
const int SIG_HORIZONTAL = 0x80000000;
const int SIG_RECT = 0xc0000000;

struct Point
{
  int x1, y1;
  bool used;

  Point() { }
  Point( const Point& pt ): x1(pt.x1), y1(pt.y1) { }
  Point( int x, int y ): x1(x), y1(y) { }

  Point& operator=( const Point& other )
  {
    x1 = other.x1;
    y1 = other.y1;
    return *this;
  }

  Point* clone() {
    return new Point(*this);
  }

  virtual int prepare_labels(int* labels)
  {
    labels[0] = y1;
    return 1;
  }
};

// Vertical Line
struct VLine : public Point
{
  int y2;

  VLine() { }
  VLine( const VLine& other ) {
    *this = other;
  }

  VLine& operator=( const VLine& other )
  {
    x1 = other.x1;
    y1 = other.y1;
    y2 = other.y2;
    return *this;
  }

  void cast(VLine** to) { *to = this; }
  
  VLine* clone() {
    return new VLine(*this);
  }
};

// Horizontal Line
struct HLine : public Point
{
  int x2;

  HLine() { }
  HLine( const HLine& other ) {
    *this = other;
  }

  HLine& operator=( const HLine& other )
  {
    x1 = other.x1;
    y1 = other.y1;
    x2 = other.x2;
    return *this;
  }

  void cast(HLine** to) { *to = this; }
  
  HLine* clone() {
    return new HLine(*this);
  }
};


/*
 * x1, y1: lower left corner
 * x2, y2: upper right corner
 */
struct Rectangle : public VLine
{
  int x2;
  
  Rectangle() { }

  Rectangle( int X1, int Y1, int X2, int Y2 )
  {
    x1 = X1;
    y1 = Y1;
    x2 = X2;
    y2 = Y2;
  }

  Rectangle( const Rectangle& other ) {
    *this = other;
  }

  Rectangle& operator=( const Rectangle& other )
  {
    x1 = other.x1;
    y1 = other.y1;
    x2 = other.x2;
    y2 = other.y2;
    return *this;
  }

  void cast(Rectangle** to) { *to = this; }

  Rectangle* clone() {
    return new Rectangle(*this);
  }

  // Virtual overloading
  int prepare_labels(int* labels)
  {
    if ( y1 == y2 ) {
      // A vertical
      labels[0] = y1 | SIG_VERTICAL;
      labels[1] = y2;
      return 2;
    }
    else if ( x1 == x2 ) {
      // A horizontal
      labels[0] = y1 | SIG_HORIZONTAL;
      labels[1] = x2;
      return 2;
    }

    labels[0] = y1 | SIG_RECT;
    labels[1] = x2;
    labels[2] = y2;
    return 3;
  }
};

// An structure used to collect the figures
struct FigureSet
{
  vector<Point*> *points;
  vector<Rectangle*> *rects;

  FigureSet()
  {
    points = new vector<Point*>;
    rects = new vector<Rectangle*>;
  }
  
  ~FigureSet()
  {
    delete points;
    delete rects;
  }

  void clear() 
  {
    points->clear();
    rects->clear();
  }
};

#endif
