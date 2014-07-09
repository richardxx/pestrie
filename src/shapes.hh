/*
 * Define the shapes.
 * Only vertical line and rectangle are defined. 
 * We represent point and horizontal line by vertical line and rectangle representatively to obtain better querying performance.
 * By richardxx, 2014.4
 */

#ifndef SHAPE_H
#define SHAPE_H

const int SIG_POINT = 0;
const int SIG_VERTICAL = 0x40000000;
const int SIG_HORIZONTAL = 0x80000000;
const int SIG_RECT = 0xc0000000;
const int SIG_FIGURE = 0xc0000000;

// Vertical Line
// It represents both points and vertical lines
struct VLine
{
  int y1, y2;

  VLine() { }
  
  VLine( int y1, int y2 ) 
  { 
    this->y1 = y1; 
    this->y2 = y2; 
  }

  VLine( const VLine& other ) 
  {
    *this = other;
  }

  VLine& operator=( const VLine& other )
  {
    y1 = other.y1;
    y2 = other.y2;
    return *this;
  }

#ifdef INDEX_UTILITY
  virtual int get_type() { return SIG_VERTICAL; }

  virtual int prepare_labels(int* labels)
  {
    if ( y1 == y2 ) {
      labels[0] = y1;
      return 1;
    }
    
    labels[0] = y1 | SIG_VERTICAL;
    labels[1] = y2;
    return 2;
  }
#endif
};

/*
 * x1, y1: lower left corner
 * x2, y2: upper right corner
 */
struct Rectangle : public VLine
{
  int x1, x2;
  
  Rectangle() { }

  Rectangle( int X1, int Y1, int X2, int Y2 )
  {
    x1 = X1;
    y1 = Y1;
    x2 = X2;
    y2 = Y2;
  }

  Rectangle( const Rectangle& other ) 
  {
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

#ifdef INDEX_UTILITY
  int get_type() { return SIG_RECT; }

  // Virtual overloading
  int prepare_labels(int* labels)
  {
    if ( y1 == y2 ) {
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
#endif
};

#endif
