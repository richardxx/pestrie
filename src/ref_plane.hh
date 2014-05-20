/*
 * A simple reference structure for rectangles.
 * A plane is discreted to a set of vertical lines. See figure 6 of our PLDI paper.
 * 
 * For every rectangle <x1, x2, y1, y2>, we refer to this rectangle for all 
 * By richardxx, 2010.6
 * Updated, 2012.9
 */

#ifndef REFPLANE_H
#define REFPLANE_H

#include "treap.hh"


// A single vertical line of references
struct RefLine
{
  int x;
  TreapNode *rects;

  RefLine(int xx) { 
    x = xx;
    rects = NULL;
  }
};


class RefPlane
{
public:
  // X axis of this vertical scan line
  int maxN;
  // all rectangles referenced all scan lines
  RefLine *ref_lines;
  
  RefPlane(int);

  bool query_point(int, int);
  
  void insert_shape( const VLine* );
  
  void visit_shapes( ShapeVisitor* );
};

#endif

