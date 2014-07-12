// Copyright 2014, Hong Kong University of Science and Technology. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

/*
 * Define the shapes.
 * Only vertical line and rectangle are defined. 
 * We represent point and horizontal line by vertical line and rectangle representatively to obtain better querying performance.
 
 * by Xiao Xiao
 * initial: 2012.7
 * improve: 2014.2
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

  VLine() 
  { 
  }
  
  VLine( int Y1, int Y2 ) 
  { 
    y1 = Y1; 
    y2 = Y2; 
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

  virtual bool merge(VLine* other)
  {
    if ( other->get_type() != SIG_VERTICAL ) return false;
    if ( other->y1 == y2 + 1 ) {
      y2 = other->y2;
      return true;
    }

    return false;
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

  Rectangle( int X1, int X2, int Y1, int Y2 )
  {
    x1 = X1;
    x2 = X2;
    y1 = Y1;
    y2 = Y2;
  }

  Rectangle( const Rectangle& other ) 
  {
    *this = other;
  }

  Rectangle& operator=( const Rectangle& other )
  {
    x1 = other.x1;
    x2 = other.x2;
    y1 = other.y1;
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

  bool merge(VLine* other)
  {
    if ( other->get_type() != SIG_RECT ) return false;

    Rectangle* o = (Rectangle*)other;
    if ( o->x1 == x1 && o->x2 == x2 
	 && o->y1 == y2 + 1 ) {
      y2 = o->y2;
      return true;
    }
    
    return false;
  }
#endif
};

#endif
