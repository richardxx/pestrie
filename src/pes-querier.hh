/*
 * The querying data structure.
 * by richardxx, 2014.6
 */

#ifndef PES_QUERIER_H
#define PES_QUERIER_H

#include "shapes.hh"

// Structure for query plane
struct QHeader
{
public:
  int x1, x2;
  QHeader *left, *right, *parent;
  bool merged;

  VECTOR(VLine*) rects;
  VECTOR(VLine*) vertis;
  //VECTOR(int) pointsto;

  QHeader()
  {
    left = right = NULL;
    parent = NULL;
    merged = false;
  }

  void add_rect( VLine* p ) 
  { 
    rects.push_back(p);
  }
  
  void add_vertis( VLine* p ) 
  { 
    vertis.push_back(p);
  }

  int n_of_rects() 
  {
    return rects.size();
  }
};


#endif
