/*
 * Unifying the building and querying structures.
 * Inserting a rectangle takes O(nlgn) time, querying is O(logn).
 * In total, this can improve the overall performance is querying is much more than insertion.
 * By richardxx, 2014.4
 */

#include "ref_plane.hh"

using namespace std;


RefPlane::RefPlane(int n)
{
  maxN = n;
  refs = new RefLine[n];

  for ( int i = 0; i < n; ++i ) {
    ref_lines[i] = new RefLine(i);
  }
}

bool RefPlane::query_point(int x, int y)
{
  RefLine* refx = ref_lines[x];
  
  VLine* r = (VLine*)find_treap( refs->rects, y );
  if ( r && y <= r->y2 ) return true;
  return false;
}


void RefPlane::insert_shape(const Rectangle& r)
{
  int x1 = r.x1;
  int x2 = r.x2;
 
  do {
    RefLine* refx = ref_line[i];
    
  }
}
