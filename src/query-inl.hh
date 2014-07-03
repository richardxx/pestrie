/*
 * Common procedures shared by all querying systems.
 * By richardxx, 2014.05
 */

#ifndef QUERY_INL_H
#define QUERY_INL_H

static int
iterate_equivalent_set( VECTOR(int) *es_set, IFilter* filter )
{
  int ans = 0;
  int size = es_set->size();
  
  for ( int i = 0; i < size; ++i ) {
    int q = es_set->at(i);
    if ( filter->validate(q) )
      ans++;
  }
  
  return ans;
}

#endif
