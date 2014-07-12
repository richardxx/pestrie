// Copyright 2014, Hong Kong University of Science and Technology. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

/*
 * Common procedures shared by all querying systems.
 *
 * by Xiao Xiao
 * initial 2014.05
 */

#ifndef QUERY_INL_H
#define QUERY_INL_H

static int
iterate_equivalent_set( VECTOR(int) &es_set, IFilter* filter )
{
  int ans = 0;
  int size = es_set.size();
  
  for ( int i = 0; i < size; ++i ) {
    int q = es_set[i];
    //    if ( filter->validate(q) )
      ans++;
  }
  
  return ans;
}

#endif
