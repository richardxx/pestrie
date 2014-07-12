// Copyright 2014, Hong Kong University of Science and Technology. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

/*
 * The library interface for constructing sparse bitmap based persistence scheme.
 
 * By Xiao Xiao
 * initial 2012.10
 */

#ifndef BIT_INDEX_H
#define BIT_INDEX_H

#include <cstdio>
#include "matrix-ops.hh"

struct BitIndexer;

typedef void (*BIT_OPTIMIZE)();
typedef void (*BIT_GENERATE_INDEX)( struct BitIndexer*, bool );
typedef void (*BIT_EXTERNALIZE_INDEX)( struct BitIndexer*, std::FILE*, bool );
typedef void (*BIT_PROFILE_INDEX)( struct BitIndexer* );

class BitIndexer
{
public:
  // short for index matrices
  Cmatrix **imats;       
  // the length of the imats array
  int n_len;
  // the statements are classified into store/load two categories
  int *distribute_map;
  // rows and columns of the original input matrix
  int n_global, m_global;
  // number of load and store statements
  int n_stores, n_loads;

  BIT_GENERATE_INDEX fp_generate_index;
  BIT_EXTERNALIZE_INDEX fp_externalize_index;
  BIT_PROFILE_INDEX fp_profile_index;

  BitIndexer()
  {
    imats = NULL;
    n_len = 0;
    distribute_map = NULL;
    n_global = m_global = 0;
    n_stores = n_loads = 0;
    fp_generate_index = NULL;
    fp_externalize_index = NULL;
  }
  
  ~BitIndexer()
  {
    if ( imats != NULL ) {
      for ( int i = 0; i < n_len; ++i )
	if ( imats[i] != NULL )
	  delete imats[i];
      delete[] imats;
    }
    
    if ( distribute_map != NULL ) 
      delete[] distribute_map;
  }
};


// Parm 1 : input file handler
// Parm 2 : input file format
extern BitIndexer*
parse_points_to_input( std::FILE*, int );

extern BitIndexer*
parse_side_effect_input( std::FILE*, int );

#endif
