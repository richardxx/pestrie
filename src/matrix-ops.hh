// Copyright 2014, Hong Kong University of Science and Technology. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

/*
 * We define the low level interface to manipulate matrix.
 * To learn how to construct the bitmap index with this interface, please see the header bit-indexer.cc
 *
 * by Xiao Xiao
 * initial: 2012.10
 */

#ifndef MATRIX_OPS_H
#define MATRIX_OPS_H

#include <cstdio>
#include <cstring>
#include "bitmap.h"

// Cmatrix = comressed matrix
// The basic element in the bitmap index
class Cmatrix
{
public:
  int n, m;            // #rows, #cols (input matrix size)
  bitmap *mat;         // the matrix itself
  int *r_reps;          // the row representatives
  int *c_reps;          // the column representatives 
  int n_r_reps, n_c_reps;    // #row reps, #col reps (current matrix size)

public:
  Cmatrix()
  {
    // others are waiting to be set
    n = m = 0;
  }
  
  Cmatrix( int row, int col, 
	   bool alloc_vectors = true, bool alloc_matrix = true )
  {
    n = row; m = col;
    
    // Allocate vectors
    if ( alloc_vectors ||
	 alloc_matrix )
      mat = new bitmap[row];

    // Allocate the matrix
    if ( alloc_matrix ) {      
      for ( int i = 0; i < row; ++i ) {
	mat[i] = BITMAP_ALLOC(NULL);
      }
    }
    else if ( alloc_vectors )
      std::memset( mat, 0, sizeof(void*) * row );

    // The non-nullness means we have compressed the matrix for its rows/columns
    r_reps = NULL;
    c_reps = NULL;
    n_r_reps = row;
    n_c_reps = col;
  }
  
  ~Cmatrix()
  {
    if ( mat != NULL ) {
      for ( int i = 0; i < n; ++i )
	BITMAP_FREE( mat[i] );
      delete[] mat;
    }

    if ( r_reps != NULL ) delete[] r_reps;
    if ( c_reps != NULL ) delete[] c_reps;
  }

  bitmap operator[]( int i ) const
  {
    return mat[i];
  }

  bitmap at( int i ) const
  {
    return mat[i];
  }

  void set( int i, bitmap bn )
  {
    mat[i] = bn;
  }
};

// The following function must be called before using the library
extern bool
__init_matrix_lib();

// We compress the rows of the matrix
extern void
compress_equivalent_rows( Cmatrix* );

// We compress the columns of the matrix
// Because we have to compute the transpose of the input matrix A first,
// therefore, we return the transpose to the user for future use 
extern Cmatrix*
compress_equivalent_columns( Cmatrix* );

// Compute a transpose of the input matrix A
extern Cmatrix*
transpose( Cmatrix* );

/*
 * We take two matrices A and B and we multiply them.
 * We don't know how to take advantage of the matrix with compressed columns yet.
 * Therefore, the input matrices A and B cannot be compressed by columns.
 */
extern Cmatrix* 
matrix_mult( Cmatrix*, Cmatrix* );

// Decide if two matrices are exactly same in the content
extern bool
matrix_equal_p( Cmatrix*, Cmatrix* );

// Show basic information of a matrix
// Returns #bits in the matrix
extern int
profile_matrix( Cmatrix*, const char*, std::FILE* );


// Persist the matrix
extern void
serialize_out( Cmatrix*, FILE*, bool do_compress = true );

#endif
