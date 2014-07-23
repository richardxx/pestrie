// Copyright 2014, Hong Kong University of Science and Technology. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

/*
 * Implementation of the basic matrix operations.
 * By Xiao Xiao
 * initial, 2012.10
 */

#include <cstring>
#include <cstdlib>
#include "matrix-ops.hh"

using namespace std;

#define N_PRIMES 6
static int primes[] = { 7901, 21977, 89107, 358907, 2002841, 4088321 };

// hash table element
// only used internally
struct hash_node
{
  int idx;
  struct hash_node *next;
};


static bool initialized = false;

// Can only be called once
bool
__init_matrix_lib()
{
  if ( initialized == true ) return true;

  bitmap_obstack_initialize(NULL);
  initialized = true;
  
  return true;
}

void
compress_equivalent_rows( Cmatrix* A )
{
  int i;
  hashval_t hv;                 // type hashval_t is defined in bitmap.h
  struct hash_node **HT;        // hash table
  struct hash_node *p, *q;

  if ( A->r_reps != NULL ) return;

  // obtain
  int n = A->n;
  int m = A->m;
  bitmap* mat = A->mat;

  // create
  int *r_reps = new int[n];
 
  // Determine the suitable hash table size
  for ( i = 0; i < N_PRIMES; ++i )
    if ( primes[i] * 5 >= n ) break;
  if ( i == N_PRIMES ) i--;
  int hash_size = primes[i];
  
  HT = new hash_node*[hash_size];
  memset( HT, 0, sizeof(void*) * hash_size );
  
  // Iterate over the rows
  for ( i = 0; i < n; ++i ) {
    if ( bitmap_empty_p( mat[i] ) ) {
      // This is an empty row
      r_reps[i] = -1;
      BITMAP_FREE( mat[i] );
      continue;
    }
    
    hv = bitmap_hash( mat[i] ) % hash_size;
    p = HT[hv];
    while ( p ) {
      if ( bitmap_equal_p( mat[i], mat[ p->idx ] ) ) break;
      p = p -> next;
    }
    
    if ( p != NULL ) {
      // Found the representative for i
      r_reps[i] = p -> idx;
      BITMAP_FREE( mat[i] );
      continue;
    }
    
    r_reps[i] = i;
    p = new hash_node();
    p -> idx = i;
    p -> next = HT[hv];
    HT[hv] = p;
  }

  // Now we cluster the rows
  int n_r_reps = 0;
  for ( i = 0; i < n; ++i ) {
    if ( r_reps[i] == i ) {
      if ( i > n_r_reps ) {
	mat[n_r_reps] = mat[i];
	mat[i] = NULL;
      }
      // We again map the original row ID to new row ID
      r_reps[i] = n_r_reps++;
    }
    else {
      if ( r_reps[i] != -1 )
	// The representative ID should also be updated
	r_reps[i] = r_reps[ r_reps[i] ];
    }
  }

  // Next we delete the hash table
  for ( i = 0; i < hash_size; ++i ) {
    p = HT[i];
    while ( p != NULL ) {
      q = p->next;
      delete p;
      p = q;
    }
  }

  // assign back
  A -> mat = mat;
  A -> r_reps = r_reps;
  A -> n_r_reps = n_r_reps;

  delete[] HT;
}

// First we compute the transpose
// Then we compress the rows of the transpose
// Finally, we cluster the column bits of the input matrix
Cmatrix*
compress_equivalent_columns( Cmatrix* A )
{
  unsigned v;
  bitmap_iterator bi; 
  
  if ( A->c_reps != NULL ) return NULL;

  // Create the column mapping from its transpose
  Cmatrix *A_T = transpose( A );
  compress_equivalent_rows( A_T );

  // copy
  int m = A->m;
  int n_c_reps = A_T->n_r_reps;
  int *c_reps = new int[m];
  memcpy( c_reps, A_T->r_reps, sizeof(int) * m );

  // Then we visit the original matrix and remove the redundant bits
  bitmap *mat = A->mat;
  int n_r_reps = A->n_r_reps;

  for ( int i = 0; i < n_r_reps; ++i ) {
    bitmap row = BITMAP_ALLOC(NULL);
    EXECUTE_IF_SET_IN_BITMAP( mat[i], 0, v, bi ) {
      if ( c_reps[v] != -1 )
	bitmap_set_bit( row, c_reps[v] );
    }
    BITMAP_FREE( mat[i] );
    mat[i] = row;
  }
  
  // Assign back
  A->c_reps = c_reps;
  A->n_c_reps = n_c_reps;
  
  return A_T;
}

Cmatrix*
transpose( Cmatrix* A )
{
  unsigned v;
  bitmap_iterator bi; 
      
  // obtain
  int n = A->n;
  int m = A->m;
  bitmap *matA = A->mat;
  int n_r_reps = A->n_r_reps;

  // create
  Cmatrix *B = new Cmatrix( m, n );
  bitmap *matB = B->mat;

  for ( int i = 0; i < n_r_reps; ++i ) {
    if ( matA[i] != NULL ) {
      EXECUTE_IF_SET_IN_BITMAP( matA[i], 0, v, bi ) {
	bitmap_set_bit( matB[v], i );
      }
    }
  }
  
  // copy the representative vectors
  if ( A->c_reps != NULL ) {
    B->r_reps = new int[m];
    B->n_r_reps = A->n_c_reps;
    memcpy( B->r_reps, A->c_reps, sizeof(int) * m );
  }
  
  if ( A->r_reps != NULL ) {
    B->c_reps = new int[n];
    B->n_c_reps = A->n_r_reps;
    memcpy( B->c_reps, A->r_reps, sizeof(int) * n );
  }

  return B;
}

Cmatrix*
matrix_mult( Cmatrix *A, Cmatrix *B )
{
  unsigned v;
  bitmap_iterator bi; 

  if ( A->m != B->n ) return NULL;

  int t_c = 0;
  if ( A->c_reps != NULL ) ++t_c;
  if ( B->r_reps != NULL ) ++t_c;

  if ( t_c != 0 && t_c != 2 ) {
    fprintf( stderr, "Multiplication requires both matrices either compressed or uncompressed.\n" );
    return NULL;
  }

  if ( t_c == 2 ) {
    if ( A->n_c_reps != B->n_r_reps ||
	 memcmp( A->c_reps, B->r_reps, A->n_c_reps) != 0 ) {
      /*
       * Why not supporting different compression strategies?
       * The simplest answer is that it CANNOT save time. Detailed analysis is as follows:
       * 1. For every row of A, we should visit all the bits BEFORE compression;
       * 2. The only time saving is if row i of A has k bits, we merge only k' (k' <= k) rows of B. But, the time saving is at the price of using an additional hash table to judge the repeat mergings;
       * 3. Maintaining the hash table is not easy because after visiting every row of A, we should clear the table or allocate a new one. This is not efficient.
       *
       * Of course, maybe it works well for certain cases. We decide to implement it in the future.
       */
      fprintf( stderr, "Multiplication requires the same compression strategy for both matrices currently.\n" );
      return NULL;
    }
  }

  // obtain
  bitmap *matA = A->mat;
  int n_r_reps = A->n_r_reps;
  bitmap *matB = B->mat;
  
  // create
  Cmatrix *C = new Cmatrix( A->n, B->m );
  bitmap *matC = C->mat;

  // compute
  for ( int i = 0; i < n_r_reps; ++i ) {
    // Because the columns of A and rows of B use the same compression strategy
    // Therefore, they always form 1-to-1 correspondence just like without compression
    EXECUTE_IF_SET_IN_BITMAP( matA[i], 0, v, bi ) {
      bitmap_ior_into( matC[i], matB[v] );
    }
  }

  // copy
  if ( A->r_reps != NULL ) {
    C->r_reps = new int[A->n];
    C->n_r_reps = A->n_r_reps;
    memcpy( C->r_reps, A->r_reps, sizeof(int) * A->n );
  }
  
  if ( B->c_reps != NULL ) {
    C->c_reps = new int[B->m];
    C->n_c_reps = B->n_c_reps;
    memcpy( C->c_reps, B->c_reps, sizeof(int) * B->m );
  }
  
  return C;
}

bool
matrix_equal_p( Cmatrix *A, Cmatrix *B )
{
  /*
  if ( A->n != B->n ||
       A->m != B-> m ) return false;
  */

  if ( A->n_r_reps != B->n_r_reps ||
       A->n_c_reps != B->n_c_reps )
    return false;

  /*
  if ( memcmp( A->r_reps, B->r_reps, A->n ) != 0 )
    return false;

  if ( memcmp( A->c_reps, B->c_reps, A->m ) != 0 )
    return false;
  */

  bitmap *matA = A->mat;
  bitmap *matB = B->mat;

  for ( int i = 0; i < A->n_r_reps; ++i )
    if ( bitmap_equal_p( matA[i], matB[i] ) == false )
      return false;

  return true;
}

int
profile_matrix( Cmatrix *A, const char *name, FILE* out )
{
  int bits = bitmap_calculate_labels( A->mat, A->n_r_reps );
  if ( bits == 0 ) bits = A->n_r_reps;

  fprintf( out, "Encoded %s Matrix: rows = %d, columns = %d, bits = %d\n", 
	   name,
	   A->n_r_reps, 
	   A->n_c_reps, 
	   bits - A->n_r_reps );
  
  return bits;
}


void
serialize_out( Cmatrix *A, FILE* fp, bool do_compress )
{
  bitmap *mat = A->mat;
  int n_r_reps = A->n_r_reps;
  int n_c_reps = A->n_c_reps;
    
  fwrite( &n_r_reps, sizeof(int), 1, fp );
  fwrite( &n_c_reps, sizeof(int), 1, fp );
  
  for ( int i = 0; i < n_r_reps; ++i )
    bitmap_write_out( mat[i], fp, 
		      do_compress ? COMPRESSED_FORMAT : UNCOMPRESSED_FORMAT );
}
