// Copyright 2014, Hong Kong University of Science and Technology. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

/*
 * Building sparse bitmap index for side-effect matrix.
 * by Xiao Xiao
 * initial: 2012.10
 */

#include <assert.h>
#include "matrix-ops.hh"
#include "bit-index.hh"
#include "profile_helper.h"
#include "constants.hh"

using namespace std;

static void
generate_index( BitIndexer* se_indexer, bool merging_eqls )
{
  // obtain the matrices
  Cmatrix **imats = se_indexer->imats;
  Cmatrix *m_store_T = imats[I_STORE_TRANS_MATRIX];
  Cmatrix *m_load_T = imats[I_LOAD_TRANS_MATRIX];

  // prepare
  Cmatrix *m_store = NULL, *m_load = NULL;
  
  if ( merging_eqls ) {
    m_store = compress_equivalent_columns( m_store_T );
    m_load = compress_equivalent_columns( m_load_T );
  }
  else {
    m_store = transpose( m_store_T );
    m_load = transpose( m_load_T );
  }
 
  // compute the store-store conflicts
  imats[I_ST_ST_MATRIX] = matrix_mult( m_store, m_store_T );

  // compute the store-load conflicts
  imats[I_ST_LD_MATRIX] = matrix_mult( m_store, m_load_T );
 
  delete m_store_T;
  delete m_load_T;
  
  // We save the store/load matrices
  imats[I_STORE_MATRIX] = m_store;
  imats[I_LOAD_MATRIX] = m_load;
  
  fprintf( stderr, "\n-------------Side-effect Index-------------\n" );
  show_res_use( "Bitmap indexing" );

  // basic profile
  fprintf( stderr, "Input Mod-Ref matrix: Stores = %d, Loads = %d, Fields = %d\n", 
	   se_indexer->n_stores, se_indexer->n_loads, 
	   se_indexer->m_global );
  
  const double int_size = sizeof(int) * 1.0;
  int labels = 0;
  double total_mem = 0.0;
  
  for ( int i = 0; i < N_OF_LOADABLE_SE_INDEX; ++i ) {
    if ( imats[i] == NULL ) continue;
    int bits = bitmap_calculate_labels( imats[i]->mat, imats[i]->n_r_reps );
    double mem = bitmap_calculate_memory( imats[i]->mat, imats[i]->n_r_reps );
    fprintf( stderr, "Encoded %s matrix: rows = %d, columns = %d, bits = %d, mem = %.0lfkb\n", 
	     se_matrix_info[i],
	     imats[i]->n_r_reps, 
	     imats[i]->n_c_reps, 
	     bits - imats[i]->n_r_reps,
	     (mem - imats[i]->n_r_reps*int_size) / 1024 );
    
    labels += bits + int_size * 3;
    total_mem += mem;
  }
  
  int n = se_indexer->n_global;

  total_mem += n * int_size;

  fprintf( stderr, "Index labels: %d\n", labels );
  fprintf( stderr, "The bitmap compressed index size is : %.0lfKb\n", total_mem / 1024 );
  //fprintf( stderr, "The uncompressed index size is : %.0lfKb\n", labels * int_size / 1024 );
}

static void
advanced_profile( BitIndexer *se_indexer )
{

}

/*
 * Store Matrix + Load Matrix + Store-Store Matrix + Store-Load Matrix:
 *
 * Output format:
 * N_p (num of loads+stores) N_o(fields)
 * Pointers map (N_p)
 * Objects map (N_o*2, for store/load respectively) (if any)
 *
 * Four matrices are given in order as follows (#rows can be inferred from the rep map):
 * matrix type (4 byte) #rows #cols
 * k(#blocks) b1(indx, data) b2 ... bk
 * .......
 */
static void
externalize_index( BitIndexer* se_indexer, FILE *fp, bool binarization )
{
  int i, k;

  // obtain
  Cmatrix **imats = se_indexer->imats;

  // Note imats[0] and imats[1] are now store and load matrices
  int n = se_indexer->n_global;
  int m = se_indexer->m_global;
  
  // Write magic number
  if ( binarization == false )
    fwrite( BITMAP_SE_1, sizeof(char), 4, fp );
  
  // Write N_p
  fwrite( &n, sizeof(int), 1, fp );
  // Write N_o
  fwrite( &m, sizeof(int), 1, fp );
  
  // Write the pointer remap functions
  if ( binarization == false ) {
    int *distribute_map = se_indexer->distribute_map;
    Cmatrix *stM = imats[I_STORE_MATRIX];
    Cmatrix *ldM = imats[I_LOAD_MATRIX];

    for ( i = 0; i < n; ++i ) {
      int v = distribute_map[i];

      if ( v < n ) {
	// store matrix
	v = stM->r_reps[v];
	assert( v < stM->n_r_reps && "oops, store map exceeds" );
      }
      else {
	v = ldM->r_reps[v-n];
	assert( v < ldM->n_r_reps && "oops, load map exceeds" );
	// we still use n as the separator
	if ( v != -1 ) 
	  v += n;
	// else
	//   // Empty set for load is set to be n + ldM->n_r_reps, which is different to store
	//   v = n + ldM->n_r_reps;
      }

      distribute_map[i] = v;
    }
    
    fwrite( distribute_map, sizeof(int), n, fp );
  }

  /*
  // Write the object remap functions
  for ( i = 0; i < 2; ++i )
    fwrite( imats[i]->c_reps, sizeof(int), m, fp );
  */
  
  // Write the matrices
  for ( k = 0; k < N_OF_LOADABLE_SE_INDEX; ++k ) {
    if ( imats[k] == NULL ) continue;
    fwrite( &k, sizeof(int), 1, fp );
    serialize_out( imats[k], fp, !binarization );
  }
}

BitIndexer*
parse_side_effect_input( FILE *fp, int fmt )
{
  int i, k;
  int dst, type;
  int n_ld, n_st;
  int n, m;
  
  __init_matrix_lib();
  fscanf( fp, "%d %d", &n, &m );

  // Loading the transpose of the side-effect matrix
  Cmatrix *store_T = new Cmatrix( m, n );
  Cmatrix *load_T = new Cmatrix( m, n );
  bitmap *mat_store_T = store_T->mat;
  bitmap *mat_load_T = load_T->mat;

  // count #load and #store statements
  n_ld = n_st = 0;
  int *distribute_map = new int[n];

  for ( i = 0; i < n; ++i ) {
    // the store/load flag
    fscanf( fp, "%d", &type );

    if ( fmt == INPUT_START_BY_SIZE ) {
      fscanf( fp, "%d", &k );
    }
    else
      // set a large value to k
      k = ((unsigned)(-1)>>1);    

    while ( k > 0 ) {
      fscanf( fp, "%d", &dst );
      if ( fmt == INPUT_END_BY_MINUS_ONE ) {
	if ( dst == -1 ) 
	  break;
      }
      
      if ( type == SE_LOAD )
	bitmap_set_bit( mat_load_T[dst], n_ld );
      else
	bitmap_set_bit( mat_store_T[dst], n_st );

      k--;
    }

    // we offset n in order to distinguish the type of statement i
    if ( type == SE_LOAD ) {
      distribute_map[i] = n + n_ld;
      ++n_ld;
    }
    else {
      distribute_map[i] = n_st;
      ++n_st;
    }
  }

  // Assign the real columns to these matrices
  store_T -> m = store_T -> n_c_reps = n_st;
  load_T -> m = load_T -> n_c_reps = n_ld;

  // 0 is for transpose store matrix
  // 1 is fore transpose load matrix
  // 2 is for store-store conflict matrix
  // 3 is for store-load conflict matrix
  Cmatrix** se_set = new Cmatrix*[N_OF_SE_INDEX];
  for ( int i = 0 ; i < N_OF_SE_INDEX; ++i ) se_set[i] = NULL;
  se_set[I_STORE_TRANS_MATRIX] = store_T;
  se_set[I_LOAD_TRANS_MATRIX] = load_T;
  
  BitIndexer *se_indexer = new BitIndexer;
  se_indexer->imats = se_set;
  se_indexer->n_len = N_OF_LOADABLE_SE_INDEX;
  se_indexer->distribute_map = distribute_map;
  se_indexer->n_global = n;
  se_indexer->m_global = m;
  se_indexer->n_stores = n_st;
  se_indexer->n_loads = n_ld;
  se_indexer->fp_generate_index = generate_index;
  se_indexer->fp_externalize_index = externalize_index;
  se_indexer->fp_profile_index = advanced_profile;

  return se_indexer;
}
