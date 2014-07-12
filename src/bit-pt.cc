// Copyright 2014, Hong Kong University of Science and Technology. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

/*
 * The sparse bitmap based points-to persistence constructor.
 * By Xiao Xiao
 * initial: 2012.10
 */

#include <cstring>
#include "bit-index.hh"
#include "profile_helper.h"
#include "constants.hh"
#include "histogram.hh"

using namespace std;

static void
compute_alias_matrix( Cmatrix** imats, bool merging_eqls )
{
  Cmatrix *ptm = imats[I_PT_MATRIX];
  Cmatrix *ptm_T = NULL;

  if ( merging_eqls ) {
    compress_equivalent_rows( ptm );
    ptm_T = compress_equivalent_columns( ptm );
  }
  else {
    ptm_T = transpose( ptm );
  }

  imats[I_ALIAS_MATRIX] = matrix_mult( ptm, ptm_T );
  imats[I_PTED_MATRIX] = ptm_T;
}


static void
generate_index( BitIndexer* pt_indexer, bool merging_eqls )
{
  Cmatrix **imats = pt_indexer->imats;

  // Generate
  compute_alias_matrix( imats, merging_eqls );
  
  fprintf( stderr, "\n-----------Points-to Index-------------\n" );
  show_res_use( "Bitmap indexing" );

  // Simple profile
  fprintf( stderr, "Input points-to matrix: Pointers = %d, Objects = %d\n", 
	   pt_indexer->n_global, pt_indexer->m_global );

  const double int_size = sizeof(int) * 1.0;
  int labels = 0;
  double total_mem = 0.0;

  for ( int i = 0; i < N_OF_LOADABLE_PT_INDEX; ++i ) {
    if ( imats[i] == NULL ) continue;
    int bits = bitmap_calculate_labels( imats[i]->mat, imats[i]->n_r_reps );
    double mem = bitmap_calculate_memory( imats[i]->mat, imats[i]->n_r_reps );
    fprintf( stderr, "Encoded %s matrix: rows = %d, columns = %d, bits = %d, mem = %.0lfkb\n", 
	     pt_matrix_info[i],
	     imats[i]->n_r_reps, 
	     imats[i]->n_c_reps, 
	     bits - imats[i]->n_r_reps,
	     (mem - imats[i]->n_r_reps*int_size) / 1024 );

    labels += bits + int_size * 3;
    total_mem += mem;
  }


  int n = pt_indexer->n_global;
  int m = pt_indexer->m_global;
  
  total_mem += (n+m) * int_size;
  labels += (n+m);

  fprintf( stderr, "Index labels: %d\n", labels );
  fprintf( stderr, "The bitmap compressed index size is : %.0lfKb\n", total_mem / 1024 );
  //fprintf( stderr, "The uncompressed index size is : %.0lfKb\n", labels * int_size / 1024 );
}


static void
advanced_profile( BitIndexer* pt_indexer )
{
  Cmatrix **imats = pt_indexer->imats;

  // Points-to matrix
  histogram pt_skew;
  long pt_scales[] = { 3, 7, 17, 45 };
  pt_skew.push_scales( pt_scales, 4 );
  Cmatrix *pm = imats[I_PT_MATRIX];
  for ( int i = 0; i < pm->n_r_reps; ++i ) {
    int sz = bitmap_count_bits( pm->mat[i] );
    pt_skew.add_sample( sz );
  }

  // Alias matrix
  histogram alias_skew;
  long as_scales[] = { 5, 17, 57, 97 };
  alias_skew.push_scales( as_scales, 4 );
  
  Cmatrix *am = imats[I_ALIAS_MATRIX];
  for ( int i = 0; i < am->n_r_reps; ++i ) {
    int sz = bitmap_count_bits( am->mat[i] );
    alias_skew.add_sample( sz );
  }

  // Pointed-to matrix
  histogram pted_skew;
  long pted_scales[] = { 3, 20, 80, 200 };
  pted_skew.push_scales( pted_scales, 4 );
  
  Cmatrix *pted = imats[I_PTED_MATRIX];
  for ( int i = 0; i < pted->n_r_reps; ++i ) {
    int sz = bitmap_count_bits( pted->mat[i] );
    pted_skew.add_sample( sz );
  }
  
  // Output
  pt_skew.print_result( stderr, "Points-to matrix size distribution", false );
  alias_skew.print_result( stderr, "Alias matrix size distribution", false );
  pted_skew.print_result( stderr, "Pointed-to matrix size distribution", false );
}


/*
 * Output format:
 * Points-to and Alias Matrix:
 *
 * N_p (pointers) N_o(objects)
 * Pointers representative map (N_p)
 * Objects representative map (N_o)
 *
 * Points-to and alias matrices are encoded as follows in order:
 * matrix type (4 byte) #rows #cols
 * k(#blocks) b1(indx, data) b2 ... bk
 */
static void
externalize_index( BitIndexer* pt_indexer, FILE *fp, bool binarization )
{
  // obtain
  Cmatrix **imats = pt_indexer->imats;
  int n = pt_indexer->n_global;
  int m = pt_indexer->m_global;

  // Write magic number
  if ( binarization == false )
    fwrite( BITMAP_PT_1, sizeof(char), 4, fp );

  // Write N_p
  fwrite( &n, sizeof(int), 1, fp );
  // Write N_o
  fwrite( &m, sizeof(int), 1, fp );

  // Write the pointer/object representatitives
  if ( binarization == false ) {
    Cmatrix *ptm = imats[I_PT_MATRIX];
    fwrite( ptm->r_reps, sizeof(int), ptm->n, fp );
    fwrite( ptm->c_reps, sizeof(int), ptm->m, fp );
  }

  // Write the points-to and alias matrix
  for ( int k = 0; k < pt_indexer->n_len; ++k ) {
    if ( imats[k] == NULL ) continue;
    fwrite( &k, sizeof(int), 1, fp );
    serialize_out( imats[k], fp, !binarization );
  }
}

BitIndexer*
parse_points_to_input( FILE *fp, int fmt )
{
  int i, k;
  int n, m;
  int dst;

  __init_matrix_lib();
  fscanf( fp, "%d %d", &n, &m );

  // create
  Cmatrix *ptm = new Cmatrix( n, m );
  bitmap *mat = ptm->mat;

  for ( i = 0; i < n; ++i ) {
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
      
      k--;
      bitmap_set_bit( mat[i], dst );
    }
  }
  
  // 0 is points-to matrix
  // 1 is reserved for alias matrix
  // 2 is reserved for pointed-to matrix
  Cmatrix **pt_set = new Cmatrix*[N_OF_PT_INDEX];
  pt_set[I_PT_MATRIX] = ptm;
  pt_set[I_ALIAS_MATRIX] = NULL;
  pt_set[I_PTED_MATRIX] = NULL;

  BitIndexer *pt_indexer = new BitIndexer;
  pt_indexer->imats = pt_set;
  pt_indexer->n_len = N_OF_LOADABLE_PT_INDEX;
  pt_indexer->n_global = n;
  pt_indexer->m_global = m;
  pt_indexer->fp_generate_index = generate_index;
  pt_indexer->fp_externalize_index = externalize_index;
  pt_indexer->fp_profile_index = advanced_profile;

  return pt_indexer;
}
