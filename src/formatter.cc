// Copyright 2014, Hong Kong University of Science and Technology. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

/*
 * Converts textual matrix input to binary format.
 * By Xiao Xiao
 * initial: 2014.2
 */

#include <cstdio>
#include <unistd.h>
#include <cstdlib>
#include "matrix_ops.hh"
#include "options.h"

using namespace std;

int matrix_type = PT_MATRIX;
int input_format = INPUT_START_BY_SIZE;
const char* input_file = NULL;
const char* output_file = NULL;


static bool 
parse_options( int argc, char **argv )
{
  int c;

  while ( (c = getopt( argc, argv, "e:" ) ) != -1 ) {
    switch ( c ) {
    case 'e':
      matrix_type = atoi( optarg );
      break;

    case 'F':
      input_format = atoi( optarg );
      break;
      
    default:
      printf( "This program doesn't support this argument.\n" );
      break;
    }
  }

  if ( optind + 2 > argc ) {
    fprintf( stderr, "Usage : %s [options] input_file output_file\n", argv[0] );
    fprintf( stderr, "Options  : \n" );
    fprintf( stderr, "-e       : Specify the input matrix type\n" );
    fprintf( stderr, "       0 : Points-to matrix (default).\n" );
    fprintf( stderr, "       1 : Side-effect matrix.\n" );
    fprintf( stderr, "-F       : Specify the format of the input file\n" );
    fprintf( stderr, "       0 : Each line starts with the number of the following elements (default);\n" );
    fprintf( stderr, "       1 : Each line ends with -1.\n" );
    return false;
  }

  if ( matrix_type != PT_MATRIX &&
       matrix_type != SE_MATRIX ) {
    fprintf( stderr, "Wrong matrix type.\n" );
    return false;
  }

  input_file = argv[optind++];
  output_file = argv[optind];
  return true;
}


Cmatrix*
read_points_to_input()
{
  int i, k;
  int n, m;
  int dst;

  FILE* fp = fopen( input_file, "r" );

  __init_matrix_lib();
  fscanf( fp, "%d %d", &n, &m );

  // create
  Cmatrix *ptm = new Cmatrix( n, m );
  bitmap *mat = ptm->mat;

  for ( i = 0; i < n; ++i ) {
    if ( input_format == INPUT_START_BY_SIZE ) {
      fscanf( fp, "%d", &k );
    }
    else
      // set a large value to k
      k = ((unsigned)(-1)>>1);
    
    while ( k > 0 ) {
      fscanf( fp, "%d", &dst );
      if ( input_format == INPUT_END_BY_MINUS_ONE ) {
	if ( dst == -1 ) 
	  break;
      }
      
      k--;
      bitmap_set_bit( mat[i], dst );
    }
  }
  
  fclose(fp);
  return ptm;
}


int main(int argc, char** argv)
{
  if ( parse_options(argc, argv) == false )
    return -1;

  Cmatrix* ptm = read_points_to_input();
  FILE* out_fp = fopen( output_file, "wb" );
  //compress_equivalent_rows(ptm);
  serialize_out( ptm, out_fp );
  fclose(out_fp);
  delete ptm;

  return 0;
}
