// Copyright 2014, Hong Kong University of Science and Technology. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

/*
 * An sparse bitmap based persistence scheme for pointer information.
 * The equivalent pointers and objects are detected by hashing.
 *
 * by Xiao Xiao
 * initial: 2009.9
 * modified: 2012.7
 * refactored & modified: 2012.10
 */

#include <cstdio>
#include <unistd.h>
#include <cstdlib>
#include "bit-index.hh"
#include "constants.hh"
#include "profile_helper.h"

static char* input_file = NULL;
static char* output_file = NULL;
static int input_format = INPUT_START_BY_SIZE;
static int matrix_type = PT_MATRIX;
static bool profile_in_detail = false;
static bool binarization = false;
static bool merging_eqls = true;


// Program options
static void
print_help(const char* prog_name) 
{
  printf( "Bitmap indexer version %s\n", BIT_VERSION );
  fprintf( stderr, "Usage : %s [options] input_file [output_file]\n", prog_name );
  fprintf( stderr, "Options  : \n" );
  fprintf( stderr, "-e       : Specify the input matrix type\n" );
  fprintf( stderr, "       0 : Points-to matrix (default).\n" );
  fprintf( stderr, "       1 : Side-effect matrix.\n" );
  fprintf( stderr, "-j       : Do not merge the equivalent pointers/objects.\n" );
  fprintf( stderr, "-B       : Directly output the input matrix in binary format. Don't make index.\n" );
  fprintf( stderr, "-F       : Specify the format of the input file\n" );
  fprintf( stderr, "       0 : Each line starts with the number of the following elements (default);\n" );
  fprintf( stderr, "       1 : Each line ends with -1.\n" );
  fprintf( stderr, "-g       : Give a comprehensive profiling of the intermediate results.\n" );
  fprintf( stderr, "-h       : Show this help.\n" );
}

static bool 
parse_options( int argc, char **argv )
{
  int c;

  while ( (c = getopt( argc, argv, "e:jF:gBh" ) ) != -1 ) {
    switch ( c ) {
    case 'e':
      matrix_type = atoi( optarg );
      break;

    case 'F':
      input_format = atoi( optarg );
      break;

    case 'j':
      merging_eqls = false;
      break;

    case 'g':
      profile_in_detail = true;
      break;
      
    case 'B':
      binarization = true;
      break;

    case 'h':
      print_help( argv[0] );
      return false;
      
    default:
      printf( "This program doesn't support this argument.\n" );
      break;
    }
  }

  if ( optind == argc ) {
    //fprintf( stderr, "Missing input file, cannot proceed. \n" );
    print_help( argv[0] );
    return false;
  }

  if ( matrix_type != PT_MATRIX &&
       matrix_type != SE_MATRIX ) {
    fprintf( stderr, "Wrong matrix type.\n" );
    return false;
  }

  input_file = argv[optind];
  output_file = NULL;

  ++optind;
  if ( optind < argc ) {
    output_file = argv[optind];
  }

  return true;
}

static BitIndexer*
input()
{
  FILE *fp;

  fp = fopen( input_file, "r" );
  if ( fp == NULL ) {
    fprintf( stderr, "Loading file failed.\n" );
    return NULL;
  }

  BitIndexer *indexer = NULL;
  fprintf( stderr, "\n------------Input:%s-----------\n", input_file );

  if ( matrix_type == PT_MATRIX ) {
    indexer = parse_points_to_input( fp, input_format );
  }
  else {
    indexer = parse_side_effect_input( fp, input_format );
  }

  fclose( fp );
  show_res_use( "Input" );

  return indexer;
}

int main( int argc, char **argv )
{
  //fprintf( stderr, "B = %d\n",BLOCK_SIZE );
  if ( parse_options( argc, argv ) == false ) 
    return -1;
  
  BitIndexer *indexer = input();
  if ( indexer == NULL ) return -1;

  if ( binarization == false ) {
    // Generate persistence
    indexer->fp_generate_index( indexer, merging_eqls );

    // Profile the result
    if ( profile_in_detail )
      indexer->fp_profile_index( indexer );
  }

  // Output
  if ( output_file != NULL ) {
    FILE *fp = fopen( output_file, "wb" );
    if ( fp == NULL ) {
      fprintf( stderr, "Cannot write to the file: %s\n", output_file );
      return -1;
    }
    indexer->fp_externalize_index( indexer, fp, binarization );
  }

  delete indexer;
  return 0;
}
