// Copyright 2014, Hong Kong University of Science and Technology. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

/*
 * A driver for constructing Pestrie index.
 *
 * by Xiao Xiao,
 * initial: 2012.9
 */

#include <cstdio>
#include <unistd.h>
#include <cstdlib>
#include "constants.hh"
#include "pestrie.hh"
#include "profile_helper.h"
#include "segtree.hh"

using namespace std;

static bool interactive_query = false;
static int input_format = 0;
static char *input_file = NULL;
static char *output_file = NULL; 
static int matrix_type = 0;
static const char* magic_numbers[] = { PESTRIE_PT_1, PESTRIE_SE_1 };


// The options for indexing programs
static void
print_help(const char* prog_name)
{
  printf( "Pestrie version %s\n", PES_VERSION );
  printf( "Usage : %s [options] input_file [output_file]\n", prog_name );
  printf( "Options  : \n" );
  printf( "-d       : Draw Pes-Trie in graphviz (default = false).\n" );
  printf( "-b [num] : Permutation of source nodes in the way of\n" );
  printf( "       0 : Sort by size;\n" );
  printf( "       1 : Sort by hub degrees (default);\n" );
  printf( "       2 : Random;\n" );
  printf( "-e [num] : Specify the format of the input matrix\n" );
  printf( "       0 : Points-to matrix (default);\n" );
  printf( "       1 : Side-effect matrix.\n" );
  printf( "-g       : Give the details of pestrie (default = false).\n" );
  printf( "-i       : interactive query.\n" );
  printf( "-m       : Disable indistinguishable objects merging.\n" );
  printf( "-F       : Specify the format of the input file\n" );
  printf( "       0 : Each line starts with the number of the following elements (default);\n" );
  printf( "       1 : Each line ends with -1.\n" );
  printf( "-l       : The input points-to information is produced by LLVM.\n" );
}

static PesOpts* 
parse_options( int argc, char **argv )
{
  int c;

  PesOpts* pes_opts = new PesOpts();
  
  while ( (c = getopt( argc, argv, "b:de:F:ighml" ) ) != -1 ) {
    switch ( c ) {
    case 'b':
      pes_opts->permute_way = atoi( optarg );
      break;

    case 'e':
      matrix_type = atoi( optarg );	
      break;

    case 'i':
      interactive_query = true;
      break;

    case 'g':
      pes_opts->profile_in_detail = true;
      break;

    case 'l':
      pes_opts->llvm_input = true;
      break;

    case 'F':
      pes_opts->input_format = atoi( optarg );
      break;

    case 'd':
      pes_opts->pestrie_draw = true;
      break;

    case 'm':
      pes_opts->obj_merge = false;
      break;

    case 'h':
      print_help(argv[0]);
      delete pes_opts;
      return NULL;
      
    default:
      printf( "This program doesn't support this argument.\n" );
      break;
    }
  }

  if ( optind == argc ) {
    //printf( "Missing input file, cannot proceed. \n" );
    print_help(argv[0]);
    delete pes_opts;
    return NULL;
  }

  if ( matrix_type != PT_MATRIX && 
       matrix_type != SE_MATRIX ) {
    printf( "Wrong input matrix type. \n" );
    delete pes_opts;
    return NULL;
  }  

  input_file = argv[optind];
  output_file = NULL;
  
  ++optind;
  if ( optind < argc ) {
    output_file = argv[optind];
  }

  return pes_opts;
}

// Read the input points-to/side-effect matrix
// We directly reverse it to obtain the pointed-to/moded-used by matrix
PesTrie* input_matrix( const PesOpts* pes_opts )
{
  FILE *fp;

  fp = fopen( input_file, "r" );
  if ( fp == NULL ) return NULL;
  fprintf( stderr, "\n---------Input: %s---------\n", input_file );

  PesTrie* pestrie = NULL;

  if ( matrix_type == PT_MATRIX )
    pestrie = self_parse_input(fp, pes_opts);
  else
    pestrie = dual_parse_input(fp, pes_opts);

  fclose( fp );    

  show_res_use( "Input" );

  return pestrie;
}


static void
execute_query( PesTrie *pestrie )
{
  int query_type;
  int x, y;
  
  SegTree *seg_tree = pestrie->seg_tree;
  int *bl = pestrie->bl;
  int *pes = pestrie->pes;
  int *preV = pestrie->preV;

  printf( "\nInput queries in the console (-1 for exit):\n" );
  printf( "Format: x y, we output IsAlias(x, y)\n" );

  while ( true ) {
    printf( ">>> ");
    
    scanf( "%d %d", &x, &y );
    if ( x == -1 ) break;

    x = bl[x];
    y = bl[y];
    int ans = false;
    if ( x != -1 && y != -1 ) {
      if ( pes[x] == pes[y] )
	ans = true;
      else {
	x = preV[x];
	y = preV[y];
	if ( x > y ) { x ^= y; y ^= x; x ^= y; }
	ans = seg_tree->query_point( x, y );
      }
    }
    printf( "(%d, %d) = %s\n", x, y, (ans == true ? "true" : "false") );
  }
}

int 
main( int argc, char** argv )
{
  PesOpts* pes_opts = NULL;
  PesTrie* pestrie = NULL;

  if ( (pes_opts = parse_options( argc, argv )) == NULL ) 
    return -1;
  
  init_pestrie();
  pestrie = input_matrix(pes_opts);
  if ( pestrie == NULL ) return -1;

  // Now we build & output the index
  build_index_with_pestrie( pestrie );
  if ( output_file != NULL ) {
    FILE *fp = fopen( output_file, "wb" );
    if ( fp == NULL )
      fprintf( stderr, "Cannot write to the file: %s\n", output_file );
    else {
      pestrie->externalize_index( fp, magic_numbers[matrix_type] );
      fclose(fp);
    }
  }

  if ( interactive_query )
    execute_query( pestrie );

  delete pes_opts;
  delete pestrie;

  return 0;
}
