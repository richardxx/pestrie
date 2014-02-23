/*
 * A header that defines the interface for answering queries.
 *
 * by richardxx, 2011.8
 */

#ifndef QUERY_H
#define QUERY_H

#include <unistd.h>
#include <cstdlib>
#include "options.h"

// Query Types
#define N_QURIES 7
#define QT_RANDOM   0              // We choose one randomly
#define IS_ALIAS        1
#define LIST_POINTS_TO  2
#define LIST_POINTED_TO 3
#define LIST_ALIASES    4
#define LIST_ACC_VARS   5
#define LIST_CONFLICTS  6

static const char* query_strs[] = {
  "random",
  "IsAlias",
  "ListPointsTo",
  "ListPontedTo",
  "ListAliases",
  "ListModRef",
  "ListStores" };


static int query_type = IS_ALIAS;
static bool do_profile = false;
static bool print_answers = false;
static bool slow_mode = false;

static const char* input_file = NULL;
static const char* query_plan = NULL;
static const char* check_file = NULL;

// Program options
static void print_help( const char* prog_name )
{
  printf( "Usage : %s [options] input_file [query_plan]\n", prog_name );
  printf( "Options  : \n" );
  printf( "-c cfile : Perform sanity check by reconstructing the index with cfile.\n" );
  printf( "-p       : Print answer to the queries to stdout.\n" );
  printf( "-g       : Profile the input index.\n" );
  printf( "-t [num] : Specify the query type:\n" );
  printf( "    0    : randomly choose any of the following (default);\n" );
  printf( "    1    : alias query\n" );
  printf( "    2    : list points-to\n" );
  printf( "    3    : list pointed-to\n" );
  printf( "    4    : list aliases\n" );
  printf( "    5    : list mod/ref vars\n" );
  printf( "    6    : list store/load conflicts\n" );
  printf( "-s       : Run queries in slow mode (only for bitmap index).\n" );
}

static bool 
parse_options( int argc, char **argv )
{
  int c;

  while ( (c = getopt( argc, argv, "c:gpst:h" ) ) != -1 ) {
    switch ( c ) {
    case 'c':
      check_file = optarg;
      break;
            
    case 'g':
      do_profile = true;
      break;

    case 'p':
      print_answers = true;
      break;

    case 's':
      slow_mode = true;
      break;

    case 't':
      query_type = std::atoi( optarg );
      if ( query_type < 0 || query_type > 9 )
	query_type = QT_RANDOM;
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
    //fprintf( stderr, "Missing index file, cannot proceed. \n" );
    print_help( argv[0] );
    return false;
  }

  input_file = argv[optind];
  query_plan = NULL;

  ++optind;
  if ( optind < argc ) {
    query_plan = argv[optind];
  }

  return true;
}

#endif
