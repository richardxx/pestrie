/*
 * This section defines the common procedures for all the indexing algorithms.
 * The functionality ranges from generating random queries, parsing input options to ...
 * by richardxx, 2011.7
 */

#ifndef COMMON_H
#define COMMON_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include "profile_helper.h"


static int prt_pairs = 0;
static int n_rand_query = 0;
static char *input_file = NULL;
static char *output_file = NULL;
static int pes_trie_draw = 0;
static int prune_before = -1;
static int pos_opt = 1;
static int build_mode = 0;
static int n, m;


static int parse_options( int argc, char **argv )
{
  int c;

  while ( (c = getopt( argc, argv, "pb:P:q:doh" ) ) != -1 ) {
    switch ( c ) {
    case 'p':
      prt_pairs = 1;
      break;
      
    case 'q':
      n_rand_query = atoi(optarg);
      if ( n_rand_query < 0 ) n_rand_query = 100000;
      break;

#ifdef PES_TRIE_BASED
    case 'd':
      pes_trie_draw = 1;
      break;

    case 'b':
      build_mode = atoi( optarg );
      break;

    case 'P':
      prune_before = atoi( optarg );
      break;
#endif

    case 'h':
      printf( "Usage : %s [options] input_file [output_file]\n", argv[0] );
      printf( "Options  : \n" );
#ifdef PES_TRIE_BASED
      printf( "-d       : Draw Pes-Trie in graphviz\n" );
      printf( "-b [num] : Permutation of source nodes\n" );
      printf( "       0 : Sort by size (default);\n" );
      printf( "       1 : Random;\n" );
      printf( "-P [num] : Apply the prunning strategy before [num]'th tree (default: K/10)\n" );
#endif
#ifndef TC_NAIVE_BASED
      printf( "-p       : Print all co-reachable pairs\n" );
#endif
      printf( "-q [num] : Generate [num] random queries. Measure both table lookup time and query on demand time.\n" );

      exit(0);
      break;
      
    default:
      printf( "This program doesn't support this argument.\n" );
      break;
    }
  }

  if ( optind == argc ) {
    fprintf( stderr, "Missing input file, cannot proceed. Please specify.\n" );
    return 0;
  }

  input_file = argv[optind];
  output_file = NULL;

  ++optind;
  if ( optind < argc ) {
    output_file = argv[optind];
  }

  return 1;
}


// A binary research procedure for answering the bulk ListAliases Query
// At each segment tree node, we binary search for the two pointers tha are closeset to the middle line
// then we only use those two pointers to search for the intersected rectangles
// It would be significantly faster than the naive implementation that calls ListAliases for all queried pointers.
static void 
recursive_enumerate_alias_or( segTreeNode* p, int s, int e, kvec_t(int)& query_points )
{
  Rectangle* r;
  int pi, pos;
  Rectangle **it, **ie;

  if ( p == NULL || s >= e ) return;

  it = p->sortByX1.begin();
  ie = p->sortByX1.end();

  if ( e - s > 22 ) {
    // We use binary search
    int begin = s, end = e;
    while ( end - begin > 1 ) {
      pi = (begin+end) / 2;
      if (query_points[pi] <= p->mid)
	begin = pi;
      else
	end = pi - 1;
    }
    
    if ( query_points[begin] <= p->mid )
      pi = begin + 1;
    else
      pi = begin;
  }
  else {
    // Linear scan may be faster
    for ( pi = s; pi < e && query_points[pi] <= p->mid; ++pi );
  }

  if ( pi > s ) {
    // At least one of the querying point is to the left of the middle line
    pos = query_points[pi-1];
    while ( it != ie ) {
      r = *it;
      if ( r->x1 > pos ) break;
      internal_container.push_back(r);
      r->selected = true;
      it++;
    }

    recursive_enumerate_alias_or( p->left, s, pi, query_points );

    // Special case
    if ( pos == p->mid ) {
      recursive_enumerate_alias_or( p->right, pi+1, e, query_points );
      return;
    }
  }
  
  if ( pi < e ) {
    // Now we should go over the sorted list to see if we miss something
    it = p->sortByX2.begin();
    ie = p->sortByX2.end();
    pos = query_points[pi];
    while ( it != ie ) {
      r = *it;
      if ( r->x2 < pos ) break;
      if ( r->selected == false ) {
	internal_container.push_back(r);
	r->selected = true;
      }
      it++;
    }
  
    recursive_enumerate_alias_or( p->right, pi, e, query_points );
  }
}

int ListAliasesBulk( kvec_t(int)& query_points )
{
  int ans;

  sort( query_points.begin(), query_points.end() );
  internal_container.clear();
  recursive_enumerate_alias_or( segRoot, 0, query_points.size(), query_points );

  // Decode the answers
  ans = 0;
  memset( vis, 0, sizeof(int) * vertex_num ); 
  
  for ( int i = 0; i < internal_container.size(); ++i ) {
    Rectangle *r = internal_container[i];
    r->selected = false;
    for ( int j = r->y1; j <= r->y2; ++j )
      if ( vis[j] == 0 ) {
	++ans;
	vis[j] = 1;
      }
  }

  return ans;
}

#endif
