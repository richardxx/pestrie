// Copyright 2014, Hong Kong University of Science and Technology. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

/*
 * Implementation of generating Pestrie index for side-effect information (experimental).
 * 
 * by Xiao Xiao
 * initial, 2012.9
 */

#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include <climits>
#include "pestrie.hh"
#include "profile_helper.h"

using namespace std;

/*
 * We sort the first m rows by the skewness degree.
 * We directly copy the sorted order to the rest m rows.
 */
void 
PesTrieDual::dual_permute_rows()
{
  int i, k;
  unsigned x;
  bitmap_iterator bi;
  
  int n = this->n;
  int m = this->m;
  int half_m = m / 2;
  bitmap* mat_T = this->mat_T;
  MatrixRow* r_order = this->r_order;
  int *r_count = this->r_count;

  int permute_way = this->pes_opts->permute_way;

  if ( permute_way != SORT_BY_RANDOM ) {
    if ( permute_way == SORT_BY_HUB_DEGREE ) {
      for ( i = 0; i < half_m; ++i ) {
	long wt = 0;
	EXECUTE_IF_SET_IN_BITMAP( mat_T[i], 0, x, bi ) {
	  long c = r_count[x];
	  wt += c * c;
	}
	EXECUTE_IF_SET_IN_BITMAP( mat_T[i+half_m], 0, x, bi ) {
	  long c = r_count[x];
	  wt += c * c;
	}
	r_order[i].wt = wt;
      }
    }
    else if ( permute_way == SORT_BY_SIZE ) {
      for ( i = 0; i < half_m; ++i ) {
	// number of elements for each row in the input matrix
	r_order[i].wt = bitmap_count_bits( mat_T[i] );
	r_order[i].wt += bitmap_count_bits( mat_T[i+half_m] );
      }
    }
    
    // Sort the rows of the input matrix by the weights from largest to smallest
    stable_sort( r_order, r_order + half_m );
  }
  else {
    //random order
    for ( i = 0; i < half_m - 1; ++i ) {
      MatrixRow a = r_order[i];
      int inx = rand() % ( half_m - i ) + i;
      r_order[i] = r_order[inx];
      r_order[inx] = a;
    }
  }

  /*
  if ( pes_opts -> obj_merge == true ) {
    // We sort the store matrix part here
    for ( i = half_m; i < m; ++i )
      r_order[i].wt = bitmap_count_bits( mat_T[i] );

    // must be stable
    stable_sort( r_order, r_order + half_m );
  }
  else {
    
  }
  */

  // For the second half, we directly copy the first half order
  for ( i = half_m; i < m; ++i )
    r_order[i].id = r_order[i-half_m].id + half_m;
}

// Collect input matrix information to optimize the PesTrie construction and querying efficiency
void 
PesTrieDual::preprocess()
{
  dual_permute_rows();
}


void 
PesTrieDual::basic_profile_pestrie()
{
  fprintf( stderr, "\n----------Pestrie Profile------------\n" );
  show_res_use( "PesTrie indexing" );

  int m = this->m;
  int half_m = m / 2;
  int vn = this->vn;
  int *es_size = this->es_size;
  int *lastV = this->lastV;
  vector<CrossEdgeRep*> *cross_edges = this->cross_edges;
  int n_cross = 0;

  // Equivalent sets for stores
  int n_es_stores = lastV[half_m-1] + 1;
  for ( int i = 0; i < half_m; ++i ) {
    int sz = es_size[i];
    // Some root nodes contain only the objects 
    if ( sz == 1 ) n_es_stores--;
    n_cross += cross_edges[i].size();
  }

  // Equivalent sets for loads
  int n_es_loads = lastV[m-1] - lastV[half_m-1];
  for ( int i = half_m; i < m; ++i ) {
    int sz = es_size[i];
    // Some root nodes contain only the objects 
    if ( sz == 1 ) n_es_loads--;
    n_cross += cross_edges[i].size();
  }
  
  fprintf( stderr, "PesTrie : Trees = %d, Nodes = %d, Edges (Cross Edges) = %d (%d)\n",
	   m,
	   vn, 
	   n_cross + vn - m, n_cross );

  fprintf( stderr, "PesTrie : ES of stores = %d, ES of loads = %d\n",
	   n_es_stores, n_es_loads );
}


/*
 * We pair up the cross edges from the two halves of the input PesTrie.
 * Then, for every root r in the first half with corresponding r' in the second half, we pair up their cross edges.  
 */
int 
PesTrieDual::build_index()
{
  int i, j, k;
  struct Rectangle r, rr;
  struct CrossEdgeRep *p, *q;

  // We first retrieve the PesTrie information
  int m = this->m;
  int half_m = m / 2;
  vector<int> *tree_edges = this->tree_edges;
  vector<CrossEdgeRep*> *cross_edges = this->cross_edges;
  int *pes = this->pes;
  int *preV = this->preV;
  int *lastV = this->lastV;
  
  // Then, the auxiliary data structures
  SegTree* seg_tree = build_segtree( 0, this->vn );

  // For statistics use
  int n_gen_rects = 0;

  // We iteratively insert all rectangles
  for ( k = 0; k < half_m; ++k ) {
    int trA = k;
    int trB = k + half_m;
    
    // Test if we have paired these two trees
    //if ( trA == -1 || trB == -1 ) continue;
    /*
    r.x1 = preV[trA];
    r.y1 = preV[trB];
    if ( query_point( seg_tree, r.x1, r.y1 ) == true )
      continue;
    */

    // We first update the \xi conditions of all the cross edges of trA and trB
    int trees[] = {trA, trB};
    for ( i = 0; i < 2; ++i ) {
      int tr = trees[i];
      for ( vector<CrossEdgeRep*>::iterator it = cross_edges[tr].begin(),
	      ie = cross_edges[tr].end(); it != ie; ++it ) {
	p = *it;
	if ( p -> start == tree_edges[p->t].size() ) {
	  // In this case, we cannot walk down from p->t
	  p->start = preV[p->t];
	}
	else 
	  p->start = lastV[ tree_edges[p->t][p->start] ];
      }
    }

    int size1 = cross_edges[trA].size();
    int size2 = cross_edges[trB].size();
 
    // We generate the store-load conflicts
    for ( i = -1; i < size1; ++i ) {
      if ( i == -1 ) {
	// We use the subtree A to pair up
	r.x1 = preV[trA];
	r.x2 = lastV[trA];
      }
      else {
	p = cross_edges[trA][i];
	r.x1 = preV[ p->t ];
	r.x2 = p->start;
      }

      // Enumerate a cross edge on tree B
      for ( j = -1; j < size2; ++j ) {
	if ( j == -1 ) {
	  // We use the subtree B to pair up
	  r.y1 = preV[trB];
	  r.y2 = lastV[trB];
	}
	else {
	  p = cross_edges[trB][j];
	  r.y1 = preV[ p->t ];
	  r.y2 = p->start;
	}
	
	// Insert
	++n_gen_rects;
	if ( seg_tree->query_point( r.x1, r.y1 ) == false ) {
	  seg_tree->insert_segtree( r );
	}
      }
    }
    
    // Then we generate the store store conflicts
    for ( vector<CrossEdgeRep*>::iterator it1 = cross_edges[trA].begin(),
	    ie = cross_edges[trA].end(); it1 != ie; ++it1 ) {
      p = *it1;
      int targetT1 = pes[p->t];

      // We fill two templates
      r.x1 = preV[p->t];
      r.x2 = p->start;
      rr.y1 = r.x1;
      rr.y2 = r.x2;

      // case-1
      r.y1 = preV[trA];
      r.y2 = lastV[trA];
      seg_tree->insert_segtree( r );
      ++n_gen_rects;

      // case-2
      for ( vector<CrossEdgeRep*>::iterator it2 = it1 + 1; it2 != ie; ++it2 ) {
	p = *it2;
	int targetT2 = pes[p->t];
	if ( targetT1 == targetT2 ) continue;
	
	Rectangle *pr;
	if ( targetT1 < targetT2 ) {
	  r.y1 = preV[p->t];
	  r.y2 = p->start;
	  pr = &r;
	}
	else {
	  rr.x1 = preV[p->t];
	  rr.x2 = p->start;
	  pr = &rr;
	}

	// Insert
	++n_gen_rects;
	if ( seg_tree->query_point( pr->x1, pr->y1 ) == false )
	  seg_tree->insert_segtree( *pr );
      }      
    }
  }

  this->seg_tree = seg_tree;
  this->n_gen_rects = n_gen_rects;
  return 0;
}

/*
 * We hide the differences between points-to matrix and side-effect matrix via input transformation.
 * Specifics for side-effect matrix:
 * 
 * We duplicate an object to two shadow copies:
 * Om: the copy that is writen by store statements;
 * Or: the copy that is read by load statements.
 *
 * The input matrix is divided into two submatrices:
 * 1. The first m rows correspond to Om.
 * 2. The next m rows correspond to Or.
 */
PesTrie* 
dual_parse_input( FILE* fp, const PesOpts* pes_opts )
{
  int i, k;
  int n, m;
  int nl, ns;
  int type, dst;
  
  // n is the number of pointers, m is the number of objects
  fscanf( fp, "%d %d", &n, &m );
  PesTrie* pestrie = new PesTrieDual( n, m + m, pes_opts );

  // Auxiliary data structures
  int *r_count = new int[n];
  memset( r_count, 0, sizeof(int) * n );
  
  /*
   * Format 0: every line is leaded by a number indicating the number of integers in this row.
   * Format 1: every line is ended by -1
   */
  int input_format = pes_opts->input_format;
  nl = ns = 0;

  for ( i = 0; i < n; ++i ) {
    // the mod/ref flag
    fscanf( fp, "%d", &type );
    if ( type == SE_STORE ) ++ns;
    else ++nl;

    if ( input_format == INPUT_START_BY_SIZE ) {
      fscanf( fp, "%d", &k );
      r_count[i] = k;
    }
    else
      // set a large value to k
      k = INT_MAX;
    
    while ( k > 0 ) {
      fscanf( fp, "%d", &dst );
      if ( dst == -1 ) break;

      k--;
      // We distinguish the same memory location via the MOD/REF flag
      if ( type == SE_LOAD ) dst += m;
      pestrie->write_bit( i, dst );
    }

    if ( input_format == INPUT_END_BY_MINUS_ONE )
      r_count[i] = INT_MAX - k;
  }

  // Output statistics
  fprintf( stderr, 
	   "Input matrix : Stores = %d, Loads = %d, Fields = %d\n", 
	   ns, nl, m );

  // We set up the processing functions
  pestrie->index_type = SE_MATRIX;
  pestrie->r_count = r_count;
  return pestrie;
}
