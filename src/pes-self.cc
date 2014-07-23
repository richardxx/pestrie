// Copyright 2014, Hong Kong University of Science and Technology. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

/*
 * Implementation of generating Pestrie index for points-to information.
 *
 * by Xiao Xiao
 * initial: 2012.9
 */
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <climits>
#include "pestrie.hh"
#include "profile_helper.h"

using namespace std;

// Try to find a good ordering of the matrix rows in order to minimize PesTrie size
void 
PesTrieSelf::self_permute_rows()
{
  int i, k;
  unsigned x;
  bitmap_iterator bi;

  // Read the address of the data structures
  //int n = pestrie->n;
  int cm = this->cm;
  bitmap* mat_T = this->mat_T;
  MatrixRow* r_order = this->r_order;
  int *r_count = this->r_count;
  
  int permute_way = this->pes_opts->permute_way;

  if ( permute_way != SORT_BY_RANDOM ) {
    if ( permute_way == SORT_BY_HUB_DEGREE ) {      
      for ( i = 0; i < cm; ++i ) {
	long wt = 0;
	EXECUTE_IF_SET_IN_BITMAP( mat_T[i], 0, x, bi ) {
	  // We square the points-to size of x
	  long c = r_count[x];
	  wt += c * c;
	}
	
	r_order[i].wt = wt;
      }
    }
    else if ( permute_way == SORT_BY_SIZE ) {
      for ( i = 0; i < cm; ++i ) {
	// number of elements for each row in the input matrix
	r_order[i].wt = bitmap_count_bits( mat_T[i] );
      }
    }
    
    // Sort the rows of the input matrix by the weights from largest to smallest
    stable_sort( r_order, r_order + cm );
  }
  else {
    //random order
    for ( i = 0; i < cm - 1; ++i ) {
      MatrixRow a = r_order[i];
      int inx = rand() % ( cm - i ) + i;
      r_order[i] = r_order[inx];
      r_order[inx] = a;
    }
  }
}

// Executed before PesTrie construction
void
PesTrieSelf::preprocess()
{
  self_permute_rows();
}


void 
PesTrieSelf::basic_profile_pestrie()
{
  fprintf( stderr, "\n----------Pestrie Profile------------\n" );
  show_res_use( "PesTrie indexing" );

  int cm = this->cm;
  int vn = this->vn;
  int *es_size = this->es_size;
  int *pes = this->pes;
  vector<CrossEdgeRep*> *cross_edges = this->cross_edges;

  // We count the number of equivalent sets
  int n_es_pointers = this->vn;
  int n_cross = 0;
  for ( int i = 0; i < cm; ++i ) {
    int sz = es_size[i];
    // Some root nodes contain only the objects 
    if ( sz == 1 ) n_es_pointers--;
    n_cross += cross_edges[i].size();
  }
  
  fprintf( stderr, "PesTrie : Trees = %d, Nodes (Contain Pointers) = %d (%d), Edges (Cross Edges) = %d (%d)\n",
	   cm, 
	   vn, n_es_pointers, 
	   n_cross + vn - cm, n_cross );
}

/*
 * We pair up the cross edges to generate the index rectangles.
 * We store the index figures in the segment tree.
 */
int 
PesTrieSelf::build_index()
{
  int i, j, k;
  int sPrev, tr, tail;
  struct Rectangle r;
  struct CrossEdgeRep *p, *q;

  // We first obtain the pestrie information
  //int m = this->m;
  int cm = this->cm;
  int vn = this->vn;
  int *pes = this->pes;
  int *preV = this->preV;
  int *lastV = this->lastV;
  vector<int> *tree_edges = this->tree_edges;
  vector<CrossEdgeRep*> *cross_edges = this->cross_edges;

  // Then, the auxiliary data structures
  int *vis = new int[cm];
  int *Queue = new int[cm];
  CrossEdgeRep **groups = new CrossEdgeRep*[cm];   // Help classify the cross edges of the same root
  SegTree* seg_tree = build_segtree( 0, vn );
  
  memset( groups, 0, sizeof(void*) * cm );
  memset( vis, 0, sizeof(int) * cm );
  
  // For statistics use
  int n_gen_rects = 0;

  // We iteratively insert all rectangles
  for ( k = 1; k < cm; ++k ) {
    vector<CrossEdgeRep*> &treeK = cross_edges[k];
    int size = treeK.size();

    // Pair up the cross pointers and local pointers
    r.y1 = preV[k];
    r.y2 = lastV[k];
    tail = 0;

    for ( i = 0; i < size; ++i ) {
      p = treeK[i];
      sPrev = preV[ p->t ];
      r.x1 = sPrev;
      if ( p -> start == tree_edges[p->t].size() ) {
	// In this case, we cannot walk down from p->t
	// We directly set p->start to be the pre-order of p->t 
	p -> start = sPrev;
      }
      else {
	p -> start = lastV[ tree_edges[p->t][p->start] ];
      }
      r.x2 = p->start;
      seg_tree->insert_segtree(r);
      
      ++n_gen_rects;

      // Group the cross edges according to their tree values
      // Here, vis is used to mark if a particular tree has been visited
      // Queue records the set of trees makred in vis
      tr = pes[ p->t ];
      if ( vis[tr] == 0 ) {
	vis[tr] = 1;
	Queue[tail++] = tr;
      }
      p -> next = groups[ tr ];
      groups[ tr ] = p;
    }

    // We visit the cross edges by their pes labels from small to large
    // TODO: a hybrid solution that only sorts the root with more than XXX cross edges
    if ( tail == 1 ) {
      // The largest case
      // tr is still that value
      vis[tr] = 0;
      groups[tr] = NULL;
    }
    else {
      sort( Queue, Queue + tail );

      /*
       * We have the property here:
       * For two PES nodes x and y, if pes[x] < pes[y] we must have preV[x] < preV[y].
       */    
      for ( i = 0; i < tail; ++i ) {
	tr = Queue[i];
	p = groups[ tr ];
	groups[ tr ] = NULL;
	vis[tr] = 0;
	
	// Pair up every two cross edges
	while ( p != NULL ) {
	  r.x1 = preV[ p->t ];
	  r.x2 = p->start;
	  for ( j = i + 1; j < tail; ++j ) {
	    q = groups[ Queue[j] ];
	    while ( q != NULL ) {
	      r.y1 = preV[ q->t ];
	      r.y2 = q->start;
	      // Insert
	      if ( seg_tree->query_point( r.x1, r.y1 ) == false ) {
		seg_tree->insert_segtree( r );
	      }
	      ++n_gen_rects;
	      q = q -> next;
	    }
	  }
	  p = p -> next;
	}
      }
    }
  }
  
  // Assign back
  this->seg_tree = seg_tree;
  this->n_gen_rects = n_gen_rects;
  
  delete[] groups;
  delete[] Queue;
  delete[] vis;

  return 0;
}

/*
 * We hide the differences between points-to matrix and side-effect matrix via input specilization.
 * Therefore, parsing input should be specific to different matrices.
 */
PesTrie* 
self_parse_input( FILE* fp, const PesOpts* pes_opts )
{
  int i, k;
  int n, m;
  int dst;
  
  // n is the number of pointers, m is the number of objects
  fscanf( fp, "%d %d", &n, &m );
  PesTrie* pestrie = new PesTrieSelf( n, m, pes_opts );

  // Auxiliary data structures
  int *r_count = new int[n];
  memset( r_count, 0, sizeof(int) * n );

  /*
   * Format 0: every line is leaded by a number indicating the number of integers in this row.
   * Format 1: every line is ended by -1
   */
  int input_format = pes_opts->input_format;

  for ( i = 0; i < n; ++i ) {
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
      pestrie->write_bit( i, dst );
    }
    
    if ( input_format == INPUT_END_BY_MINUS_ONE )
      r_count[i] = INT_MAX - k;
  }

  // Output statistics
  fprintf( stderr, 
	   "Input matrix : Pointers = %d, Objects = %d\n", 
	   n, m );

  // We set up the processing functions
  pestrie->index_type = PT_MATRIX;
  pestrie->r_count = r_count;
  return pestrie;
}
