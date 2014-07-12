// Copyright 2014, Hong Kong University of Science and Technology. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

/*
 * The common procedures for building Pestrie persistence.
 *
 * By Xiao Xiao
 * initial, 2012.9
 */

#include <vector>
#include <cmath>
#include <ctime>
#include <cstring>
#include "segtree.hh"
#include "histogram.hh"
#include "pestrie.hh"
#include "profile_helper.h"
#include "matrix-ops.hh"

using namespace std;

/*
 * Two objects are equivalent if they are always pointed to by the same pointers.
 * We merge them in order to build less PesTrie subtrees.
 */
void 
PesTrie::merge_equivalent_rows()
{
  // obtain
  int n = this->n;
  int m = this->m;
  bitmap* mat_T = this->mat_T;
  MatrixRow *r_order = this->r_order;

  // create
  int *m_rep = NULL;
  int n_reps = m;

  if ( this->pes_opts->obj_merge == true &&
       this->index_type != SE_MATRIX ) {

    // m_rep is a mapping from raw_id to aggregated_id.
    // appendix: raw_id --(many-to-1)-> aggregated_id --(1-to-1)-> sorted_id
    
    // We use the matrix library to identify the representatives
    Cmatrix *cmat = new Cmatrix( m, n, false );
    cmat->mat = mat_T;
    compress_equivalent_rows( cmat );

    // copy back
    m_rep = cmat->r_reps;
    n_reps = cmat->n_r_reps;
    
    // delete
    cmat->mat = NULL;
    cmat->r_reps = NULL;
    delete cmat;       
  }
  
  // modify
  this->cm = n_reps;
  this->m_rep = m_rep;
}


/*
 * A 3-pass scan algorithm to build the PesTrie.
 * 2-pass is also possible, but it requires one more linear space vector.
 */
void 
PesTrie::build_pestrie_core()
{
  int i, j, k;
  int es, last_vertex_num, Q_end;
  unsigned x, y;
  bitmap_iterator bi;
  
  // Obtain existing data
  int n = this->n;
  int cm = this->cm;
  bitmap* mat_T = this->mat_T;
  MatrixRow *r_order = this->r_order;

  // Allocate auxiliary data structures
  vector<int> *tree_edges = new vector<int>[n+cm];
  vector<CrossEdgeRep*> *cross_edges = new vector<CrossEdgeRep*>[cm];
  int *bl = new int[n];
  int *pes = new int[n+cm];
  int *es_size = new int[n+cm];
  int *Queue = new int[n];
  int *split = new int[n+cm];
  
  // belong is initialized by -1 to indicate those pointers point to nothing
  memset( bl, -1, sizeof(int) * n );
  // aux_array records the number of pointers that a pestrie node represents
  memset( es_size, 0, sizeof(int) * (n+cm) );
  // indicate if a pestrie node has been splitted
  memset( split, -1, sizeof(int) * (n+cm) );

  // First cm entries are reserved for the PesTrie subtree roots
  // But not every root has a non-empty subtree
  int vertex_num = cm;
  for ( k = 0; k < cm; ++k ) {
   
    i = r_order[k].id;
    // We use this lower-bound to distinguish the phase
    last_vertex_num = vertex_num;

    // First pass, scan all reachable pointers
    Q_end = 0;
    pes[k] = k;
    EXECUTE_IF_SET_IN_BITMAP(mat_T[i],0,x,bi) {
      Queue[Q_end++] = x;
      es = bl[x];
      if ( es != -1 ) {
	es_size[es]--;
      }
      else {
	// Node x is put into the same ES with the root k
	bl[x] = k;
      }
    }
    
    // the matrix is of no use any more
    // BITMAP_FREE( mat_T[i] );

    // Second pass, produce new pes-nodes and tree edges
    for ( j = 0; j < Q_end; ++j ) {
      x = Queue[j];
      es = bl[x];
      // es_size[es] == 0 means this es has nothing left, so we reuse it
      if ( es_size[es] > 0 ) {
	// This pes-node has NOT been splitted
	if ( split[es] < last_vertex_num ) {
	  // The first time we encounter this pes-node
	  pes[ vertex_num ] = pes[ es ];
	  split[es] = vertex_num;
	  // Create a tree edge
	  tree_edges[es].push_back( vertex_num );
	  vertex_num++;
	}
	bl[x] = split[es];
      }
    }
    
    // Third pass, create the cross edges
    es_size[k] = 1;        // essential, otherwise the root can be empty in the partition process
    for ( j = 0; j < Q_end; ++j ) {
      x = Queue[j];
      es = bl[x];
      if ( es_size[es] == 0 ) {
	// Creat a cross edge if this node is newly created or an empty node
	CrossEdgeRep* tp = new CrossEdgeRep;
	tp -> t = es;
	tp -> start = tree_edges[es].size();
	cross_edges[k].push_back( tp );
      }
      es_size[es]++;
    }
  }

  // The last step, we generate the interval labels of PesTrie
  int *preV = new int[vertex_num];
  int *lastV = new int[vertex_num];

  // split is used for tracking the next walkable tree edge for every ES
  for ( i = 0; i < vertex_num; ++i )
    split[i] = tree_edges[i].size() - 1;
  
  // A non-recursive version of tree traversal
  int pre_order = 0;
  for ( i = 0; i < cm; ++i ) {
    Q_end = 0;
    Queue[Q_end++] = i;
    preV[i] = pre_order++;

    while ( Q_end > 0 ) {
      x = Queue[Q_end-1];
      j = split[x];
      if ( j > -1 ) {
	// We take this edge and walk down
	y = tree_edges[x][j];
	split[x] = j - 1;
	Queue[Q_end++] = y;
	preV[y] = pre_order++;
      }
      else {
	// All the subtree nodes of x have been visited
	--Q_end;
	lastV[x] = pre_order - 1;
      }
    }
  }

  // Now we update the PesTrie descriptor
  this->vn = vertex_num;
  this->tree_edges = tree_edges;
  this->cross_edges = cross_edges;
  this->bl = bl;
  this->pes = pes;
  this->es_size = es_size;
  this->preV = preV;
  this->lastV = lastV;

  //show_res_use( "Building PesTrie" );

  delete[] Queue;
  delete[] split;
}

void PesTrie::profile_index()
{
  int n = this->n;
  int m = this->m;
  int vn = this->vn;

  int n_points = seg_tree->n_points;
  int n_vertis = seg_tree->n_vertis;
  int n_horizs = seg_tree->n_horizs;
  int n_rects = seg_tree->n_rects;
  int n_total_stored = n_points + n_vertis + n_horizs + n_rects;

  fprintf( stderr, "\n------------Pestrie Index--------------\n" );
  fprintf( stderr, "We totally generate %d figures, %d of them are indexed.\n", 
	   n_gen_rects, n_total_stored );

  fprintf( stderr, "-->%d rectangles, percentage = %.2lf\%\n", 
	   n_rects, 
	   (double)(n_rects) / n_total_stored * 100 );

  fprintf( stderr, "-->%d vertical lines, percentage = %.2lf\%\n", 
	   n_vertis, 
	   (double)(n_vertis) / n_total_stored * 100 );

  fprintf( stderr, "-->%d horizontal lines, percentage = %.2lf\%\n", 
	   n_horizs,
	   (double)(n_horizs) / n_total_stored * 100 );

  fprintf( stderr, "-->%d points, percentage = %.2lf\%\n", 
	   n_points,
	   (double)(n_points) / n_total_stored * 100 );

  //fprintf( stderr, "nlgn roughly equals to : %.0lf\n", n * log(n) );
  
  fprintf( stderr, "Rectangle pairs : %d, on average %.3lf alias pairs per rectangle\n", 
	   seg_tree->n_pairs, 
	   (double)(seg_tree->n_pairs)/(n_total_stored) );

  show_res_use( NULL );
}

void 
PesTrie::advanced_profile_pestrie()
{
  if ( this->pes_opts->profile_in_detail == false )
    return;
  
  //fprintf( stderr, "\n" );
  //fprintf( stderr, "--------------Additional Information for Index--------------\n" );

  int n = this->n;
  int m = this->m;
  int cm = this->cm;
  int vn = this->vn;
  int *bl = this->bl;
  int *pes = this->pes;
  int *es_size = this->es_size;
  bitmap* mat_T = this->mat_T;
  int *r_count = this->r_count;

  // We profile the objects (pointed-to sizes + hub degrees)
  if ( 1 ) {
    double max_wt = 0;
    double ari_avg = 0;
    double geo_avg = 0;
    
    // Hub degrees
    histogram hubD;
    long skew_scales[] = { 10, 200, 5000, 50000 }; 
    hubD.push_scales( skew_scales, 4 );
    
    // Pointed-to size distribution
    histogram pted_sizes;
    long pted_scales[] = { 10, 30, 100, 200 };
    pted_sizes.push_scales( pted_scales, 4 );
    
    int *vis = new int[vn];
    double *wts = new double[cm];
    memset( vis, 0, sizeof(int) * vn );
    bool is_llvm_input = this->pes_opts->llvm_input;

    // We recompute the hub degrees and pointed-to size
    for ( int i = 0; i < cm; ++i ) {
      unsigned x;
      bitmap_iterator bi;

      int n_bits = 0;
      long wt = 0;

      // Count bits
      EXECUTE_IF_SET_IN_BITMAP( mat_T[i], 0, x, bi ) {
	int rep_x = bl[x];
	if ( vis[rep_x] == 0 ) {
	  ++n_bits;
	  vis[rep_x] = 1;
	}
      }

      EXECUTE_IF_SET_IN_BITMAP( mat_T[i], 0, x, bi ) {
	int rep_x = bl[x];
	if ( vis[rep_x] == 1 ) {
	  long ptsize = r_count[x];
	  wt += ptsize * ptsize;
	  vis[rep_x] = 0;
	}
      }
      
      if ( !is_llvm_input || wt > 1) {
	double c = sqrt(wt);
	ari_avg += c;
	geo_avg += log2(c);
	wts[i] = c;
	if ( c > max_wt ) max_wt = c;
	hubD.add_sample( c );
      }
      
      pted_sizes.add_sample( n_bits );
    }
    
    ari_avg /= cm;
    geo_avg /= cm;
    geo_avg = pow( 2.0, geo_avg );
    
    delete[] vis;
    delete[] wts;
    
    fprintf( stderr, "\n" );
    fprintf( stderr, "Max hub degree is %.1lf.\n", max_wt );
    fprintf( stderr, "Arithmetic mean is %.1lf.\n", ari_avg ); 
    fprintf( stderr, "Geometric mean is %.1lf.\n", geo_avg );
    hubD.print_result( stderr, "Hub degrees Distribution", false );
    pted_sizes.print_result( stderr, "Pointed-to-by Matrix", false );
  }
      
  // We profile the cross edges (#cross edges for each subtree)
  int tot_cross_edges = 0;
  histogram cross_edge_size;
  long cross_scales[] = {1, 3, 17, 67};
  cross_edge_size.push_scales( cross_scales, 4 );
  
  for ( int i = 0; i < cm; ++i ) {
    int sz = this->cross_edges[i].size();
    cross_edge_size.add_sample( sz );
      tot_cross_edges += sz;
  }
  
  fprintf( stderr, "Total cross edges = %d\n", tot_cross_edges );
  cross_edge_size.print_result( stderr, "PesTrie Cross Edge Distribution", false );
  
  // Output
  //tree_size.print_result( stderr, "PesTrie Tree Size Distribution", false );
  show_res_use( NULL );
}

/*
 * Now we traverse the segment tree to generate the index file.
 * The index file is in binary form and the format is shown below:
 * 
 * Magic Number (4 bytes)
 * N_p(pointer) N_o(object) N_vn(ES)
 * preV labels (N_p+N_o)
 * (Followed is a matrix of the figures. Line X is the group figures where X1 == X)
 * n_bytes1, the labels of the figures where X1 = 1
 * n_bytes2, ....
 * .....
 * n_bytesk, ....
 */
void 
PesTrie::externalize_index( FILE* fp, const char* magic_number)
{
  int n = this->n;
  int m = this->m;
  int cm = this->cm;
  int vn = this->vn;
  int *m_rep = this->m_rep;
  int *preV = this->preV;
  int *bl = this->bl;
  MatrixRow *r_order = this->r_order;
  SegTree* seg_tree = this->seg_tree;

  // Fill up the preV mapping for input pointers and objects
  int *pre_aux = new int[n+m];
  int *obj_pos = new int[cm];

  // First are the pointers
  for ( int i = 0; i < n; ++i ) {
    int x = bl[i];
    // a pointer may not have a timestamp
    pre_aux[i] = ( x == -1 ? -1: preV[x] );
  }

  // Second are the objects
  for ( int i = 0; i < cm; ++i ) {
    int j = r_order[i].id;
    obj_pos[j] = i;
  }

  for ( int i = 0; i < m; ++i ) {
    // i'th object is aggregated into the j'th group (pes node)
    int j = ( m_rep == NULL ? i : m_rep[i] );
    // The sorted id of group j
    int k = ( j == -1 ? -1 : obj_pos[j] );
    // set the timestamp of object i
    pre_aux[i+n] = ( k == -1 ? -1 : preV[k] );   
  }
  
  // Now we start to output the index
  // Write the magic number
  fwrite( magic_number, sizeof(char), 4, fp );

  // Write N_p
  fwrite( &n, sizeof(int), 1, fp );
  // Write N_o
  fwrite( &m, sizeof(int), 1, fp );
  // Write N_vn
  fwrite( &vn, sizeof(int), 1, fp );
  // Write the preV mappings for pointers+objects
  fwrite( pre_aux, sizeof(int), n + m, fp );

  int n_labels = 3 + n + m;

  // Write the figures
  seg_tree->flush_left_shapes();
  n_labels += seg_tree->dump_figures( fp );

  // Profile
  int n_points = seg_tree->n_out_points;
  int n_vertis = seg_tree->n_out_vertis;
  int n_horizs = seg_tree->n_out_horizs;
  int n_rects = seg_tree->n_out_rects;
  int n_total_stored = n_points + n_vertis + n_horizs + n_rects;

  fprintf( stderr, "\n--------------Persistence Generation---------------\n" );
  fprintf( stderr, "We persist %d figures.\n", 
	   n_total_stored );

  fprintf( stderr, "-->%d rectangles, percentage = %.2lf\%\n", 
	   n_rects, 
	   (double)(n_rects) / n_total_stored * 100 );

  fprintf( stderr, "-->%d vertical lines, percentage = %.2lf\%\n", 
	   n_vertis, 
	   (double)(n_vertis) / n_total_stored * 100 );

  fprintf( stderr, "-->%d horizontal lines, percentage = %.2lf\%\n", 
	   n_horizs,
	   (double)(n_horizs) / n_total_stored * 100 );

  fprintf( stderr, "-->%d points, percentage = %.2lf\%\n", 
	   n_points,
	   (double)(n_points) / n_total_stored * 100 );

  // Estimate the index size (almost real)
  double intsize = sizeof(int);
  fprintf( stderr, "Index labels : %d\n", n_labels );
  fprintf( stderr, "The PesTrie index size is : %.0lfKb\n", n_labels * intsize / 1024 );

  delete[] pre_aux;
  delete[] obj_pos;
}


// ---------------------------------------------------------
// Public interface for building pestrie encoding
// ---------------------------------------------------------

// We configure the global environment
void 
init_pestrie()
{
  bitmap_obstack_initialize(NULL);
  srand( time(NULL) );
}

// Driver function
void 
build_index_with_pestrie( PesTrie* pestrie )
{
  // Merge the equivalent objects
  pestrie->merge_equivalent_rows();
  
  // First we generate a proper processing order
  pestrie->preprocess();

  // Then we construct the PesTrie
  pestrie->build_pestrie_core();

  // Finally we decompose the PesTrie and generate the index
  pestrie->build_index();

  // Output the statistics of PesTrie index
  pestrie->basic_profile_pestrie();
  pestrie->advanced_profile_pestrie();
  pestrie->profile_index();
}
