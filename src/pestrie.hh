// Copyright 2014, Hong Kong University of Science and Technology. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

/*
 * Implementation of the PesTrie based algorithm.
 * The PesTrie algorithm is basically designed for computing and querying the result of A*trans(A), when the data have skewness property.
 * It can also be extended to compute general matrix multiplication A*B.
 * We try our best to optimize the code in terms of performance and code length.
 *
 * by Xiao Xiao
 * initial 2012.7: idea and prototype
 * improve 2012.9: modularize
 */
#ifndef PESTRIE_H
#define PESTRIE_H

#include "bitmap.h"
#include "segtree.hh"
#include "constants.hh"
#include <vector>

struct PesOpts 
{
  // The format for the input matrix, read self_parse_input for details
  int input_format;   
  // The way to obtain the matrix processing order
  int permute_way;
  // Merge the indistinguishable objects?
  bool obj_merge;
  // Profiling PesTrie or not
  bool profile_in_detail;
  // Output graphviz format for PesTrie visualization
  bool pestrie_draw;
  //
  bool llvm_input;

  PesOpts()
  {
    // Setup the default values
    input_format = INPUT_START_BY_SIZE;
    permute_way = SORT_BY_HUB_DEGREE;
    obj_merge = true;
    profile_in_detail = false;
    pestrie_draw = false;
    llvm_input = false;
  }
};

struct CrossEdgeRep
{
  int t;	// the other side of the edge
  int start;	// birthday
  struct CrossEdgeRep *next;

  CrossEdgeRep()
  {
    next = NULL;
  }
};

// Every aspects of a row of the input matrix
struct MatrixRow
{
  int id;
  long wt;       // the weight of this row, used for row permutation 
  
  bool operator<( const MatrixRow& o ) const
  {
    return wt > o.wt;
  }
};


// The base class PesTrie
class PesTrie
{
public:
  // Used for self-identification
  int index_type;

  // Input matrix and its descriptions
  int n, m;                  // #rows, #columns (pt-matrix)
  bitmap *mat_T;             // transpose of the input matrix (pted-matrix)
  MatrixRow *r_order;        // the processing order of the pted-matrix
  int *m_rep, cm;            // representatives of rows of the pted-matrix 
  int *r_count;              // #non zero columns for each row of pt-matrix

  // PesTrie and its descriptions
  int vn;                  // #vertex (ES)
  std::vector<int> *tree_edges;
  std::vector<CrossEdgeRep*> *cross_edges;
  int *bl, *pes;           // The ES label and PES label of the pointers and ESes
  int *es_size;            // #pointers for every ES
  int *preV, *lastV;       // Interval labels

  // Index and descriptions
  SegTree *seg_tree;
  int n_gen_rects;

  // User provided constrols
  const PesOpts* pes_opts;

public:
  // Initialize the necessary configurations
  PesTrie( int row, int col, const PesOpts* opts ) 
  { 
    n = row; m = col; cm = col;

    // Allocate the matrix (the matrix is born in its transpose form)
    mat_T = new bitmap[col];
    r_order = new MatrixRow[col];
    for ( int i = 0; i < col; ++i ) {
      r_order[i].id = i;
      r_order[i].wt = 0;
      mat_T[i] = BITMAP_ALLOC(NULL);
    }
    
    tree_edges = NULL;
    cross_edges = NULL;
    bl = NULL;
    pes = NULL;
    preV = NULL;
    lastV = NULL;
    seg_tree = NULL;
    pes_opts = opts;
  }
  
  ~PesTrie()
  {
    if ( mat_T != NULL ) {
      for ( int i = 0; i < m; ++i )
	BITMAP_FREE( mat_T[i] );
      delete[] mat_T;
    }

    if ( r_order != NULL ) delete[] r_order;
    
    if ( tree_edges != NULL ) delete[] tree_edges;
    if ( cross_edges != NULL ) delete[] cross_edges;
    if ( bl != NULL ) delete[] bl;
    if ( pes != NULL ) delete[] pes;
    if ( es_size != NULL ) delete[] es_size;
    if ( preV != NULL ) delete[] preV;
    if ( lastV != NULL ) delete[] lastV;
    
    if ( seg_tree != NULL ) delete seg_tree;
    pes_opts = NULL;
  }

  // The clients do not need to know the internal matrix is transposed
  void write_bit( int r, int c )
  {
    bitmap_set_bit( mat_T[c], r );
  }

public:
  // Collapse the rows filled with same data for input matrix
  void merge_equivalent_rows();

  // Construct PesTrie from input matrix
  // This is common to both points-to and side-effect matrices
  void build_pestrie_core();

  void externalize_index( FILE* fp, const char* magic_number);

  void profile_index();

  void advanced_profile_pestrie();

public:
  // PesTrie specialized processing functions
  // They should be implemented sub-classes
  virtual void preprocess() = 0;
  virtual int build_index() = 0;
  virtual void basic_profile_pestrie() = 0;
};


// This class is for indexing points-to information
// Computing A*A^T
class PesTrieSelf : public PesTrie 
{
public:
  PesTrieSelf( int row, int col, const PesOpts* opts ):
    PesTrie(row, col, opts) { }

public:
  void preprocess();
  int build_index();
  void basic_profile_pestrie();

private:
  void self_permute_rows();
};


// This class is for indexing mod-ref information
// Computing A*B
class PesTrieDual : public PesTrie 
{
public:
  PesTrieDual( int row, int col, const PesOpts* opts ):
    PesTrie(row, col, opts) { }
  
public:
  void preprocess();
  int build_index();
  void basic_profile_pestrie(); 
  
private:
  void dual_permute_rows();
};


// This function should be called first!
extern void init_pestrie();

//
extern PesTrie* self_parse_input( FILE*, const PesOpts* );
extern PesTrie* dual_parse_input( FILE*, const PesOpts* );
extern void build_index_with_pestrie( PesTrie* );

#endif
