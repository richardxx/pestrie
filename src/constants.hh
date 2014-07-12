// Copyright 2014, Hong Kong University of Science and Technology. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

/*
 * Define the global options for both bitmap and PesTrie index.
 * by Xiao Xiao
 * initial: 2012.9
 */

#ifndef CONSTANTS_H
#define CONSTANTS_H

#define PES_VERSION "2.1"
#define BIT_VERSION "1.3"

// Magic numbers
#define PESTRIE_PT_1 "PTP1"  
#define PESTRIE_SE_1 "SEP1"
#define BITMAP_PT_1 "PTB1"
#define BITMAP_SE_1 "SEB1"

// Categories of the input matrix
#define UNDEFINED_MATRIX -1
#define PT_MATRIX 0
#define SE_MATRIX 1

// Categories of the statements in side-effect matrix
#define SE_STORE 1
#define SE_LOAD 2

// Categories of the input matrix representation format
#define INPUT_START_BY_SIZE 0
#define INPUT_END_BY_MINUS_ONE 1

// Ways to sort the rows for PesTrie construction.
#define SORT_BY_SIZE 0
#define SORT_BY_HUB_DEGREE 1
#define SORT_BY_RANDOM 2

// Categories of matrices in points-to/side-effect bitmap index
#define N_OF_PT_INDEX 3
#define N_OF_LOADABLE_PT_INDEX 2
#define I_PT_MATRIX 0
#define I_ALIAS_MATRIX 1
#define I_PTED_MATRIX 2           // pointed-to-by matrix

static const char* pt_matrix_info[] = { 
  "Points-to",
  "Alias",
  "Pointed-to-by" 
};
  
#define N_OF_SE_INDEX  7
#define N_OF_LOADABLE_SE_INDEX 4
#define I_STORE_MATRIX 0
#define I_LOAD_MATRIX 1
#define I_ST_ST_MATRIX 2
#define I_ST_LD_MATRIX 3
#define I_LD_ST_MATRIX 4
#define I_STORE_TRANS_MATRIX 5
#define I_LOAD_TRANS_MATRIX 6

static const char* se_matrix_info[] = { 
  "Store", 
  "Load", 
  "Store-Store", 
  "Store-Load",
  "Load-Store",
  "Store-Trans",
  "Load-Trans" 
};

#endif
