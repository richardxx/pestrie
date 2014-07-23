// Copyright 2014, Hong Kong University of Science and Technology. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

/*
 * The querying system for the sparse bitmap points-to/side-effect index.
 *
 * by Xiao Xiao,
 * initial: 2011.9
 * enhancement: 2012.7
 * refactor: 2012.10
 */

#include <cstdio>
#include <cstring>
#include <set>
#include "options.hh"
#include "profile_helper.h"
#include "matrix-ops.hh"
#include "query.hh"
#include "query-inl.hh"

using namespace std;

class BitQS : public IQuery
{
public:
  // Interface functions
  bool IsAlias( int x, int y);
  int ListPointsTo( int x, IFilter* filter );
  int ListAliases( int x, IFilter* filter );
  int ListPointedBy( int o, IFilter* filter );
  int ListModRefVars( int x, IFilter* filter );
  int ListConflicts( int x, IFilter* filter );

public:
  int getPtrEqID(int x) { return pt_map[x]; }
  int getObjEqID(int x) { return obj_map[x]; }
  int nOfPtrs() { return n; }
  int nOfObjs() { return m; }
  int getIndexType() { return index_type; }

public:
  BitQS(int n_ptrs, int n_objs, int type, int mode)
  {
    n = n_ptrs; m = n_objs;
    index_type = type;
    trad_mode = mode;

    n_of_mat = 0;
    n_ld = n_st = 0;
    n_es = 0;

    if ( type == PT_MATRIX ) {
      n_of_mat = N_OF_PT_INDEX;
      obj_map = new int[n_objs];
      es2objs = new VECTOR(int)[n_objs];
    }
    else {
      n_of_mat = N_OF_SE_INDEX;
      obj_map = NULL;
      es2objs = NULL;
    }
    
    qmats = new Cmatrix*[n_of_mat];
    es2ptrs = new VECTOR(int)[n_ptrs];
    pt_map = new int[n_ptrs];
  }

  ~BitQS()
  {
    if ( qmats != NULL ) {
      for ( int i = 0; i < n_of_mat; ++i )
	if ( qmats[i] != NULL )
	  delete qmats[i];
      delete[] qmats;
    }
    
    if ( pt_map != NULL ) delete[] pt_map;
    if ( obj_map != NULL ) delete[] obj_map;
    if ( es2ptrs != NULL ) delete[] es2ptrs;
    if ( es2objs != NULL ) delete[] es2objs;
  }
  
public:
  void load_pt_index( FILE *fp );
  void load_se_index( FILE *fp );
  void rebuild_eq_groups();

private:
  int ListStores( int x, IFilter* filter );
  int ListLoads( int x, IFilter* filter );

private:
  // Input matrices
  Cmatrix **qmats;
  int n_of_mat;

  int n, m, n_es;       // #pointers, #objects, #pointer equivalent sets
  int n_ld, n_st;       // #loads, #stores

  // Mapping pointers/objects to Eq IDs
  int *pt_map, *obj_map;

  // Mapping Eq IDs to pointers/objects
  VECTOR(int) *es2ptrs, *es2objs;

  // Options
  int index_type;
  int trad_mode;
};

//static int cnt_same_es = 0;

void
BitQS::load_pt_index( FILE *fp )
{
  // Load the mapping info
  fread( pt_map, sizeof(int), n, fp );
  fread( obj_map, sizeof(int), m, fp );

  // Load the index
  int i = 0;
  while ( i < N_OF_LOADABLE_PT_INDEX ) {
    int m_type, dim_r, dim_c;

    fread( &m_type, sizeof(int), 1, fp );
    fread( &dim_r, sizeof(int), 1, fp );
    fread( &dim_c, sizeof(int), 1, fp );
    if ( m_type != i ) i = m_type;
    
    Cmatrix *cm = new Cmatrix( dim_r, dim_c, true, false );
    bool skip = (trad_mode && m_type > I_PT_MATRIX);
    for ( int k = 0; k < dim_r; ++k ) {
      bitmap r = bitmap_read_row( fp, COMPRESSED_FORMAT, skip );
      cm->set(k, r);
    }

    profile_matrix( cm, pt_matrix_info[i], stderr );
    qmats[i++] = cm;
  }
 
  // Compute the pointed-by matrix
  qmats[I_PTED_MATRIX] = transpose( qmats[I_PT_MATRIX] );
}

void
BitQS::load_se_index( FILE *fp )
{
  int i;

  // Load the mapping info
  fread( pt_map, sizeof(int), n, fp );  

  // Profile the mapping info
  n_ld = 0;
  n_st = 0;
  for ( i = 0; i < n; ++i ) {
    if ( pt_map[i] >= n ) ++n_ld;
    else ++n_st;
  }

  // Load the index
  i = 0;
  while ( i < N_OF_LOADABLE_SE_INDEX ) {
    int m_type, dim_r, dim_c;

    fread( &m_type, sizeof(int), 1, fp );
    fread( &dim_r, sizeof(int), 1, fp );
    fread( &dim_c, sizeof(int), 1, fp );
    if ( m_type != i ) i = m_type;
    
    Cmatrix *cm = new Cmatrix( dim_r, dim_c, true, false );
    bool skip = (trad_mode && m_type > I_LOAD_MATRIX);
    for ( int k = 0; k < dim_r; ++k ) {
      bitmap r = bitmap_read_row( fp, COMPRESSED_FORMAT, skip );
      cm->set(k, r);
    }

    profile_matrix( cm, se_matrix_info[i], stderr );
    qmats[i] = cm;
    ++i;
  }
  
  // Compute the transposed store and load matrix
  qmats[I_STORE_TRANS_MATRIX] = transpose( qmats[I_STORE_MATRIX] );
  qmats[I_LOAD_TRANS_MATRIX] = transpose( qmats[I_LOAD_MATRIX] );
  qmats[I_LD_ST_MATRIX] = transpose( qmats[I_ST_LD_MATRIX] );
}

void 
BitQS::rebuild_eq_groups()
{
  // Now we create the group -> pointers mapping
  // The statements in side-effect analysis are treated specially
  for ( int i = 0; i < n; ++i ) {
    int es = pt_map[i];
    if ( es == -1 ) continue;

    if ( es >= n ) {
      // A load statement
      // We reallocate its position to make the mapping compact
      es = es - n + n_st;
      pt_map[i] = es;
    }

    if ( es > n_es ) n_es = es;
    es2ptrs[es].push_back(i);
  }

  if ( index_type == PT_MATRIX ) {
    // Decompress the equivalent objects
    for ( int i = 0; i < m; ++i ) {
      int es = obj_map[i];
      if ( es != -1 )
	es2objs[es].push_back(i);
    }
  }

  // We set the number of equivalent pointer sets to be the largest es ID + 1
  n_es = n_es + 1;
}

bool 
BitQS::IsAlias( int x, int y )
{
  x = pt_map[x];
  if ( x == -1 ) return false;
  y = pt_map[y];
  if ( y == -1 ) return false;
  if ( x == y ) return true;

  int ans = 0;
  
  if ( trad_mode == true ) {
    Cmatrix* ptm = qmats[I_PT_MATRIX];
    ans = bitmap_same_bit_p( ptm->at(x), ptm->at(y) );
  }
  else {
    // We lookup the result in alias matrix
    Cmatrix *am = qmats[I_ALIAS_MATRIX];
    bitmap amx = am->at(x);
    ans = bitmap_bit_p( amx, y );
  }
  
  return ans != 0;
}

int 
BitQS::ListPointsTo( int x, IFilter* filter )
{
  int ans = 0;
  unsigned o;
  bitmap_iterator bi;

  x = pt_map[x];

  if ( x != -1 ) {
    bitmap ptx = qmats[I_PT_MATRIX]->at(x);
    EXECUTE_IF_SET_IN_BITMAP( ptx, 0, o, bi ) {
      ans += iterate_equivalent_set( es2objs[o], filter );
    }
  }

  return ans;
}

int 
BitQS::ListPointedBy( int o, IFilter* filter )
{
  int ans = 0;
  unsigned p;
  bitmap_iterator bi;
  
  o = obj_map[o];

  if ( o != -1 ) {
    bitmap pto = qmats[I_PTED_MATRIX]->at(o);
    EXECUTE_IF_SET_IN_BITMAP( pto, 0, p, bi ) {
      ans += iterate_equivalent_set( es2ptrs[p], filter );
    }
  }

  return ans;
}

int 
BitQS::ListAliases( int x, IFilter* filter )
{
  int ans = 0;
  unsigned q, o;

  // Translate base pointer to its rerepsentative
  x = pt_map[x];

  // Continue if x has points-to information
  if ( x != -1 ) {
    bitmap res = NULL;
    bitmap_iterator bi;

    res = qmats[I_ALIAS_MATRIX]->at(x);
    
    // We compute the result immediately
    if ( res == NULL ) {
      Cmatrix *ptm = qmats[I_PT_MATRIX];
      bitmap ptx = ptm->at(x);
      res = BITMAP_ALLOC(NULL);
      
      Cmatrix *ptedm = qmats[I_PTED_MATRIX];
      EXECUTE_IF_SET_IN_BITMAP( ptx, 0, o, bi ) {
    	bitmap_ior_into( res, ptedm->at(o) );
      }
      
      qmats[I_ALIAS_MATRIX]->set(x, res);
    }
    
    // Extract the base pointers as the answer
    EXECUTE_IF_SET_IN_BITMAP( res, 0, q, bi ) {
      ans += iterate_equivalent_set( es2ptrs[q], filter );
    }
  }
  
  return ans;
}

int 
BitQS::ListModRefVars( int x, IFilter* filter )
{
  int ans = 0;
  unsigned o;
  bitmap_iterator bi;

  x = pt_map[x];

  if ( x != -1 ) {
    bool is_load = ( x >= n_st );
    if ( is_load ) x = x - n_st;
    Cmatrix *mat = qmats[is_load ? I_LOAD_MATRIX : I_STORE_MATRIX];
    bitmap accx = mat->at(x);
    EXECUTE_IF_SET_IN_BITMAP( accx, 0, o, bi ) {
      if ( filter->validate(o) )
	++ans;
    }
  }

  return ans;
}

int
BitQS::ListLoads( int x, IFilter* filter )
{
  int ans = 0;
  unsigned v;
  bitmap_iterator bi;
  bitmap res_other = NULL;
  
  // load-store or store-load conflicts
  res_other = qmats[I_LD_ST_MATRIX]->at(x);

  // Compute and cache it on demand
  if ( res_other == NULL ) {
    bitmap ptx = qmats[I_LOAD_MATRIX]->at(x);
    res_other = BITMAP_ALLOC(NULL);

    Cmatrix* store_matrix = qmats[I_STORE_MATRIX];
    int size = store_matrix->n_r_reps;
    for ( int i = 0; i < size; ++i )
      if ( bitmap_same_bit_p( ptx, store_matrix->at(i) ) != 0 )
	bitmap_set_bit( res_other, i );
    
    qmats[I_LD_ST_MATRIX]->set(x, res_other);
  }

  // visit
  EXECUTE_IF_SET_IN_BITMAP( res_other, 0, v, bi ) {
    ans += iterate_equivalent_set( es2objs[v], filter );
  }
  
  return ans;
}

int
BitQS::ListStores( int x, IFilter* filter )
{
  int ans = 0;
  unsigned v;
  bitmap_iterator bi;
  bitmap res_other = NULL, res_self = NULL;
  bitmap ptx = qmats[I_STORE_MATRIX]->at(x);
  
  // Store-load conflicts
  res_other = qmats[I_ST_LD_MATRIX]->at(x);  
  if ( res_other == NULL ) {
    res_other = BITMAP_ALLOC(NULL);

    Cmatrix* load_matrix = qmats[I_LOAD_MATRIX];
    int size = load_matrix->n_r_reps;
    for ( int i = 0; i < size; ++i )
      if ( bitmap_same_bit_p( ptx, load_matrix->at(i) ) != 0 )
    	bitmap_set_bit( res_other, i );

    qmats[I_ST_LD_MATRIX]->set(x, res_other);
  }
  
  EXECUTE_IF_SET_IN_BITMAP( res_other, 0, v, bi ) {
    ans += iterate_equivalent_set( es2objs[v+n_st], filter );
  }

  // Store-store conflicts
  res_self = qmats[I_ST_ST_MATRIX]->at(x);
  if ( res_self == NULL ) {
    res_self = BITMAP_ALLOC(NULL);
    
    Cmatrix *store_trans_matrix = qmats[I_STORE_TRANS_MATRIX];
    EXECUTE_IF_SET_IN_BITMAP( ptx, 0, v, bi ) {
      bitmap_ior_into( res_self, store_trans_matrix->at(v) );
    }
        
    qmats[I_ST_ST_MATRIX]->set( x, res_self );
  }
  
  EXECUTE_IF_SET_IN_BITMAP( res_self, 0, v, bi ) {
    ans += iterate_equivalent_set( es2objs[v], filter );
  }

  return ans;
}

int
BitQS::ListConflicts( int x, IFilter* filter )
{
  int xx = pt_map[x];
  if ( xx == -1 ) return 0;

  int ans = 0;

  if ( xx < n_st ) { 
    // Store
    ans = ListStores( x, filter );
  }
  else {
    // Load
    ans = ListLoads( x, filter );
  }

  return ans;
}


/*
// We rebuild the index and compre it to the loaded index
bool
BitQS::sanity_check( const char* check_file )
{
  FILE *fp;
  
  fp = fopen( check_file, "r" );
  if ( fp == NULL ) {
    fprintf( stderr, "Loading verify file failed.\n" );
    return false;
  }

  BitIndexer *indexer = parse_points_to_input( fp, 
					       INPUT_START_BY_SIZE );
  fclose( fp );
  
  indexer->fp_generate_index( indexer, true );
  
  // line by line compare
  Cmatrix **mat_set_query = bitqs->qmats;
  Cmatrix **mat_set_index = indexer->imats;
  bool ret = true;

  if ( matrix_equal_p( mat_set_query[I_PT_MATRIX],
		       mat_set_index[I_PT_MATRIX] ) == true ) {
    
    if ( matrix_equal_p( mat_set_query[I_ALIAS_MATRIX],
			 mat_set_index[I_ALIAS_MATRIX] ) == true ) {
      fprintf( stderr, "Verify successfully.\n" );
    }
    else {
      fprintf( stderr, "Verify alias matrix failed.\n" );
      ret = false;
    }
  }
  else {
    fprintf( stderr, "Verify points-to matrix failed.\n" );
    ret = false;
  }

  delete indexer;

  return ret;
}
*/


IQuery*
load_bitmap_index( FILE* fp, int index_type, bool t_mode )
{
  __init_matrix_lib();
  
  // Load header info
  int n, m;
  fread( &n, sizeof(int), 1, fp );
  fread( &m, sizeof(int), 1, fp );
  
  // Process the index body
  BitQS *bitqs = new BitQS(n, m, index_type, t_mode);
  fprintf( stderr, "----------Index File Info----------\n" );

  if ( index_type == PT_MATRIX )
    bitqs->load_pt_index( fp );
  else
    bitqs->load_se_index( fp );

  bitqs->rebuild_eq_groups();
  return bitqs;
}
