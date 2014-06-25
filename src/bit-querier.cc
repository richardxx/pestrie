/*
 * The querying system for the sparse bitmap points-to/side-effect index.
 *
 * by richardxx,
 * initial, 2011.9
 * enhancement, 2012.7
 * refactor, 2012.10
 */

#include <cstdio>
#include <cstring>
#include <vector>
#include <set>
#include "query.hh"
#include "profile_helper.h"
#include "matrix-ops.hh"
#include "kvec.hh"
#include "bit-index.hh"

using namespace std;

class BitQS : public IQuery
{
public:
  // Interface functions
  bool IsAlias( int x, int y);
  int ListPointsTo( int x, const IFilter* filter );
  int ListAliases( int x, const IFilter* filter );
  int ListPointedBy( int o, const IFilter* filter );
  int ListModRefVars( int x, const IFilter* filter );
  int ListConflicts( int x, const IFilter* filter );

public:
  BitQS(int n_ptrs, int n_objs, int type)
  {
    n_of_mat = 0;
    n_ld = n_st = 0;
    n_es = 0;
    pt_map = obj_map = NULL;
    
    n = n_ptrs; m = n_objs;
    index_type = type;

    if ( type == PT_MATRIX ) {
      n_of_mat = N_OF_PT_INDEX;
      obj_map = new int[n_objs];
      es2objs = new VECTOR(int)[n_objs];
    }
    else {
      n_of_mat = N_OF_SE_INDEX;
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
  // Input matrices
  Cmatrix **qmats;
  int n_of_mat;

  int n, m;             // #pointers, #objects
  int n_ld, n_st;       // #loads, #stores
  int n_es;             // #pointer equivalent sets

  // An one-level index to navigate the search the locations of compressed pointers/object
  int *pt_map, *obj_map;
  int index_type;

  // Members for each pointer (object) representative
  VECTOR(int) *es2ptrs = NULL;
  VECTOR(int) *es2objs = NULL;
};

//static int cnt_same_es = 0;

static BitQS*
load_pt_index( FILE *fp )
{
  // Load header info
  int n, m;
  fread( &n, sizeof(int), 1, fp );
  fread( &m, sizeof(int), 1, fp );
  
  // Load mapping info
  BitQS *bitqs = new BitQS(n, m, PT_MATRIX);
  fread( bitqs->pt_map, sizeof(int), n, fp );
  fread( bitqs->obj_map, sizeof(int), m, fp );

  // Load the index
  int i = 0;
  while ( i < N_OF_LOADABLE_PT_INDEX ) {
    int m_type, dim_r, dim_c;

    fread( &m_type, sizeof(int), 1, fp );
    fread( &dim_r, sizeof(int), 1, fp );
    fread( &dim_c, sizeof(int), 1, fp );
    if ( m_type != i ) i = m_type;
    
    Cmatrix *cm = new Cmatrix( dim_r, dim_c, true, false );
    bitmap *mat = cm->mat;
    bool skip = (trad_mode && m_type > I_PT_MATRIX);
    for ( int k = 0; k < dim_r; ++k )
      mat[k] = bitmap_read_row( fp, COMPRESSED_FORMAT, skip );
    
    profile_matrix( cm, pt_matrix_info[i], stderr );
    bitqs->qmats[i++] = cm;
  }
 
  // Compute the pointed-by matrix
  bitqs->qmats[I_PTED_MATRIX] = transpose( bitqs->qmats[I_PT_MATRIX] );

  return bitqs;
}

static BitQS*
load_se_index( FILE *fp )
{
  int i;

  // Load header info
  int n, m;
  fread( &n, sizeof(int), 1, fp );
  fread( &m, sizeof(int), 1, fp );

  // Load mapping info
  BitQS *bitqs = new BitQS(n, m, SE_MATRIX);
  int *pt_map = new int[n];
  fread( pt_map, sizeof(int), n, fp );  

  // Profile the mapping information
  int n_ld = 0, n_st = 0;
  for ( i = 0; i < n; ++i ) {
    if ( pt_map[i] >= n ) ++n_ld;
    else ++n_st;
  }

  bitqs->n_st = n_st;
  bitqs->n_ld = n_ld;

  // Load the index
  i = 0;
  while ( i < N_OF_LOADABLE_SE_INDEX ) {
    int m_type, dim_r, dim_c;

    fread( &m_type, sizeof(int), 1, fp );
    fread( &dim_r, sizeof(int), 1, fp );
    fread( &dim_c, sizeof(int), 1, fp );
    if ( m_type != i ) i = m_type;
    
    Cmatrix *cm = new Cmatrix( dim_r, dim_c, true, false );
    bitmap *mat = cm->mat;
    bool skip = (trad_mode && m_type > I_LOAD_MATRIX);
    for ( int k = 0; k < dim_r; ++k )
      mat[k] = bitmap_read_row( fp, COMPRESSED_FORMAT, skip );

    profile_matrix( cm, se_matrix_info[i], stderr );
    bitqs->qmats[i] = cm;
    ++i;
  }
  
  // Compute the transposed store and load matrix
  bitqs->qmats[I_STORE_TRANS_MATRIX] = transpose( bitqs->qmats[I_STORE_MATRIX] );
  bitqs->qmats[I_LOAD_TRANS_MATRIX] = transpose( bitqs->qmats[I_LOAD_MATRIX] );
  bitqs->qmats[I_LD_ST_MATRIX] = transpose( bitqs->qmats[I_ST_LD_MATRIX] );

  return bitqs;
}

BitQS*
load_bitmap_index( FILE* fp, int index_type )
{
  __init_matrix_lib();
  fprintf( stderr, "\n-------Input: %s-------\n", input_file );

  BitQS *bitqs = NULL;
  if ( index_type == PT_MATRIX )
    bitqs = load_pt_index( fp );
  else
    bitqs = load_se_index( fp );

  // Now we create the group -> pointers mapping
  // The statements in side-effect analysis are treated specially
  int n = bitqs->n;
  int m = bitqs->m;
  int n_st = bitqs->n_st;
  int n_ld = bitqs->n_ld;
  int n_es = -1;

  for ( int i = 0; i < n; ++i ) {
    int es = pt_map[i];
    if ( es == -1 ) continue;

    if ( es >= n ) {
      // A load statement
      // We reallocate its position to make the mapping compact
      es = es - n + n_st;
      bitqs->pt_map[i] = es;
    }

    if ( es > n_es ) n_es = es;
    bitqs->es2ptrs[es].push_back(i);
  }

  if ( index_type == PT_MATRIX ) {
    // Decompress the equivalent objects
    for ( int i = 0; i < m; ++i ) {
      int es = bitqs->obj_map[i];
      if ( es != -1 )
	bitqs->es2objs[es].push_back(i);
    }
  }

  // We set the number of equivalent pointer sets to be the largest es ID + 1
  bitqs->n_es = n_es + 1;
  
  show_res_use( "Index loading" );
  return bitqs;
}

static int
iterate_equivalent_set( VECTOR(int) *es_set, const IFilter* filter )
{
  int ans = 0;
  int size = es_set->size();
  
  for ( int i = 0; i < size; ++i ) {
    int q = es_set->at(i);
    if ( filter->validate(q) )
      ans++;
  }
  
  return ans;
}

bool 
BitQS::IsAlias( int x, int y )
{
  x = pt_map[x];
  if ( x == -1 ) return false;
  y = pt_map[y];
  if ( y == -1 ) return false;

  if ( x == y ) {
    //++cnt_same_es;
    return true;
  }

  int ret = 0;

  if ( trad_mode == true ) {
    Cmatrix* ptm = qmats[I_PT_MATRIX];
    ret = bitmap_same_bit_p( ptm->at(x), ptm->at(y) );
  }
  else {
    // We lookup the result in alias matrix
    Cmatrix *am = qmats[I_ALIAS_MATRIX];
    bitmap amx = am->at(x);
    ret = bitmap_bit_p( amx, y );
  }
  
  return ret != 0;
}

int 
BitQS::ListPointsTo( int x )
{
  int ans = 0;
  unsigned o;
  bitmap_iterator bi;

  x = pt_map[x];

  if ( x != -1 ) {
    bitmap ptx = qmats[I_PT_MATRIX]->at(x);
    EXECUTE_IF_SET_IN_BITMAP( ptx, 0, o, bi ) {
      VECTOR(int) *objs = &es2objs[o];
      ans += iterate_equivalent_set( objs );
    }
  }

  return ans;
}

int 
BitQS::ListPointedTo( int o, const IFilter* filter )
{
  int ans = 0;
  unsigned p;
  bitmap_iterator bi;
  
  o = obj_map[o];

  if ( o != -1 ) {
    bitmap pto = qmats[I_PTED_MATRIX]->at(o);
    EXECUTE_IF_SET_IN_BITMAP( pto, 0, p, bi ) {
      VECTOR(int) *ptrs = &es2pointers[p];
      ans += iterate_equivalent_set( ptrs, filter );
    }
  }

  return ans;
}

int 
BitQS::ListAliases( int x, const IFilter* filter )
{
  int ans = 0;
  unsigned q, o;
  bitmap_iterator bi;
  bitmap res;

  // Translate base pointer to its rerepsentative
  x = pt_map[x];

  // Continue if x has points-to information
  if ( x != -1 ) {
    Cmatrix **qmats = qmats;
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
      VECTOR(int) *ptrs = &es2pointers[q];
      ans += iterate_equivalent_set( ptrs, filter );
    }
  }
  
  return ans;
}
/*
int
ListAliases_by_pointers( BitQS *bitqs, int i, VECTOR(int) &pointers )
{
  int ans = 0;
  unsigned q;
  bitmap_iterator bi;
  bitmap res;

  int x = pointers[i];
  x = bitqs->pt_map[x];

  if ( x != -1 ) {
    Cmatrix **qmats = bitqs->qmats;
    res = qmats[I_ALIAS_MATRIX]->at(x);

    // We compute the result immediately
    if ( res == NULL ) {
      Cmatrix *ptm = qmats[I_PT_MATRIX];
      bitmap ptx = ptm->at(x);
      res = BITMAP_ALLOC(NULL);

      int size = pointers.size();
      for ( int j = i + 1; j < size; ++j ) {
	int y = pointers[j];
	y = bitqs->pt_map[y];
	if ( y == -1 ) continue;
	
	bitmap pty = ptm->at(y);
	if ( bitmap_same_bit_p( ptx, pty ) != 0 )
	  bitmap_set_bit( res, y );
      }
      
      // The bits set in res are not the representatives but the individual pointers
      qmats[I_ALIAS_MATRIX]->set(x, res);
    }
    
    // visit
    EXECUTE_IF_SET_IN_BITMAP( res, 0, q, bi ) {
      if ( q >= 0 )
	++ans;
    }
  }
  
  return ans;
}
*/

int 
BitQS::ListModRefVars( int x, const IFilter* filter )
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
BitQS::ListLoads( int x, const IFilter* filter )
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
    VECTOR(int) *stores = &es2pointers[v];
    ans += iterate_equivalent_set( stores, filter );
  }
  
  return ans;
}

int
BitQS::ListStores( int x, const IFilter* filter )
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
    VECTOR(int) *stmts = &es2pointers[v+n_st];
    ans += iterate_equivalent_set( stmts, filter );
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
    VECTOR(int) *stmts = &es2pointers[v];
    ans += iterate_equivalent_set( stmts, filter );
  }

  return ans;
}

int
BitQS::ListConflicts( int x, const IFilter* filter )
{
  int xx = bitqs->pt_map[x];
  if ( xx == -1 ) return 0;

  int ans = 0;

  if ( xx < bitqs->n_st ) { 
    // Store
    ans = ListStores( x, filter );
  }
  else {
    // Load
    ans = ListLoads( x, filter );
  }

  return ans;
}
