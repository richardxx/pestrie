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
#include "options.h"
#include "profile_helper.h"
#include "matrix_ops.hh"
#include "kvec.hh"
#include "query.h"
#include "bit_index.hh"

using namespace std;

struct BitQS
{
  Cmatrix **qmats;
  int n_of_mat;
  int n, m;
  int n_ld, n_st;       // #of load, #of store
  int n_es;             // #of ES of pointers

  // An one-level index to navigate the search the locations of compressed pointers/object
  int *pt_map, *obj_map;

  BitQS()
  {
    qmats = NULL;
    n_of_mat = 0;
    n = m = 0;
    n_ld = n_st = 0;
    n_es = 0;
    pt_map = obj_map = NULL;
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
  }
};


static int index_type;
static int cnt_same_es = 0;

// Members for each pointer (object) representative
static kvec_t(int) *es2pointers = NULL;
static kvec_t(int) *es2objs = NULL;


static BitQS*
load_pt_index( FILE *fp )
{
  // load header info
  int n, m;
  fread( &n, sizeof(int), 1, fp );
  fread( &m, sizeof(int), 1, fp );

  int *pt_map = new int[n];
  int *obj_map = new int[m];
  fread( pt_map, sizeof(int), n, fp );
  fread( obj_map, sizeof(int), m, fp );
  
  // create
  BitQS *bitqs = new BitQS;
  bitqs->qmats = new Cmatrix*[N_OF_PT_INDEX];
  bitqs->n_of_mat = N_OF_PT_INDEX;
  bitqs->n = n;
  bitqs->m = m;
  bitqs->pt_map = pt_map;
  bitqs->obj_map = obj_map;

  // load the index
  int i = 0;
  while ( i < N_OF_LOADABLE_PT_INDEX ) {
    int m_type, dim_r, dim_c;

    fread( &m_type, sizeof(int), 1, fp );
    fread( &dim_r, sizeof(int), 1, fp );
    fread( &dim_c, sizeof(int), 1, fp );
    if ( m_type != i ) i = m_type;
    
    Cmatrix *cm = new Cmatrix( dim_r, dim_c, true, false );
    bitmap *mat = cm->mat;
    bool skip = (slow_mode && m_type > I_PT_MATRIX);
    for ( int k = 0; k < dim_r; ++k )
      mat[k] = bitmap_read_row( fp, COMPRESSED_FORMAT, skip );
    
    profile_matrix( cm, pt_matrix_info[i], stderr );
    bitqs->qmats[i++] = cm;
  }
 
  // We compute the pointed-to-by matrix
  bitqs->qmats[I_PTED_MATRIX] = transpose( bitqs->qmats[I_PT_MATRIX] );
 
  return bitqs;
}

static BitQS*
load_se_index( FILE *fp )
{
  int i;

  // load header info
  int n, m;
  fread( &n, sizeof(int), 1, fp );
  fread( &m, sizeof(int), 1, fp );

  // load the mapping information
  int *pt_map = new int[n];
  fread( pt_map, sizeof(int), n, fp );  

  int n_ld = 0, n_st = 0;
  for ( i = 0; i < n; ++i ) {
    if ( pt_map[i] >= n ) ++n_ld;
    else ++n_st;
  }
  
  // create
  BitQS *bitqs = new BitQS;
  Cmatrix **qmats = new Cmatrix*[N_OF_SE_INDEX];
  bitqs->qmats = qmats;
  bitqs->n_of_mat = N_OF_SE_INDEX;
  bitqs->n = n;
  bitqs->m = m;
  bitqs->n_st = n_st;
  bitqs->n_ld = n_ld;
  bitqs->pt_map = pt_map;

  // load the index
  // We set the dimensions of the matrices
  i = 0;
  while ( i < N_OF_LOADABLE_SE_INDEX ) {
    int m_type, dim_r, dim_c;

    fread( &m_type, sizeof(int), 1, fp );
    fread( &dim_r, sizeof(int), 1, fp );
    fread( &dim_c, sizeof(int), 1, fp );
    if ( m_type != i ) i = m_type;
    
    Cmatrix *cm = new Cmatrix( dim_r, dim_c, true, false );
    bitmap *mat = cm->mat;
    bool skip = (slow_mode && m_type > I_LOAD_MATRIX);
    for ( int k = 0; k < dim_r; ++k )
      mat[k] = bitmap_read_row( fp, COMPRESSED_FORMAT, skip );

    profile_matrix( cm, se_matrix_info[i], stderr );
    qmats[i] = cm;
    ++i;
  }
  
  // The transposed store and load matrix
  qmats[I_STORE_TRANS_MATRIX] = transpose( qmats[I_STORE_MATRIX] );
  qmats[I_LOAD_TRANS_MATRIX] = transpose( qmats[I_LOAD_MATRIX] );
  qmats[I_LD_ST_MATRIX] = transpose( qmats[I_ST_LD_MATRIX] );

  return bitqs;
}

BitQS*
load_index()
{
  __init_matrix_lib();

  FILE *fp = fopen( input_file, "rb" );
  if ( fp == NULL ) return NULL;

  // We first validate the index file
  char magic_code[7];
  fread( magic_code, sizeof(char), 4, fp );
  magic_code[4] = 0;
  
  index_type = UNDEFINED_MATRIX;

  if ( strcmp( magic_code, BITMAP_PT_1 ) == 0)
    index_type = PT_MATRIX;
  else if ( strcmp( magic_code, BITMAP_SE_1 ) == 0 )
    index_type = SE_MATRIX;
  
  if ( index_type == UNDEFINED_MATRIX ) {
    fprintf( stderr, "This is an INVALID bitmap index file.\n" );
    return false;
  }

  // Load
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
  int *pt_map = bitqs->pt_map;
  int n_es = -1;
  
  es2pointers = new kvec_t(int)[n];
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
    es2pointers[es].push_back(i);
  }

  if ( index_type == PT_MATRIX ) {
    // Particular for points-to matrix
    int *obj_map = bitqs->obj_map;
    es2objs = new kvec_t(int)[m];

    for ( int i = 0; i < m; ++i ) {
      int es = obj_map[i];
      if ( es != -1 )
	es2objs[es].push_back(i);
    }
  }

  // We set the number of equivalent pointer sets to be the largest es ID + 1
  bitqs->n_es = n_es + 1;

  show_res_use( "Index loading" );
  return bitqs;
}


bool 
IsAlias( BitQS* bitqs, int x, int y )
{
  x = bitqs->pt_map[x];
  if ( x == -1 ) return false;
  y = bitqs->pt_map[y];
  if ( y == -1 ) return false;

  if ( x == y ) {
    ++cnt_same_es;
    return true;
  }

  int ret = 0;

  if ( slow_mode == true ) {
    Cmatrix* ptm = bitqs->qmats[I_PT_MATRIX];
    ret = bitmap_same_bit_p( ptm->at(x), ptm->at(y) );
  }
  else {
    // We lookup the result in alias matrix
    Cmatrix *am = bitqs->qmats[I_ALIAS_MATRIX];
    bitmap amx = am->at(x);
    ret = bitmap_bit_p( amx, y );
  }
  
  return ret != 0;
}

int
iterate_equivalent_set( kvec_t(int) *es_set )
{
  int ans = 0;
  int size = es_set->size();

  for ( int i = 0; i < size; ++i ) {
    int q = es_set->at(i);
    ++ans;
  }

  return ans;
}

bool 
IsPointsTo( BitQS *bitqs, int x, int o )
{
  x = bitqs->pt_map[x];
  if ( x == -1 ) return false;
  o = bitqs->pt_map[o];
  Cmatrix *cm = bitqs->qmats[I_PT_MATRIX];
  return bitmap_bit_p( cm->at(x), o );
}

int 
ListPointsTo( BitQS *bitqs, int x )
{
  int ans = 0;
  unsigned o;
  bitmap_iterator bi;

  x = bitqs->pt_map[x];

  if ( x != -1 ) {
    bitmap ptx = bitqs->qmats[I_PT_MATRIX]->at(x);
    EXECUTE_IF_SET_IN_BITMAP( ptx, 0, o, bi ) {
      kvec_t(int) *objs = &es2objs[o];
      ans += iterate_equivalent_set( objs );
    }
  }

  return ans;
}

int 
ListPointedTo( BitQS *bitqs, int o )
{
  int ans = 0;
  unsigned p;
  bitmap_iterator bi;

  o = bitqs->obj_map[o];
  if ( o != -1 ) {
    bitmap pto = bitqs->qmats[I_PTED_MATRIX]->at(o);
    EXECUTE_IF_SET_IN_BITMAP( pto, 0, p, bi ) {
      kvec_t(int) *ptrs = &es2pointers[p];
      ans += iterate_equivalent_set( ptrs );
    }
  }

  return ans;
}

int 
ListAliases_by_representatives( BitQS *bitqs, int x, kvec_t(int) *es2baseptrs )
{
  int ans = 0;
  unsigned q, o;
  bitmap_iterator bi;
  bitmap res;

  // Translate base pointer to its rerepsentative
  x = bitqs->pt_map[x];

  // Continue if x has points-to information
  if ( x != -1 ) {
    Cmatrix **qmats = bitqs->qmats;
    res = qmats[I_ALIAS_MATRIX]->at(x);

    // We compute the result immediately
    // if ( res == NULL ) {
    //   Cmatrix *ptm = qmats[I_PT_MATRIX];
    //   bitmap ptx = ptm->at(x);
    //   res = BITMAP_ALLOC(NULL);
      
    //   Cmatrix *ptedm = qmats[I_PTED_MATRIX];
    //   EXECUTE_IF_SET_IN_BITMAP( ptx, 0, o, bi ) {
    // 	bitmap_ior_into( res, ptedm->at(o) );
    //   }
      
    //   qmats[I_ALIAS_MATRIX]->set(x, res);
    // }
    
    // Extract the base pointers as the answer
    EXECUTE_IF_SET_IN_BITMAP( res, 0, q, bi ) {
      kvec_t(int) *ptrs = ( es2baseptrs == NULL ? &es2pointers[q] : &es2baseptrs[q] );
      ans += iterate_equivalent_set( ptrs );
    }

    // if ( slow_mode == true )
    //   BITMAP_FREE(res);
  }
  
  return ans;
}

int
ListAliases_by_pointers( BitQS *bitqs, int i, kvec_t(int) &pointers )
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


static int 
ListModRefVars( BitQS *bitqs, int x )
{
  int ans = 0;
  unsigned o;
  bitmap_iterator bi;

  x = bitqs->pt_map[x];

  if ( x != -1 ) {
    bool is_load = ( x >= bitqs->n_st );
    if ( is_load ) x = x - bitqs->n_st;
    Cmatrix *mat = bitqs->qmats[is_load ? I_LOAD_MATRIX : I_STORE_MATRIX];
    bitmap accx = mat->at(x);
    EXECUTE_IF_SET_IN_BITMAP( accx, 0, o, bi ) {
      if ( o >= 0 )
	++ans;
    }
  }

  return ans;
}

static int
ListLoads( BitQS *bitqs, int x )
{
  int ans = 0;
  unsigned v;
  bitmap_iterator bi;
  bitmap res_other = NULL;
  
  // load-store or store-load conflicts
  Cmatrix **qmats = bitqs->qmats;
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
  int n_st = bitqs->n_st;
  EXECUTE_IF_SET_IN_BITMAP( res_other, 0, v, bi ) {
    kvec_t(int) *stores = &es2pointers[v];
    ans += iterate_equivalent_set( stores );
  }
  
  return ans;
}

static int
ListStores( BitQS *bitqs, int x )
{
  int ans = 0;
  unsigned v;
  bitmap_iterator bi;
  bitmap res_other = NULL, res_self = NULL;
  Cmatrix **qmats = bitqs->qmats;
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
  
  int n_st = bitqs->n_st;
  EXECUTE_IF_SET_IN_BITMAP( res_other, 0, v, bi ) {
    kvec_t(int) *stmts = &es2pointers[v+n_st];
    ans += iterate_equivalent_set( stmts );
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
    kvec_t(int) *stmts = &es2pointers[v];
    ans += iterate_equivalent_set( stmts );
  }

  return ans;
}

int
ListConflicts( BitQS *bitqs, int x )
{
  x = bitqs->pt_map[x];
  if ( x == -1 ) return 0;

  int ans = 0;

  if ( x < bitqs->n_st ) { 
    // Store
    ans = ListStores( bitqs, x );
  }
  // else {
  //   // Load
  //   ans = ListLoads( bitqs, x );
  // }

  return ans;
}

void
execute_query_plan( BitQS *bitqs )
{
  int x, y;

  FILE *fp = fopen( query_plan, "r" );
  if ( fp == NULL ) {
    fprintf( stderr, "Cannot open the query plan file. Simulation exits.\n" );
    return;
  }

  // Read the base pointers for query evaluation
  // We also aggregate the base pointers by their representatives
  int n_es = bitqs->n_es;
  kvec_t(int) pointers;
  kvec_t(int) *es2baseptrs = new kvec_t(int)[n_es];
  
  kv_init(int, pointers);
  while ( fscanf( fp, "%d", &x ) != EOF ) {
    pointers.push_back( x );
    int es = bitqs->pt_map[x];
    //if ( es > n_es ) printf( "%d\n", x );
    if ( es != -1 )
      es2baseptrs[es].push_back(x);
  }
  
  fclose( fp );

  int n_query = pointers.size();
  //fprintf( stderr, "Query plan loaded : %d entries.\n", n_query );
  show_res_use( NULL );


  // Execute
  for ( int i = 0; i < n_query; ++i ) {
    switch ( query_type ) {
    case IS_ALIAS:
      {
	x = pointers[i];
	// No difference in PesTrie case
	for ( int j = i + 1; j < n_query; ++j ) {
	  y = pointers[j];
	  bool ans = IsAlias( bitqs, x, y );
	  if ( print_answers )
	    printf( "(%d, %d) : %s\n", x, y, ans == true ? "true" : "false" );
	}
      }
      break;
      
    case LIST_POINTS_TO:
      {
	x = pointers[i];
	int ans = ListPointsTo( bitqs, x );
	if ( print_answers )
	  printf( "%d : %d\n", x, ans );
      }
      break;
      
    case LIST_ALIASES:
      {
	int ans = 0;
	
	if ( slow_mode )
	  ans = ListAliases_by_pointers( bitqs, i, pointers );
	else {
	  x = pointers[i];
	  ans = ListAliases_by_representatives( bitqs, x, es2baseptrs );
	}
	
	if ( print_answers )
	  printf( "%d : %d\n", x, ans );
      }
      break;
    }
  }

  delete[] es2baseptrs;
}

// We generate random pointers for evaluation
void 
traverse_result( BitQS *bitqs )
{
  int x, y;
  int ans = 0;

  if ( (index_type == PT_MATRIX && query_type >= LIST_ACC_VARS) ||
       (index_type == SE_MATRIX && query_type < LIST_ACC_VARS) ) {
    fprintf( stderr, "The query commands are supported by input index file.\n" );
    return;
  }

  int n = bitqs->n;
  int m = bitqs->m;
  int n_query = ( query_type == LIST_POINTED_TO ? m : n );

  for ( int i = 0; i < n_query; ++i ) {   
    switch ( query_type ) {
    case IS_ALIAS:
      x = rand() % n; y = rand() % n;
      ans += IsAlias( bitqs, x, y );
      break;
      
    case LIST_POINTS_TO:
      ans += ListPointsTo( bitqs, i );
      break;
      
    case LIST_POINTED_TO:
      ans += ListPointedTo( bitqs, i );
      break;
      
    case LIST_ALIASES:
      ans += ListAliases_by_representatives( bitqs, i, NULL );
      break;
  
    case LIST_ACC_VARS:
      ans += ans += ListModRefVars(bitqs, i);
      break;
      
    case LIST_CONFLICTS:
      ans += ListConflicts(bitqs, i);
      break;
    }
  }
  
  //fprintf( stderr, "Total Conflicts = %d\n", ans );
}

// We rebuild the index currently
static bool
sanity_check( BitQS *bitqs )
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
  
  indexer->fp_generate_index( indexer );
  
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

int main( int argc, char** argv )
{
  if ( parse_options( argc, argv ) == 0 )
    return -1;

  BitQS *bitqs = load_index();
  if ( bitqs == NULL ) return -1;
  if ( check_file != NULL &&
       sanity_check( bitqs ) == false )
    return -1;

  query_plan != NULL ? execute_query_plan(bitqs) : traverse_result(bitqs);
  
  // if ( query_type == IS_ALIAS ) {
  //   long total = (long)n_query * (n_query -1);
  //   fprintf ( stderr, "Direct answers = %d, total = %lld\n", cnt_same_es, total/2 );
  // }

  delete bitqs;

  char buf[128];
  sprintf( buf, "%s querying (%s)", query_strs[query_type],
	   slow_mode == true ? "on-demand" : "use-index" );
  show_res_use( buf );

  return 0;
}
