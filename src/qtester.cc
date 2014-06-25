/*
 * A test driver of querying system.
 *
 */

#include "query.hh"
#include <unistd.h>
#include <cstdlib>
#include <cstdio>

using namespace std;

static const char* query_strs[] = {
  "random",
  "IsAlias",
  "ListPointsTo",
  "ListPontedTo",
  "ListAliases",
  "ListModRef",
  "ListStores" 
};

static int query_type = IS_ALIAS;
static bool do_profile = false;
static bool print_answers = false;
static bool trad_mode = false;
static bool demand_merging = false;

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
  printf( "-s       : Use only points-to matrix for querying (Bitmap ONLY).\n" );
  printf( "-d       : Merging the figures up-to-root before querying (Pestrie ONLY).\n" );
}

static bool 
parse_options( int argc, char **argv )
{
  int c;

  while ( (c = getopt( argc, argv, "c:dgpst:h" ) ) != -1 ) {
    switch ( c ) {
    case 'c':
      check_file = optarg;
      break;

    case 'd':
      demand_merging = true;
      break;
            
    case 'g':
      do_profile = true;
      break;

    case 'p':
      print_answers = true;
      break;

    case 's':
      trad_mode = true;
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
  VECTOR(int) pointers;
  VECTOR(int) *es2baseptrs = new VECTOR(int)[n_es];
  
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
	
	if ( trad_mode )
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

  // We permute the pointers/objects
  /*
  int *queries = new int[n_query];
  for ( int i = 0; i < n_query; ++i ) queries[i] = i;

  for ( int i = 0; i < n_query - 1; ++i ) {
    int j = queries[i];
    int k = rand() % ( n_query - i ) + i;
    queries[i] = queries[k];
    queries[k] = j;
  }
  */

  for ( int i = 0; i < n_query; ++i ) {
    //x = queries[i];
    x = i;

    switch ( query_type ) {
    case IS_ALIAS:
      y = n_query - x;
      ans += (IsAlias( bitqs, x, y ) == true ? 1 : 0);
      break;
      
    case LIST_POINTS_TO:
      ans += ListPointsTo( bitqs, x );
      break;
      
    case LIST_POINTED_TO:
      ans += ListPointedTo( bitqs, x );
      break;
      
    case LIST_ALIASES:
      ans += ListAliases_by_representatives( bitqs, x, NULL );
      break;
  
    case LIST_ACC_VARS:
      ans += ans += ListModRefVars(bitqs, x);
      break;
      
    case LIST_CONFLICTS:
      ans += ListConflicts(bitqs, x);
      break;
    }
  }
  
  fprintf( stderr, "\nReference answer = %d\n", ans );
  //delete[] queries;
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

IQuery*
load_index()
{
  FILE *fp = fopen( input_file, "rb" );
  if ( fp == NULL ) return NULL;

  // We first validate the index file
  char magic_code[8];
  fread( magic_code, sizeof(char), 4, fp );
  magic_code[4] = 0;
  
  int index_type = UNDEFINED_MATRIX;
  IQuery* qs = NULL;

  if ( strcmp( magic_code, BITMAP_PT_1 ) == 0)
    qs = load_bitmap_index( fp, PT_MATRIX );
  else if ( strcmp( magic_code, BITMAP_SE_1 ) == 0 )
    qs = load_bitmap_index( fp, SE_MATRIX );
  else if ( strcmp( magic_code, PESTRIE_PT_1 ) == 0)
    qs = load_pestrie_index( fp, PT_MATRIX );
  else if ( strcmp( magic_code, PESTRIE_SE_1 ) == 0 )
    qs = load_pestrie_index( fp, SE_MATRIX );

  fclose( fp );

  if ( qs == NULL ) {
    fprintf( stderr, "This is an INVALID index file.\n" );    
    return NULL;
  }
  
  return qs;
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
	   trad_mode == true ? "on-demand" : "use-index" );
  show_res_use( buf );

  return 0;
}

