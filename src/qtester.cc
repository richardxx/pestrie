/*
 * A test driver of querying system.
 *
 */

#include "query.hh"
#include <unistd.h>
#include <cstdlib>
#include <cstdio>

using namespace std;

// Accepts all pointers
class AllAcceptFilter 
  : public IFilter
{
public:
  bool validate(int x) { return true; }
};

// Only accepts the base pointers
class BasePtrFilter :
  public IFilter
{
public:
  bool validate(int x)
  {
    return valid_ptrs.find(x) != valid_ptrs.end();
  }
  
  void add_ptr(int x)
  {
    valid_ptrs.insert(x);
  }

private:
  set<int> valid_ptrs;
};

static const char* query_strs[] = {
  "random",
  "IsAlias",
  "ListPointsTo",
  "ListPontedTo",
  "ListAliases",
  "ListModRef",
  "ListStores" 
};

struct QueryOpts
{
  int query_type;
  bool print_answers;
  bool trad_mode;
  bool demand_merging;
  const char* input_file;
  const char* query_plan;

  QueryOpts()
  {
    query_type = IS_ALIAS;
    print_answers = false;
    trad_mode = true;
    demand_merging = false;
    input_file = NULL;
    query_plan = NULL;
  }
}
query_opts;

// Program options
static void print_help( const char* prog_name )
{
  printf( "Usage : %s [options] input_file [query_plan]\n", prog_name );
  printf( "Options  : \n" );
  printf( "-p       : Print answer to the queries to stdout.\n" );
  printf( "-t [num] : Specify the query type:\n" );
  printf( "    0    : randomly choose any of the following queries (default);\n" );
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

  while ( (c = getopt( argc, argv, "dpst:h" ) ) != -1 ) {
    switch ( c ) {
    case 'd':
      query_opts.demand_merging = true;
      break;
      
    case 'p':
      query_opts.print_answers = true;
      break;

    case 's':
      query_opts.trad_mode = true;
      break;

    case 't':
      {
	int query_type = std::atoi( optarg );
	if ( query_type < 0 || query_type > 9 )
	  query_type = QT_RANDOM;
	query_opts.query_type = query_type;
      }
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

  query_opts.input_file = argv[optind];
  
  ++optind;
  if ( optind < argc ) {
    query_opts.query_plan = argv[optind];
  }
  
  return true;
}

void
execute_query_plan( IQuery *qs )
{
  int x, y;

  FILE *fp = fopen( query_plan, "r" );
  if ( fp == NULL ) {
    fprintf( stderr, "Cannot open the query plan file. Simulation exits.\n" );
    return;
  }

  // Read the base pointers for query evaluation
  // We also aggregate the base pointers by their representatives
  VECTOR(int) pointers;
  BasePtrFilter* ptr_filter = new BasePtrFilter;
  
  while ( fscanf( fp, "%d", &x ) != EOF ) {
    pointers.push_back( x );
    int es = qs->getPtrEqID(x);
    if ( es != -1 )
      ptr_filter->add_ptr(x);
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
	  bool ans = qs->IsAlias( x, y );
	  if ( print_answers )
	    printf( "(%d, %d) : %s\n", x, y, ans == true ? "true" : "false" );
	}
      }
      break;
      
    case LIST_POINTS_TO:
      {
	x = pointers[i];
	int ans = qs->ListPointsTo( x, ptr_filter );
	if ( print_answers )
	  printf( "%d : %d\n", x, ans );
      }
      break;
      
    case LIST_ALIASES:
      {
	int ans = qs->ListAliases( i, ptr_filters );
	if ( print_answers )
	  printf( "%d : %d\n", x, ans );
      }
      break;
    }
  }

  delete ptr_filters;
}

// We generate random pointers for evaluation
void 
traverse_result( IQuery *qs )
{
  int x, y;
  int ans = 0;

  int index_type = qs->get_query_type();

  if ( (index_type == PT_MATRIX && query_type >= LIST_ACC_VARS) ||
       (index_type == SE_MATRIX && query_type < LIST_ACC_VARS) ) {
    fprintf( stderr, "The query command is not supported by the input index file.\n" );
    return;
  }

  int n = qs->n_of_ptrs();
  int m = qs->n_of_objs();
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

  AllAcceptFilter* ptr_filter = new AllAcceptFilter;

  for ( int i = 0; i < n_query; ++i ) {
    //x = queries[i];
    x = i;

    switch ( query_type ) {
    case IS_ALIAS:
      y = n_query - x;
      ans += (qs->IsAlias( x, y ) == true ? 1 : 0);
      break;
      
    case LIST_POINTS_TO:
      ans += qs->ListPointsTo( x, ptr_filter );
      break;
      
    case LIST_POINTED_TO:
      ans += qs->ListPointedTo( x, ptr_filter );
      break;
      
    case LIST_ALIASES:
      ans += qs->ListAliases( x, ptr_filter );
      break;
  
    case LIST_ACC_VARS:
      ans += qs->ListModRefVars( x, ptr_filter );
      break;
      
    case LIST_CONFLICTS:
      ans += qs->ListConflicts( x, ptr_filter );
      break;
    }
  }
  
  fprintf( stderr, "\nReference answer = %d\n", ans );
  delete ptr_filer;
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
    qs = load_bitmap_index( fp, PT_MATRIX, query_opts.trad_mode );
  else if ( strcmp( magic_code, BITMAP_SE_1 ) == 0 )
    qs = load_bitmap_index( fp, SE_MATRIX, query_opts.trad_mode );
  else if ( strcmp( magic_code, PESTRIE_PT_1 ) == 0)
    qs = load_pestrie_index( fp, PT_MATRIX, query_opts.demand_merging );
  else if ( strcmp( magic_code, PESTRIE_SE_1 ) == 0 )x
    qs = load_pestrie_index( fp, SE_MATRIX, query_opts.demand_merging );

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

  IQuery *qs = load_index();
  if ( qs == NULL ) return -1;

  query_plan != NULL ? 
    execute_query_plan(qs) : traverse_result(qs);
  
  delete qs;

  char buf[128];
  sprintf( buf, "%s querying (%s)", query_strs[query_type],
	   trad_mode == true ? "on-demand" : "use-index" );
  show_res_use( buf );

  return 0;
}

