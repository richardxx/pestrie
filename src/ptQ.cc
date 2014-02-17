/*
 * The querying simulation system for the sparse bitmap points-to index.
 * by richardxx, 2011.9
 */

#include <cstdio>
#include <cstring>
#include <vector>
#include "options.h"
#include "profile_helper.h"
#include "bitmap.h"

using namespace std;

int *es_map, *re_num, *re_num_rev, es_num;
bitmap *ptm, *ptm_rev, *alm;
vector<int> query_points;


bool read_index()
{
  char magic_code[8];
  int open_mode;
  FILE *fp;

  fp = fopen( input_file, "rb" );
  if ( fp == NULL )
    return false;

  fread( magic_code, sizeof(char), 4, fp );
  magic_code[4] = 0;
  open_mode = 0;

  if ( strcmp( magic_code, "PTS1" ) == 0 )
    open_mode = 1;
  else if ( strcmp( magic_code, "PTS2" ) == 0 )
    open_mode = 2;
  
  if ( open_mode != 1 && open_mode != 2 ) {
    fprintf( stderr, "The input file is not valid.\n" );
    return false;
  }
  
  fread( &n, sizeof(int), 1, fp );
  fread( &m, sizeof(int), 1, fp );
  es_map = new int[n];
  re_num = new int[n];
  re_num_rev = new int[n];

  fread( es_map, sizeof(int), n, fp );
  es_num = 0;

  for ( int i = 0; i < n; ++i )
    if ( es_map[i] == i ) {
      re_num[i] = es_num;
      re_num_rev[es_num] = i;
      ++es_num;
    }
  for ( int i = 0; i < n; ++i )
    if ( es_map[i] != i )
      re_num[i] = re_num[ es_map[i] ];

  // Now we read in the points-to and alias matrix
  bitmap_obstack_initialize(NULL);
  ptm = new bitmap[es_num];
  ptm_rev = new bitmap[m];
  alm = new bitmap[es_num];
  for ( int i = 0; i < es_num; ++i ) {
    ptm[i] = BITMAP_ALLOC(NULL);
    alm[i] = BITMAP_ALLOC(NULL);
  }
  for ( int i = 0; i < m; ++i )
    ptm_rev[i] = BITMAP_ALLOC(NULL);


  if ( open_mode == 1 ) {
    
  }
  else if ( open_mode == 2 ) {
    // The points-to matrix
    for ( int i = 0; i < es_num; ++i ) {
      int k, o;

      fread( &k, sizeof(int), 1, fp );
      while (k--) {
	fread( &o, sizeof(int), 1, fp );
	bitmap_set_bit( ptm[i], o );
	bitmap_set_bit( ptm_rev[o], i );
      }
    }
    
    // The alias matrix
    for ( int i = 0; i < es_num; ++i ) {
      int k, q;
      fread( &k, sizeof(int), 1, fp );
      while (k--) {
	fread( &q, sizeof(int), 1, fp );
	bitmap_set_bit( alm[i], q );
	bitmap_set_bit( alm[q], i );
      }
    }
  }

  show_res_use( "Input" );
  return true;
}

bool IsAlias( int x, int y )
{
  x = re_num[x];
  y = re_num[y];
  
  return x != y ? bitmap_same_bit_p( ptm[x], ptm[y] ) : true;
}

bool IsPointsTo( int x, int o )
{
  x = re_num[x];
  return bitmap_bit_p( ptm[x], o );
}

void ListPointsTo( int x, vector<int>& pts )
{
  unsigned o;
  bitmap_iterator bi;

  x = re_num[x];
  EXECUTE_IF_SET_IN_BITMAP( ptm[x], 0, o, bi ) {
    pts.push_back(o);
  }
}

void ListPointedTo( int o, vector<int>& ans )
{
  unsigned p;
  bitmap_iterator bi;

  EXECUTE_IF_SET_IN_BITMAP( ptm_rev[o], 0, p, bi ) {
    ans.push_back( p );
  }
}

void ListAllAliasTo( int x, vector<int>& alias_to )
{
  unsigned q;
  bitmap_iterator bi;

  x = re_num[x];
  EXECUTE_IF_SET_IN_BITMAP( alm[x], 0, q, bi ) {
    alias_to.push_back( re_num_rev[q] );
  }
}

void ListAliasOr( vector<int>& ans )
{
  bitmap temp;
  unsigned p, q;
  bitmap_iterator bi;

  temp = BITMAP_ALLOC(NULL);
  for ( int i = 0; i < query_points.size(); ++i ) {
    p = query_points[i];
    p = re_num[p];
    bitmap_ior_into( temp, alm[p] );
  }

  EXECUTE_IF_SET_IN_BITMAP( temp, 0, q, bi ) {
    ans.push_back( re_num_rev[q] );
  }

  BITMAP_FREE( temp );
}

void ListAliasAnd( vector<int>& ans )
{
  bitmap temp;
  unsigned p, q;
  bitmap_iterator bi;

  temp = BITMAP_ALLOC(NULL);
  p = query_points[0];
  p = re_num[p];
  bitmap_copy( temp, alm[p] );
  
  for ( int i = 1; i < query_points.size(); ++i ) {
    p = query_points[i];
    p = re_num[p];
    bitmap_and_into( temp, alm[p] );
  }
  
  EXECUTE_IF_SET_IN_BITMAP( temp, 0, q, bi ) {
    ans.push_back( re_num_rev[q] );
  }
  
  BITMAP_FREE( temp );
}

void simulate_queries()
{
  int x, y;
  vector<int> res;

  for ( int i = 0; i < n_rand_query; ++i ) {
    switch ( query_type ) {
    case 0:
      x = rand() % n;
      y = rand() % n;
      IsAlias( x, y );
      break;

    case 1:
      x = rand() % n;
      y = rand() % m;
      IsPointsTo( x, y );
      break;

    case 2:
      x = rand() %n;
      res.clear();
      ListPointsTo( x, res );
      break;

    case 3:
      x = rand() % n;
      res.clear();
      ListAllAliasTo( x, res );
      break;

    case 4:
      y = rand() % m;
      res.clear();
      ListPointedTo( y, res );
      break;

    case 5:
      query_points.clear();
      res.clear();
      for ( int j = 0; j < 100; ++j )
	query_points.push_back( rand() % n );
      ListAliasOr( res );
      break;
      
    case 6:
      query_points.clear();
      res.clear();
      for ( int j = 0; j < 5; ++j )
	query_points.push_back( rand() % n );
      ListAliasAnd( res );
      break;
    }
  }

  show_res_use( "Querying" );
}

int main( int argc, char** argv )
{
  if ( parse_options_querier( argc, argv ) == 0 )
    return -1;

  if ( read_index() == false )
    return -1;

  simulate_queries();
  return 0;
}
