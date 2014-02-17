/*
 * An implementation of William Weihl's algorithm presented in his POPL81 paper for equivalent sets detection.
 * Alias matrix is computed by sets intersection, which is fairly fast in practice.
 * Modified by richardxx, 2009.9
 */

#include "options.h"
#include "profile_helper.h"

struct hash_node
{
  int idx;
  struct hash_node *next;
};

struct hash_node **HT;
bitmap *ptm, *alm;     	  // points-to matrix and alias matrix
int *es_map, *renum;      // equivalent set mapping and renumbering for alias matrix
int n_es;
int hash_size;
int edge_num;

void input()
{
  int i, k;
  int t;
  FILE *fp;

  fp = fopen( input_file, "r" );
  if ( fp == NULL ) return;
  
  fscanf( fp, "%d %d", &n, &m );

  ptm = new bitmap[n];
  for ( i = 0; i < n; ++i )
    ptm[i] = BITMAP_ALLOC(NULL);

  if ( input_format == 0 ) {
    for ( i = 0; i < n; ++i ) {
      fscanf( fp, "%d", &k );
      edge_num += k;
      while ( k > 0 ) {
	fscanf( fp, "%d", &t );
	bitmap_set_bit( ptm[i], t );
	k--;
      }
    }
  }
  else {
    for ( i = 0; i < n; ++i ) {
      while (1) {
	fscanf( fp, "%d", &t );
	if ( t == -1 ) break;
	bitmap_set_bit( ptm[i], t );
	edge_num++;
      }
    }
  }
  
  fclose( fp );

  fprintf( stderr, "\n------------Input:%s---------\n", input_file );
  show_res_use( "Input" );
}

void estimate_index_size()
{
  double intsize;
  double mem;
  int labels;
  
  //intsize = (sizeof(int)*8 - __builtin_clz(n)) * 1.0 / 8;
  intsize = sizeof(int);

  // The mapping vector
  mem = n * intsize;
  mem += bitmap_calculate_memory( ptm, n );
  mem += bitmap_calculate_memory( alm, n_es );

  labels = n;
  labels += bitmap_calculate_labels( ptm, n );
  labels += bitmap_calculate_labels( alm, n_es );

  fprintf( stderr, "\n" );
  fprintf( stderr, "Index labels: %d\n", labels );
  fprintf( stderr, "The bitmap compressed index size is roughly : %.0lfKb\n", mem / 1024 );
  fprintf( stderr, "The uncompressed index size is : %.0lfKb\n", labels * intsize / 1024 );
}

// Equivalent Sets Identification
void make_points_to_index()
{
  int i;
  hashval_t hv;
  struct hash_node *p, *q;
  int s;

  es_map = new int[n];
  
  // Determine the suitable hash size
  for ( i = 0; i < N_PRIMES; ++i )
    if ( primes[i] * 5 >= n ) break;
  if ( i == N_PRIMES ) i--;
  hash_size = primes[i];
  
  HT = new hash_node*[hash_size];
  memset( HT, 0, sizeof(void*) * hash_size );
  
  for ( i = 0; i < n; ++i ) {
    hv = bitmap_hash( ptm[i] ) % hash_size;
    p = HT[hv];
    while ( p ) {
      if ( bitmap_equal_p( ptm[i], ptm[ p->idx ] ) ) break;
      p = p -> next;
    }
    
    if ( p ) {
      // Found identical set
      es_map[i] = p -> idx;
      BITMAP_FREE( ptm[i] );
      continue;
    }
    
    es_map[i] = i;
    p = new hash_node();
    p -> idx = i;
    p -> next = HT[hv];
    HT[hv] = p;
    ++n_es;
  }
  
  fprintf( stderr, "\n--------------Indexing------------\n" );
  show_res_use( "Indexing" );
  fprintf( stderr, "\n" );
  fprintf( stderr, "Nodes = %d, Edges = %d, Equivalent Sets = %d\n", n, edge_num, n_es );
}

void make_alias_matrix()
{
  int i, j, k;
  int s, t;
  int *holder;
  int alias_num;

  renum = new int[n];
  holder = new int[n];
  for ( i = k = 0; i < n; ++i )
    if ( es_map[i] == i ) {
      renum[i] = k;
      holder[k] = i;
      ++k;
    }
  
  alm = new bitmap[k];
  alias_num = 0;
  for ( i = 0; i < k; ++i )
    alm[i] = BITMAP_ALLOC(NULL);

  for ( i = 0; i < k; ++i ) {
    s = holder[i];    
    for ( j = i + 1; j < k; ++j ) {
      t = holder[j];
      if ( bitmap_same_bit_p( ptm[s], ptm[t] ) ) {
	bitmap_set_bit( alm[i], j );
	//bitmap_set_bit( alm[j], i );
	++alias_num;
      }
    }
  }
  
  fprintf( stderr, "\n---------------Make Alias Matrix-----------------\n" );
  show_res_use( "Alias Matrix" );
  //fprintf( stderr, "\n" );
  fprintf( stderr, "Alias Pairs: %d\n", alias_num );
}

/*
 * Output format:
 * Points-to and Alias Matrix:
 *
 * N_p (pointers) N_o(objects)
 * Pointer representative map (N_p)
 *
 * First the points-to matrix is given:
 * k(points-to size) o1 o2 ... ok
 *
 * Then, the alias-to matrix:
 * k(alias-to size) p1 p2 p3 ... pk
 */
void externalize_uncompressed_index()
{
  int i, k;
  unsigned j;
  bitmap_iterator bi;
  FILE *fp = NULL;

  if ( output_file == NULL )
    return;

  fp = fopen( output_file, "wb" );
  if ( fp == NULL ) {
    fprintf( stderr, "Cannot write to the file: %s\n", output_file );
    return;
  }

  // Write magic number
  fwrite( "PTS1", sizeof(char), 4, fp );

  // Write N_p
  fwrite( &n, sizeof(int), 1, fp );
  // Write N_o
  fwrite( &m, sizeof(int), 1, fp );
  // Write the pointer representatitives
  fwrite( es_map, sizeof(int), n, fp );

  // Write the points-to matrix
  for ( i = 0; i < n; ++i ) {
    if ( ptm[i] == NULL ) continue;
    // output number of blocks for this bitmap
    k = bitmap_count_bits( ptm[i] );
    fwrite( &k, sizeof(int), 1, fp );
    
    EXECUTE_IF_SET_IN_BITMAP( ptm[i], 0, j, bi ) {
      fwrite( &j, sizeof(int), 1, fp );
    }
  }
  
  // Write the alias-to matrix
  for ( i = 0; i < n_es; ++i ) {
    if ( alm[i] == NULL ) continue;
    
    k = bitmap_count_bits( alm[i] );
    fwrite( &k, sizeof(int), 1, fp );
    
    EXECUTE_IF_SET_IN_BITMAP( alm[i], 0, j, bi ) {
      fwrite( &j, sizeof(int), 1, fp );
    }
  }

  fclose( fp );
}

/*
 * Output format:
 * Points-to and Alias Matrix:
 *
 * N_p (pointers) N_o(objects)
 * Pointer representative map (N_p)
 *
 * First the points-to matrix is given:
 * k(#blocks) b1(indx, data) b2 ... bk
 *
 * Then, the alias-to matrix:
 * k(#blocks) b1 b2 ... bk
 */
void externalize_compressed_index()
{
  int i, k;
  unsigned j;
  bitmap_iterator bi;
  FILE *fp = NULL;

  if ( output_file == NULL )
    return;

  fp = fopen( output_file, "wb" );
  if ( fp == NULL ) {
    fprintf( stderr, "Cannot write to the file: %s\n", output_file );
    return;
  }

  // Write magic number
  fwrite( "PTS2", sizeof(char), 4, fp );

  // Write N_p
  fwrite( &n, sizeof(int), 1, fp );
  // Write N_o
  fwrite( &m, sizeof(int), 1, fp );
  // Write the pointer representatitives
  fwrite( es_map, sizeof(int), n, fp );

  // Write the points-to matrix
  for ( i = 0; i < n; ++i ) {
    if ( ptm[i] == NULL ) continue;
    bitmap_write_out( ptm[i], fp );
  }
  
  // Write the alias-to matrix
  for ( i = 0; i < n_es; ++i ) {
    if ( alm[i] == NULL ) continue;
    bitmap_write_out( alm[i], fp );
  }
  
  fclose( fp );
}

int answer_query( int x, int y )
{
  x = es_map[x];
  y = es_map[y];
  
  return x != y ? bitmap_same_bit_p( ptm[x], ptm[y] ) : 1;
}

int main( int argc, char **argv )
{
  //fprintf( stderr, "B = %d\n",BLOCK_SIZE );
  if ( parse_options( argc, argv ) == 0 ) return -1;
  
  bitmap_obstack_initialize(NULL);
  input();
  make_points_to_index();
  make_alias_matrix();
  estimate_index_size();

  if ( output_format == 0 )
    externalize_compressed_index();
  else
    externalize_uncompressed_index();

  return 0;
}
