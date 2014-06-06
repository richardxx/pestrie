/*
 * The common procedures for manipulating PesTrie.
 * richardxx, 2012.9
 */

#include <vector>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <cassert>
#include <cstring>
#include "segtree.hh"
#include "histogram.hh"
#include "pestrie.hh"
#include "profile_helper.h"
#include "matrix-ops.hh"

using namespace std;

// An internal structure used to collect the index figures
struct FigureSet
{
  vector<Point*> *points;
  vector<Rectangle*> *vertis;
  vector<Rectangle*> *horizs;
  vector<Rectangle*> *rects;

  FigureSet()
  {
    points = new vector<Point*>;
    vertis = new vector<Rectangle*>;
    horizs = new vector<Rectangle*>;
    rects = new vector<Rectangle*>;
  }

  ~FigureSet()
  {
    delete points;
    delete vertis;
    delete horizs;
    delete rects;
  }
};


static 
void treap_collect_rects( TreapNode* p, void* par )
{
  struct FigureSet* figures = (struct FigureSet*)par;
  Rectangle *rect = (Rectangle*)p->data;
  
  if ( rect->x1 == rect->x2 )
    figures->vertis->push_back( rect );
  else if ( rect->y1 == rect->y2 )
    figures->horizs->push_back( rect );
  else
    figures->rects->push_back( rect );
}

static 
void treap_collect_points( TreapNode* p, void* par )
{
  struct FigureSet* figures = (struct FigureSet*)par;
  figures->points->push_back( p->data );
}

static 
void seg_tree_visitor( SegTreeNode* pseg, void* par )
{
  // We directly call the treap traversal procedure
  if ( pseg->rects != NULL )
    visit_treap( pseg->rects, treap_collect_rects, par );
  if ( pseg->points != NULL )
    visit_treap( pseg->points, treap_collect_points, par );
}


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
    /*
    if ( k == 5715 || k == 17686 || k == 5736 || k == 10130 ) {
      fprintf( stderr, "object %d ( permute id = %d), size = %d\n", Vsrc[k], k, 
	       bitmap_count_bits(mat_T[ Vsrc[k] ]) );
      EXECUTE_IF_SET_IN_BITMAP(mat_T[ Vsrc[k] ],0,x,bi) {
	fprintf( stderr, "%d  ", x );
      }
      fprintf( stderr, "\n" );
    }
    */
   
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

  /*
    int n_points = this->seg_tree->n_points;
    int n_vertis = this->seg_tree->n_vertis;
    int n_horizs = this->seg_tree->n_horizs;
    int n_rects = this->seg_tree->n_rects;
    
    assert( n_rects == this->seg_tree->n_rects );
    assert( n_vertis == this->seg_tree->n_vertis );
    assert( n_horizs == this->seg_tree->n_horizs );
    assert( n_points == this->seg_tree->n_points );
  */

  FigureSet *figures = new FigureSet;
  visit_segTree( this->seg_tree, seg_tree_visitor, figures );
  int n_rects = figures->rects->size();
  int n_vertis = figures->vertis->size();
  int n_horizs = figures->horizs->size();
  int n_points = figures->points->size();
  int n_total_stored = n_points + n_vertis + n_horizs + n_rects;
  int n_gen_rects = this->n_gen_rects;

  fprintf( stderr, "We totally generate %d rectangles, %d of them are indexed.\n", 
	   n_gen_rects, n_total_stored );

  fprintf( stderr, "\t %d vertical lines, percentage = %.2lf\%\n", 
	   n_vertis, 
	   (double)(n_vertis) / n_total_stored * 100 );

  fprintf( stderr, "\t %d horizontal lines, percentage = %.2lf\%\n", 
	   n_horizs,
	   (double)(n_horizs) / n_total_stored * 100 );

  fprintf( stderr, "\t %d points, percentage = %.2lf\%\n", 
	   n_points,
	   (double)(n_points) / n_total_stored * 100 );

  //fprintf( stderr, "nlgn roughly equals to : %.0lf\n", n * log(n) );
  
  fprintf( stderr, "Rectangle pairs : %d, on average %.3lf alias pairs per rectangle\n", 
	   this->seg_tree->n_pairs, 
	   (double)(this->seg_tree->n_pairs)/(n_total_stored) );

  // Estimate the index size (almost real)
  double intsize = sizeof(int);
  int labels = 0;

  // We precisely evaluate the storage size for points
  for ( int i = figures->points->size() - 1; i > -1; ) {    
    Point *p = figures->points->at(i);
    int X = p->x1;
    int j = i - 1;

    while ( j > -1 ) {
      p = figures->points->at(j);
      if ( p->x1 != X ) break;
      --j;
    }

    labels += 2 + i - j;
    i = j;
  }
  
  labels += 6 + n + m + n_rects * 4 + (n_vertis + n_horizs) * 3;
  fprintf( stderr, "Index labels : %d\n", labels );
  fprintf( stderr, "The PesTrie index size is : %.0lfKb\n", labels * intsize / 1024 );

  delete figures;
  show_res_use( NULL );
}

void 
PesTrie::profile_additional()
{
  if ( this->pes_opts->profile_in_detail == false )
    return;
  
  fprintf( stderr, "\n--------------Additional Information--------------\n" );

  int n = this->n;
  int m = this->m;
  int cm = this->cm;
  int vn = this->vn;
  int *bl = this->bl;
  int *pes = this->pes;
  int *es_size = this->es_size;
  bitmap* mat_T = this->mat_T;
  int *r_count = this->r_count;
  int *alias_count = this->alias_count;


  // object induced alias frequency + tree size
  if ( 0 && this->index_type == PT_MATRIX ) {
    histogram obj_freq;
    long freq_scales[] = { 10, 100, 1000, 5000 };
    obj_freq.push_scales( freq_scales, 4 );

    histogram tr_sizes;
    long size_scales[] = { 5, 50, 100, 200 };
    tr_sizes.push_scales( size_scales, 4 );

    // We first collect the #of internal pairs
    int *tr_sz = new int[cm];
    memset( tr_sz, 0, sizeof(int) * cm );
    
    for ( int i = cm; i < vn; ++i ) {
      int tr_code = pes[i];
      tr_sz[tr_code]++;
    }
    
    long tot_internals = 0;
    long tot_pairs = 0;

    for ( int i = 0; i < cm; ++i ) {
      int sz = tr_sz[i] + 1;
      
      int n_internal = sz * (sz-1) / 2;
      tot_internals += sz;
      alias_count[i] += n_internal; 
      tot_pairs += alias_count[i];

      obj_freq.add_sample( alias_count[i] );
      tr_sizes.add_sample( n_internal );
    }
  
    delete[] tr_sz;

    fprintf( stderr, "Internal pairs = %lld, All pairs = %lld, percent = %.2lf\n",
	     tot_internals, tot_pairs, (double)tot_internals/tot_pairs );

    obj_freq.print_result( stderr, "Alias Groups Sizes", false );
    obj_freq.print_weights( stderr, "Alias Groups Weighted Sizes", false );
    tr_sizes.print_result( stderr, "Internal Pairs", false );
  }


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
    fprintf( stderr, "Max hub degree is %.1lf. Arithmetic mean is %.1lf. Geometric mean is %.1lf\n", max_wt, ari_avg, geo_avg );
    hubD.print_result( stderr, "Hub degrees Distribution", false );
    pted_sizes.print_result( stderr, "Pointed-to-by Matrix", false );
  }
      
  // We profile the percentage of alias pairs introduced by first k% trees
  if ( 0 ) {
    histogram pt_seen;
    long percents[] = { cm/1000, cm/100, cm/50, cm/10, cm/5 };
    pt_seen.push_scales( percents, 5 );
    
    for ( int i = 0; i < cm; ++i ) {
      pt_seen.add_sample( i, alias_count[i] );
    }
    
    pt_seen.print_weights( stderr, "Alias pairs introduced by the first K trees.." );
  }

  // We profile the cross edges (#cross edges for each subtree)
  if ( 1 ) {
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
  }

  // Output
  //tree_size.print_result( stderr, "PesTrie Tree Size Distribution", false );
  show_res_use( NULL );
}

/*
 * Now we traverse the segment tree to generate the index file.
 * The index file is in binary form and the format is shown below:
 * 
 * Magic Number (4 bytes)
 * N_p(pointer) N_o(object) N_vn (ES) N_r(rectangles) N_v(verticals) N_h(horizontals) N_p(points)
 * preV labels (N_p+N_o)
 * R
 * V
 * H
 * P
 */
void 
PesTrie::externalize_index( FILE* fp, const char* magic_number)
{
  int i, j, k;
    
  int n = this->n;
  int m = this->m;
  int cm = this->cm;
  int vn = this->vn;
  MatrixRow *r_order = this->r_order;
  int *m_rep = this->m_rep;
  int *preV = this->preV;
  int *bl = this->bl;

  // Fill up the preV mapping for input pointers and objects
  int *pre_aux = new int[n+m];
 
  // First are the pointers
  for ( i = 0; i < n; ++i ) {
    int x = bl[i];
    // a pointer may not have a timestamp
    pre_aux[i] = ( x == -1 ? -1: preV[x] );
  }

  // Second are the objects
  int *obj_pos = new int[cm];
  for ( i = 0; i < cm; ++i ) {
    j = r_order[i].id;
    obj_pos[j] = i;
  }

  for ( i = 0; i < m; ++i ) {
    // i'th object is aggregated into the j'th group (pes node)
    j = ( m_rep == NULL ? i : m_rep[i] );
    // The sorted id of group j
    k = ( j == -1 ? -1 : obj_pos[j] );
    // set the timestamp of object i
    pre_aux[i+n] = ( k == -1 ? -1 : preV[k] );   
  }
  
  // Next, we collect the figures
  FigureSet *figures = new FigureSet;
  visit_segTree( this->seg_tree, seg_tree_visitor, figures );
  int n_rects = figures->rects->size();
  int n_vertis = figures->vertis->size();
  int n_horizs = figures->horizs->size();
  int n_points = figures->points->size();
  /*
  assert( figures->rects->size() == this->seg_tree->n_rects );
  assert( figures->vertis->size() == this->seg_tree->n_vertis );
  assert( figures->horizs->size() == this->seg_tree->n_horizs );
  assert( figures->points->size() == this->seg_tree->n_points );
  */

  // Now we start to output the index
  // Write the magic number
  fwrite( magic_number, sizeof(char), 4, fp );

  // Write N_p
  fwrite( &n, sizeof(int), 1, fp );
  // Write N_o
  fwrite( &m, sizeof(int), 1, fp );
  // Write N_vn
  fwrite( &vn, sizeof(int), 1, fp );
  // Write N_r
  fwrite( &(n_rects), sizeof(int), 1, fp );
  // Write N_v
  fwrite( &(n_vertis), sizeof(int), 1, fp );
  // Write N_h
  fwrite( &(n_horizs), sizeof(int), 1, fp );
  // Write N_p
  fwrite( &(n_points), sizeof(int), 1, fp );
  // Write the preV mappings for pointers+objects
  fwrite( pre_aux, sizeof(int), n + m, fp );

  // Write the index figures
  int labels[] = { 0, 0, 0, 0 };

  for ( i = 0; i < n_rects; ++i ) {
    Rectangle *p = figures->rects->at(i);
    labels[0] = p->x1; labels[1] = p->y1; labels[2] = p->x2; labels[3] = p->y2; 
    fwrite( labels, sizeof(int), 4, fp );
    //printf( "%d %d %d %d\n", p->x1, p->x2, p->y1, p->y2 );
  }

  for ( i = 0; i < n_vertis; ++i ) {
    Rectangle *p = figures->vertis->at(i);
    // We store the values y1, x2, y2
    labels[1] = p->y1; labels[2] = p->x2; labels[3] = p->y2; 
    fwrite( labels + 1, sizeof(int), 3, fp );
    //printf( "%d %d %d %d\n", p->x1, p->x2, p->y1, p->y2 );
  }

  for ( i = 0; i < n_horizs; ++i ) {
    Rectangle *p = figures->horizs->at(i);
    // We store the values x1, y1, x2
    labels[0] = p->x1; labels[1] = p->y1; labels[2] = p->x2;
    fwrite( labels, sizeof(int), 3, fp );
    //printf( "%d %d %d %d\n", p->x1, p->x2, p->y1, p->y2 );
  }

  /*
   * We perform a special treatment to points because they are the mass portion of the index.
   * The format for consecutive points are:
   * X n Y1 Y2 Y3 ... Yn
   * The first is the X-coordinate shared by all the following n points.
   */
  for ( i = n_points -1; i > -1; ) {
    Point *p = figures->points->at(i);
    int X = p->x1;
    j = i - 1;

    while ( j > -1 ) {
      p = figures->points->at(j);
      if ( p->x1 != X ) break;
      --j;
    }

    // The shared X
    labels[0] = X;
    // #points that share X
    labels[1] = i - j;
    fwrite( labels, sizeof(int), 2, fp );

    for ( ; i > j; --i ) {
      p = figures->points->at(i);
      int Y = p->y1;
      fwrite( &Y, sizeof(int), 1, fp );
      //printf( "%d %d %d %d\n", X, X, p->y1, p->y1 );
    }
  }

  delete pre_aux;
  delete obj_pos;
  delete figures;
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
  pestrie->profile_pestrie();  
  pestrie->profile_index();
  pestrie->profile_additional();
}
