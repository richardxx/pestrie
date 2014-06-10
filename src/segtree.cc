/*
 * Implementation of segment tree based index structure for point locations queries.
 * richardxx, 2012.9
 */
#include <cstdio>
#include <cstring>
#include <vector>
#include "segtree.hh"

using namespace std;

static int my_count = 0;
static FILE* test_file = NULL;


// Insert this rectangle into the index
// Construct new segment tree node when necessary
static void 
insert_rectangle( SegTree *p_seg, int x1, int x2, VLine* r )
{
  SegTreeNode **unitRoots = p_seg->unitRoots;

  // We binary search for the right position
  int s = 0, e = p_seg->maxN;
  int mid;
  SegTreeNode *p;

  while ( e >= s ) {
    mid = (s+e) / 2;
    p = unitRoots[mid];
    if ( p == NULL ) {
      p = new SegTreeNode;
      unitRoots[mid] = p;
    }

    // hit
    if ( x1 <= mid && mid <= x2 ) {
      p->rects = insert_treap( p->rects, r );
      break;
    }
    
    // Otherwise, follow the link to the next level
    if ( x1 > mid ) {
      // go right
      s = mid + 1;
      //p->right = true;
    }
    else {
      e = mid - 1;
      //p->left = true;
    }
  }
}

static 
void collect_figures_visitor( SegTreeNode* pseg, FigureSet* fs )
{
  // We directly call the treap traversal procedure
  if ( pseg->rects != NULL )
    visit_treap<VLine*>( pseg->rects, fs->rects );
}
//==================================================================

SegTree::SegTree( int mx )
{ 
  maxN = mx;
  unitRoots = new SegTreeNode*[mx+1];
  memset( unitRoots, 0, sizeof(void*) * (mx+1) );

  n_points = 0;
  n_horizs = 0;
  n_vertis = 0;
  n_rects = 0;
  n_pairs = 0;
}

SegTree::~SegTree()
{
  for ( int i = 0; i <= maxN; ++i )
    if ( unitRoots[i] != NULL )
      delete unitRoots[i];
  
  delete[] unitRoots;
}

struct SegTree* 
build_segtree( int s, int e )
{
  SegTree* seg_tree = new SegTree(e);
  //test_file = fopen( "rects.txt", "w" );
 
  return seg_tree;
}

/*
 * Lookup if the query point is covered by some rectangle.
 */
bool 
query_point( SegTree *seg_tree, int x, int y )
{
  int s = 0, e = seg_tree->maxN;
  int mid;
  SegTreeNode **unitRoots = seg_tree->unitRoots;

  VLine *pl;
  SegTreeNode *p;

  while ( e >= s ) {
    mid = (s+e) / 2;
    p = unitRoots[mid];
    if ( p == NULL ) break; 

    // We search the figure set
    pl = find_treap( p->rects, y );
    if ( pl ) {
      int type = pl->get_type();

      if ( type == SIG_VERTICAL &&
	   pl->y2 >= y ) return true;

      if ( type == SIG_RECT ) {
	Rectangle* r = (Rectangle*)pl;
	if ( r->x1 <= x && x <= r->x2 && y <= r->y2 ) 
	  return true;
      }
    }
    
    if ( x == mid ) break;
    
    if ( x > mid ) {
      s = mid + 1;
    }
    else {
      e = mid - 1;
    }
  }
  
  return false;
}

/*
 * We distinguish the points and rectangles.
 * We also record basic statistical information.
 */
void 
insert_segtree_wrapper( SegTree* seg_tree, const Rectangle& r )
{
  VLine* p = NULL;

  if ( r.x1 == r.x2 ) {
    p = new VLine( r );
    if (r.y1 == r.y2) (seg_tree->n_points)++;
    else (seg_tree->n_vertis)++;
  }
  else {
    p = new Rectangle( r );
    if (r.y1 == r.y2) (seg_tree->n_horizs)++;
    else (seg_tree->n_rects)++;    
  }

  insert_rectangle( seg_tree, r.x1, r.x2, p);
  seg_tree->n_pairs += (r.x2-r.x1+1) * (r.y2-r.y1+1);
  //fprintf( test_file, "(%d, %d, %d, %d)\n", r.x1, r.x2, r.y1, r.y2 );
}

// A universal traversal procedure for customized tree visit
int
dump_figures( SegTree *seg_tree, FILE* fp )
{
  SegTreeNode **unitRoots = seg_tree->unitRoots;
  int maxN = seg_tree->maxN;
  FigureSet* fs = new FigureSet;
  int *labels = new int[maxN*3];
  int total_labels = 0;

  // Visit the figures
  for ( int i = 0; i < maxN; ++i ) {
    // Reset
    int offset = 1;

    if ( unitRoots[i] != NULL ) {
      fs->clear();
      collect_figures_visitor( unitRoots[i], fs );
      vector<VLine*> *rects = fs->rects;

      // Write out
      int size = rects->size();
      for ( int j = 0; j < size; ++j ) {
	VLine *r = rects->at(j);
	offset += r->prepare_labels( labels + offset );
      }
    }

    labels[0] = offset - 1;
    total_labels += offset;
    fwrite( labels, sizeof(int), offset, fp );
  }

  delete fs;
  delete[] labels;

  return total_labels;
}

/*
 * Reinsert all the figures to the treap aligned by left point.
 */
void
flush_left_shapes( SegTree *seg_tree )
{
  // Collect the rectangles and reinsert
  int maxN = seg_tree->maxN;
  SegTreeNode **unitRoots = seg_tree->unitRoots;
  FigureSet* fs = new FigureSet;

  for ( int i = 0; i < maxN; ++i ) {
    SegTreeNode *segNode = unitRoots[i];
    if ( segNode == NULL ||
	 segNode->rects == NULL )
      continue;

    // Collect
    vector<VLine*> *rects_set = fs->rects;
    rects_set->clear();
    visit_treap<VLine*>( segNode->rects, rects_set );
    
    // Redistribute the figures again
    int size = rects_set->size();
    for ( int j = 0; j < size; ++j ) {
      VLine* r = rects_set->at(j);
      
      if ( r->get_type() == SIG_RECT ) {
	int x1 = ((Rectangle*)r)->x1;
	
	if ( x1 != i ) {
	  // We first remove this figure
	  segNode->rects = remove_treap( segNode->rects, r->y1 );
	
	  // Insert again
	  if (unitRoots[x1] == NULL ) {
	    unitRoots[x1] = new SegTreeNode;
	  }
	
	  SegTreeNode* p = unitRoots[x1];
	  p->rects = insert_treap( p->rects, r );
	}
      }
    }
  }
  
  delete fs;
}
