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
insert_rectangle( SegTree *p_seg, Point* r, bool is_point )
{
  SegTreeNode **unitRoots = p_seg->unitRoots;

  // We binary search for the right position
  int s = 0, e = p_seg->maxN;
  int mid;
  SegTreeNode *p;
  int x1, x2;
  
  x1 = r->x1;
  if ( is_point ) x2 = x1; 
  else x2 = ((Rectangle*)r)->x2;

  while ( e >= s ) {
    mid = (s+e) / 2;
    p = unitRoots[mid];
    if ( p == NULL ) {
      p = new SegTreeNode;
      unitRoots[mid] = p;
    }

    // hit
    if ( x1 <= mid && mid <= x2 ) {
      if ( is_point )
	p->points = insert_treap( p->points, r );
      else
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

// A universal traversal procedure for customized tree visit
void visit_segTree( SegTree *seg_tree, SEG_TREE_VISITOR fp_visitor, void* par )
{
  SegTreeNode **unitRoots = seg_tree->unitRoots;
  int maxN = seg_tree->maxN;
  
  for ( int i = 0; i < maxN; ++i )
    if ( unitRoots[i] != NULL )
      fp_visitor( unitRoots[i], par );
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

  Rectangle *r;
  Point *pt;
  SegTreeNode *p;

  while ( e >= s ) {
    mid = (s+e) / 2;
    p = unitRoots[mid];
    if ( p == NULL ) break; 

    // We first search the rectangle set, the cast does not fail
    r = (Rectangle*)find_treap( p->rects, y );
    if ( r &&
	 r->x1 <= x && x <= r->x2 && y <= r->y2 ) return true;
    
    if ( x == mid ) {
      // We then search the points set
      pt = find_treap( p->points, y );
      return pt && pt->y1 == y;
    }
    
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
  if ( r.x1 == r.x2 && r.y1 == r.y2 ) {
    Point *pt = new Point( r );
    insert_rectangle( seg_tree, pt, true );
    //(seg_tree->n_points)++;
  }
  else {
    Rectangle *rr = new Rectangle( r );
    insert_rectangle( seg_tree, rr, false );
    /*
    if (r.x1 == r.x2) (seg_tree->n_vertis)++;
    else if (r.y1 == r.y2) (seg_tree->n_horizs)++;
    else (seg_tree->n_rects)++;
    */
  }

  seg_tree->n_pairs += (r.x2-r.x1+1) * (r.y2-r.y1+1);
  //fprintf( test_file, "(%d, %d, %d, %d)\n", r.x1, r.x2, r.y1, r.y2 );
}
