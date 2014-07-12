// Copyright 2014, Hong Kong University of Science and Technology. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

/*
 * Implementation of segment tree based structure for point locations queries.
 *
 * by Xiao Xiao
 * initial 2012.9
 */
#include <cstdio>
#include <cstring>
#include <vector>
#include "segtree.hh"

using namespace std;

SegTreeNode*
SegTree::get_unit_node(int x)
{
  SegTreeNode *p = unitNodes[x];

  if ( p == NULL ) {
    p = new SegTreeNode;
    unitNodes[x] = p;
  }
  
  return p;
}

// Insert this rectangle into the index
// Construct new segment tree node when necessary
void 
SegTree::insert_rectangle( int x1, int x2, VLine* r )
{
  // We binary search for the right position
  int s = 0, e = maxN;
  int mid;
  SegTreeNode *p;

  while ( e > s ) {
    mid = (s+e) / 2;
    p = get_unit_node(mid);
    
    // hit
    if ( x1 <= mid && mid <= x2 ) {
      p->rects = insert_treap( p->rects, r );
      return;
    }
    
    // Otherwise, follow the link to the next level
    if ( x1 > mid ) {
      // go right
      s = mid + 1;
      //p->right = true;
    }
    else {
      e = mid;
      //p->left = true;
    }
  }

  // Unreachable
  fprintf( stderr, "Error occurred!\n" );
}

void 
SegTree::insert_unit_node(int x, VLine* pv)
{
  SegTreeNode* p = get_unit_node(x);
  p->rects = insert_treap( p->rects, pv );
}
 
bool 
SegTreeNode::collect_figures( VECTOR(VLine*) &fs )
{
  fs.clear();

  // We directly call the treap traversal procedure
  if ( rects != NULL ) {
    inorder_treap( rects, fs );
    return true;
  }
  
  return false;
}

SegTree::SegTree( int mx )
{ 
  maxN = mx;
  unitNodes = new SegTreeNode*[mx];
  memset( unitNodes, 0, sizeof(void*) * mx );

  n_points = n_horizs = n_vertis = n_rects = 0;
  n_out_points = n_out_horizs = n_out_vertis = n_out_rects = 0;
  n_pairs = 0;
}

SegTree::~SegTree()
{
  for ( int i = 0; i < maxN; ++i )
    if ( unitNodes[i] != NULL )
      delete unitNodes[i];
  
  delete[] unitNodes;
}

/*
 * Lookup if the query point is covered by some rectangle.
 */
bool 
SegTree::query_point( int x, int y )
{
  int s = 0, e = maxN;
  int mid;
  VLine *pl;
  SegTreeNode *p;

  while ( e > s ) {
    // mid is current X coordinate for testing
    mid = (s+e) / 2;
    p = unitNodes[mid];
    // Early termination: no rectangles have exercised this top-down path
    if ( p == NULL ) break;

    // We search the figures
    pl = find_treap( p->rects, y );
    if ( pl ) {
      if ( x == mid ) {
	if ( y <= pl->y2 ) return true;
      }
      else {
	// pl->get_type() == SIG_RECT
	Rectangle* r = (Rectangle*)pl;
	if ( r->x1 <= x && x <= r->x2 && y <= r->y2 ) 
	  return true;
      }
    }
    
    if ( x == mid ) 
      return false;
    
    x > mid ? 
      (s = mid + 1) : (e = mid);
  }

  // Last try
  p = unitNodes[x];
  if ( p != NULL ) {
    pl = find_treap( p->rects, y );
    if ( pl && y <= pl->y2 ) return true;
  }

  return false;
}

/*
 * We distinguish the points and rectangles.
 * We also record basic statistical information.
 */
void 
SegTree::insert_segtree( const Rectangle& r )
{
  if ( r.x1 == r.x2 ) {
    // We directly insert this figure
    VLine* p = new VLine(r);
    insert_unit_node( r.x1, p );
    
    if (r.y1 == r.y2) n_points++;
    else n_vertis++;
  }
  else {
    VLine* p = new Rectangle( r );
    insert_rectangle( r.x1, r.x2, p);
  
    if (r.y1 == r.y2) n_horizs++;
    else n_rects++;    
  }

  n_pairs += (r.x2-r.x1+1) * (r.y2-r.y1+1);
}

/*
 * Aligning the figures by their left bounds.
 */
void
SegTree::flush_left_shapes()
{
  // Collect the rectangles and reinsert
  VECTOR(VLine*) fs;
  //fprintf( stderr, "----------Flush left\n" );

  for ( int i = 0; i < maxN; ++i ) {
    SegTreeNode *segNode = unitNodes[i];
    if ( segNode == NULL ||
	 !segNode->collect_figures( fs ) )
      continue;
    
    // Redistribute the figures
    int size = fs.size();
    //if ( size >= 1000 ) fprintf( stderr, "X = %d, Size = %d\n", i, size );
    for ( int j = 0; j < size; ++j ) {
      VLine* r = fs[j];
      
      if ( r->get_type() == SIG_RECT ) {
	int x1 = ((Rectangle*)r)->x1;
	
	if ( x1 != i ) {
	  // We first remove this figure
	  segNode->rects = remove_treap( segNode->rects, r->y1 );
	
	  // Insert again
	  SegTreeNode* p = get_unit_node(x1);
	  p->rects = insert_treap( p->rects, r );
	}
      }
    }
  }
}

static void
merge_figures(VECTOR(VLine*) &fs)
{
  int size = fs.size();
  if ( size < 2 ) return;

  int last_pos = 0;
  VLine* last_fig = fs[0];
  
  for ( int i = 1; i < size; ++i ) {
    VLine* p = fs[i];
    if ( !last_fig->merge(p) ) {
      fs[last_pos++] = last_fig;
      last_fig = p;
    }
    else {
      delete p;
    }
  }
  
  fs[last_pos++] = last_fig;
  while ( fs.size() > last_pos ) fs.pop_back();
}


// Traverse and write the figures into a binary format file
int
SegTree::dump_figures( FILE* fp )
{
  VECTOR(VLine*) fs;
  int *labels = new int[maxN*3];
  int total_labels = 0;

  //fprintf( stderr, "-------------->Dump figures\n" );
  // Visit the figures
  for ( int i = 0; i < maxN; ++i ) {
    // Reset
    int offset = 1;
    SegTreeNode *segNode = unitNodes[i];
    if ( segNode != NULL &&
	 segNode->collect_figures( fs ) ) {

      // Merging adjacent figures
      merge_figures(fs);

      // Write out
      int size = fs.size();
      //if ( size >= 100 ) fprintf( stderr, "X = %d, Size = %d\n", i, size );
      //int lasty = -1;
      for ( int j = 0; j < size; ++j ) {
	VLine *r = fs[j];
	
	// Verify
	/*
	if ( r->y1 <= lasty ) {
	  fprintf( stderr, "Not sorted! X = %d, Y = %d\n", i, r->y1 );
	}
	lasty = r->y2;
	
	if ( r->get_type() == SIG_RECT &&
	     ((Rectangle*)r)->x1 != i ) {
	  fprintf( stderr, "Not aligned! X = %d, Expected = %d\n",
		   ((Rectangle*)r)->x1, i );
	}
	*/
	if ( r->get_type() == SIG_VERTICAL ) {
	  if ( r->y1 == r->y2 ) ++n_out_points;
	  else ++n_out_vertis;
	}
	else {
	  if ( r->y1 == r->y2 ) ++n_out_horizs;
	  else ++n_out_rects;
	}

	offset += r->prepare_labels( labels + offset );
      }
    }
    
    labels[0] = offset - 1;
    total_labels += offset;
    fwrite( labels, sizeof(int), offset, fp );
  }

  delete[] labels;
  return total_labels;
}


SegTree* 
build_segtree( int s, int e )
{
  SegTree* seg_tree = new SegTree(e); 
  return seg_tree;
}
