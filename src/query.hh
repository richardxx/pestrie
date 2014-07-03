/*
 * A header that defines the interface for answering queries.
 *
 * by richardxx, 2011.8
 */

#ifndef QUERY_H
#define QUERY_H

#include <cstdio>
#include "constants.hh"

// Query Types
#define N_QURIES 7
#define QT_RANDOM   0              // We choose one randomly
#define IS_ALIAS        1
#define LIST_POINTS_TO  2
#define LIST_POINTED_TO 3
#define LIST_ALIASES    4
#define LIST_ACC_VARS   5
#define LIST_CONFLICTS  6


// Querying result filter
class IFilter
{
public:
  // Keep x or not?
  virtual bool validate(int x) = 0;
};


// Querying interface
class IQuery
{
public:
  // Test if (x, y) is an alias pair
  virtual bool IsAlias( int x, int y ) = 0;

  // Test if x points to o
  virtual int ListPointsTo( int x, IFilter* filter ) = 0;

  // Show all the pointers in ptrs that are alias to x
  virtual int ListAliases( int x, IFilter* filter ) = 0;

  // Show all the pointers that point to o
  virtual int ListPointedBy( int o, IFilter* filter ) = 0;

  // Show all the variables that modify or read by statement x
  virtual int ListModRefVars( int x, IFilter* filter ) = 0;

  // Show all the statements that read/write conflict to x 
  virtual int ListConflicts( int x, IFilter* filter ) = 0;

public:
  virtual int getPtrEqID(int x);
  virtual int getObjEqID(int x);

public:
  virtual int nOfPtrs() = 0;
  virtual int nOfObjs() = 0;
  virtual int getIndexType() = 0;
};

// Generating the querying instance

extern IQuery* 
load_bitmap_index( std::FILE* fp, int index_type, bool t_mode );

extern IQuery* 
load_pestrie_index( std::FILE* fp, int index_type, bool d_mering );

#endif
