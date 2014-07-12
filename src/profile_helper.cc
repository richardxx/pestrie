// Copyright 2014, Hong Kong University of Science and Technology. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

/*
 * Obtaining time and memory usage on linux system.
 *
 * by Xiao Xiao.
 */
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <climits>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include "profile_helper.h"

using namespace std;

#define TIME() (getrusage(RUSAGE_SELF,&ru) ? 0 :			\
		(ru.ru_utime.tv_sec * 1000 + ru.ru_utime.tv_usec / 1000.0 +	\
		 ru.ru_stime.tv_sec * 1000 + ru.ru_stime.tv_usec / 1000.0))

static struct rusage ru;
static double last_tick = 0.0f;
static int last_mem = 0;


// Translate page numbers into megabytes
static int ptok ( unsigned pages )
{
  unsigned ps = 0;
  unsigned tmp;
  int size = INT_MAX;

  /* Initialization.  */
  ps = sysconf( _SC_PAGESIZE );

  /* Conversion.  */
  if (pages > (INT_MAX / ps))
    {                           /* Could overflow.  */
      tmp = pages / 1024;       /* Smaller first, */
      size = tmp * ps;          /* then larger.  */
    }
  else
    {                           /* Could underflow.  */
      tmp = pages * ps;         /* Larger first, */
      size = tmp / 1024;        /* then smaller.  */
    }

  return size / 1024;
}

// Pick memory usage from /proc/ system
static int pick_memory()
{
  int i;
  char temp[2048], label[64];
  int mem;
  FILE *fp;

  fp = fopen( "/proc/self/status", "r" );
  if ( fp == NULL ) return -1;

  for (int i = 0; i < 16; i++) { 
    fgets(temp, 1024, fp);
  }
  
  sscanf( temp, "%s %d", label, &mem );
  fclose(fp);

  return mem;
}

void show_res_use( const char* text, double time_scale )
{
  double cur_tick;
  int cur_mem;

  cur_tick = TIME();
  cur_mem = pick_memory(); 
  //cur_mem = ptok( ru.ru_minflt );
  
  if ( text != NULL ) {
    fprintf( stderr, "%s time (memory): %.0lfms (%dKb) \n", 
	     text, 
	     (cur_tick - last_tick) * time_scale,
	     cur_mem - last_mem );
  }
  
  last_tick = cur_tick;
  last_mem = cur_mem;
}

void* my_malloc( int sz )
{
  void *t;

  if ( (t = malloc( sz )) == NULL ) exit(-1);
  return t;
}
