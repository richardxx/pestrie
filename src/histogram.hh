// Copyright 2014, Hong Kong University of Science and Technology. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

/*
 * A histogram data statistics.
 * By Xiao Xiao
 * initial, 2011.9
 */

#ifndef HISTOGRAM_HH
#define HISTOGRAM_HH

#include <cstdio>
#include <vector>
#include <climits>

using std::printf;
using std::vector;

class histogram
{
public:
  // The scale all the sections
  vector<long> limits;
  // results[i] stores the values in the half-open interval: (limits[i-1], limits[i]], limits[-1] = -inf, limits[n] = inf
  vector<long> results;  
  // The accumulated weights for each interval
  vector<double> weights;
  // The maximum sample value recorded
  double maxValue;
  // #of samples
  int count;   
  // total weights of all samples
  double tot_wts;  
  
public:
  histogram()
  {
    limits.clear();
    results.clear();
    weights.clear();
    count = 0;
    tot_wts = 0.0;
    maxValue = 2.2250738585072014e-308;
  }

  ~histogram()
  {
    limits.clear();
    results.clear();
    weights.clear();
  }
  
public:
  /*
   * We support the dynamic update of the scales.
   * However, it requires that the caller should guarantee the scales are inputted in small first sorted order.
   */
  void push_scale( long scale )
  {
    limits.push_back( scale );
    results.push_back( 0 );
  }

  void push_scales( long scales[], int ns )
  {
    for ( int i = 0; i < ns; ++i ) {
      limits.push_back( scales[i] );
      results.push_back( 0 );
      weights.push_back( 0.0 );
    }
  }
  
  int find_max_nonzero_scale()
  {
    int i;

    for ( i = limits.size() - 1; i > -1; --i )
      if ( results[i] > 0 ) break;
    
    return i;
  }

  // The interval is deicided by val
  void add_sample( double val )
  {
    int i;
    int n = limits.size();

    if ( val > maxValue ) maxValue = val;

    for (i = 0; i < n; ++i )
      if ( val <= limits[i] )
	break;
    
    // This value may be larger than any of the scale
    // We guarantee the last section is there
    if ( i == n && 
	 results.size() == n ) {
      results.push_back( 0 );
      weights.push_back( 0.0 );
    }

    results[i]++;
    weights[i] += val;
    ++count;
    tot_wts += val;
  }

  // Interval and val are different
  void add_sample( int sec_val, double val )
  {
    int i;
    int n = limits.size();

    if ( val > maxValue ) maxValue = val;
    
    for (i = 0; i < n; ++i )
      if ( sec_val <= limits[i] )
	break;
    
    // This value may be larger than any of the scale
    // We guarantee the last section is there
    if ( i == n && 
	 results.size() == n ) {
      results.push_back( 0 );
      weights.push_back( 0.0 );
    }
    
    results[i]++;
    weights[i] += val;
    ++count;
    tot_wts += val;
  }

  
  // The caller should guarantee other and me has the same limits
  void merge( const histogram& other )
  {
    for ( int i = 0; i < results.size(); ++i )
      results[i] += other.results[i];

    count += other.count;
    if ( other.maxValue > maxValue ) maxValue = other.maxValue;
  }

  int get_samples_count()
  {
    return count;
  }

  // We ouput the results for the firstN bars
  void print_result( FILE* out, const char* title, bool accumulate = true, int firstN = -1 )
  {
    long num = 0;
    long maxN;

    maxN = limits.size();
    if ( firstN == -1 || firstN > maxN ) firstN = maxN; 
    
    fprintf( out, "%s\n", title );
    if ( count == 0 ) {
      fprintf( out, "No samples are inserted, no output!\n" );
      return;
    }
    else {
      fprintf( out,  "Samples : %d\n", count );
    }

    for ( int i = 0; i < firstN; ++i ) {
      if ( accumulate )
	num += results[i];
      else
	num = results[i];
      
      if ( i == 0 )
	fprintf( out, "-inf < x <= %d: %d", limits[0], num );
      else
	fprintf( out, "%d < x <= %d: %d", limits[i-1], limits[i], num );
      
      fprintf( out, ", percentage = %.3lf\%\n", (double)num / count * 100 ); 
    }
    
    // Aggreate the other items in one bar
    if ( !accumulate ) num = 0;
    int n = results.size();
    for ( int i = firstN; i < n; ++i ) {
      num += results[i];
    }

    if ( firstN == 0 )
      fprintf( out, "x > -inf: %d", num );
    else
      fprintf( out, "x > %d: %d", limits[firstN-1], num );

    fprintf( out, ", percentage = %.3lf\%\n", (double)num / count * 100 );     
    fprintf( out, "\n" );
  }


  // We ouput the results for the firstN bars
  void print_weights( FILE* out, const char* title, bool accumulate = true, int firstN = -1 )
  {
    double num = 0.0;
    long maxN;

    maxN = limits.size();
    if ( firstN == -1 || firstN > maxN ) firstN = maxN; 
    
    fprintf( out, "\n%s\n", title );
    if ( count == 0 ) {
      fprintf( out, "No samples are inserted, no output!\n" );
      return;
    }
    else {
      fprintf( out,  "Samples : %d\n", count );
    }

    for ( int i = 0; i < firstN; ++i ) {
      if ( accumulate )
	num += weights[i];
      else
	num = weights[i];
      
      if ( i == 0 )
	fprintf( out, "-inf < x <= %d: %.2lf", limits[0], num );
      else
	fprintf( out, "%d < x <= %d: %.2lf", limits[i-1], limits[i], num );
      
      fprintf( out, ", percentage = %.3lf\%\n", num / tot_wts * 100 ); 
    }
    
    // Aggreate the other items in one bar
    if ( !accumulate ) num = 0.0;
    int n = weights.size();
    for ( int i = firstN; i < n; ++i ) {
      num += weights[i];
    }

    if ( firstN == 0 )
      fprintf( out, "x > -inf: %.2lf", num );
    else
      fprintf( out, "x > %d: %.2lf", limits[firstN-1], num );

    fprintf( out, ", percentage = %.3lf\%\n", num / tot_wts * 100 );     
    fprintf( out, "\n" );
  }

  
};

#endif

