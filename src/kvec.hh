/* The MIT License

   Copyright (c) 2008, by Attractive Chaos <attractivechaos@aol.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/*
  An example:

  #include "kvec.h"
  int main() {
  kvec_t(int) array;
  kv_init(array);
  kv_push(int, array, 10); // append
  kv_a(int, array, 20) = 5; // dynamic
  kv_A(array, 20) = 4; // static
  kv_destroy(array);
  return 0;
  }
*/

/*
  2008-09-22 (0.1.0): The initial version.
  2010-04-29 (0.2.0): Better support for C++, by Xiao Xiao
  2014 (0.3): C++ version is significantly enhanced, by Xiao Xiao
*/

#ifndef AC_KVEC_H
#define AC_KVEC_H

#define INITIAL_SIZE 3

#ifdef __cplusplus

#include <cstring>

// Defined using c++ new operator
// Since we don't know the size of the memory pointed by pin, which should be passed in by the user
// sz1, size of pin;
// sz2, size of pout.
template<typename T> 
T* realloc( T* pin, int sz1, int sz2 )
{
  T* pout = NULL;
  
  if ( sz2 != 0 ) {
    pout = new T[sz2];
    if ( pout == NULL ) return NULL;
    if ( sz1 ) std::memcpy( pout, pin, sz1 < sz2 ? sz1 : sz2 );
  }

  delete pin;
  return pout;
}

template<typename T>
class fast_vec_t
{
private:
  size_t n, m;
  T *a;
  
public:
  fast_vec_t()
  {
    n = m = 0;
    a = NULL;
    //resize( INITIAL_SIZE );
  }

  fast_vec_t(int sz)
  {
    n = m = 0;
    a = NULL;
    resize( sz );
  }

  fast_vec_t( fast_vec_t<T> &other )
  {
    n = m = 0;
    a = NULL;
    copy( other );
  }
  
  ~fast_vec_t()
  {
    if ( a != NULL )
      delete[] a;
  }

  void copy( fast_vec_t<T> &other ) 
  {
    int sz = other.size();
    if ( m < sz ) resize(sz);
    
    for ( int i = 0; i < sz; ++i ) {
      T t = other[i];
      a[i] = t;
    }
  }

  // Swap the contents with another vector
  void swap( fast_vec_t<T> &other )
  {
    int nn = other.n;
    int mm = other.m;
    T *aa = other.a;
    
    other.n = n;
    other.m = m;
    other.a = a;
    
    n = nn; m = mm; a = aa;
  }
  
  void resize( int size )
  {
    T* b = new T[size];

    if ( a == NULL ) {
      a = b;
    }
    else {
      int cp_size = m;
      if ( m > size ) { cp_size = size; n = size; }
      memcpy( b, a, sizeof(T) * cp_size );
      delete[] a;
      a = b;
    }
    
    m = size;
  }

  // We don't check the range
  T& operator[](int inx) { return a[inx]; }
  T& at( int inx ) { return a[inx]; }
  T& back() { return a[n-1]; }
  
  void set( int inx, T& v ) { a[inx] = v; }
  void reset_end( int inx ) {
    if ( inx <= m )
      n = inx;
  }

  T* begin() { return a; }
  T* end() { return a + n; }
  int size() { return n; }
  
  void push_back( const T &x )
  {
    if ( n == m ) {
      int s = (m << 1) + 5;
      resize(s);
    }
    
    a[n++] = x;
  }

  T& pop_back() { return a[--n]; }
 
  void add_all( fast_vec_t<T> &other ) 
  {
    int sz = other.size();
    for ( int i = 0; i < sz; ++i ) {
      T e = other[i];
      push_back(e);
    }
  }
  
  void clear() { n = 0; }
};
  
#else
#include <stdlib.h>
#endif

#define kv_roundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#define kv_A(v, i) ((v).a[(i)])
#define kv_pop(v) ((v).a[--(v).n])
#define kv_size(v) ((v).n)
#define kv_max(v) ((v).m)
#define kv_clear(v) ((v).n = 0)
#define kv_begin(v) ((v).a)
#define kv_end(v) ((v).a+(v).n)
 

/* C++ interface */
#ifdef __cplusplus

#define kvec_t(type) fast_vec_t<type>
#define kv_init(type, v) (v.resize(INITIAL_SIZE))
#define kv_destroy(v) ( delete[] (v).a )

#define kv_resize(type, v, s)						\
  ( (v).a = (type*)realloc((v).a, sizeof(type) * (v).m, sizeof(type) * (s)), (v).m = (s) )

#define kv_push(type, v, x) do {		\
    if ((v).n == (v).m) {			\
      int s = ( (v).m ? ((v).m << 1) + 1 : INITIAL_SIZE );	\
      kv_resize( type, v, s );			\
    }						\
    (v).a[(v).n++] = (x);			\
  } while (0)


#else   /* C interfaces */

#define kvec_t(type) struct { size_t n,m; type *a }
#define kv_init(type, v) ((v).n = 0, (v).m = INITIAL_SIZE, (v).a = (type*)malloc( sizeof(type) * INITIAL_SIZE ) )
#define kv_destroy(v) free((v).a)

#define kv_resize(type, v, s)  ((v).m = (s), (v).a = (type*)realloc((v).a, sizeof(type) * (v).m))

#define kv_push(type, v, x) do {				\
    if ((v).n == (v).m) {					\
      (v).m = (v).m? (v).m<<1 : 2;				\
      (v).a = (type*)realloc((v).a, sizeof(type) * (v).m);	\
    }								\
    (v).a[(v).n++] = (x);					\
  } while (0)

#define kv_a(type, v, i) ((v).m <= (size_t)(i)?				\
			  ((v).m = (v).n = (i) + 1, kv_roundup32((v).m), \
			   (v).a = (type*)realloc((v).a, sizeof(type) * (v).m), 0) \
			  : (v).n <= (size_t)(i)? (v).n = (i)		\
			  : 0), (v).a[(i)]

/*
  #define kv_pushp(type, v) (((v).n == (v).m)?				\
  ((v).m = ((v).m? (v).m<<1 : 2),		\
  (v).a = (type*)realloc((v).a, sizeof(type) * (v).m), 0) \
  : 0), ((v).a + ((v).n++))
*/

#endif

#endif
