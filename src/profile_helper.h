// Copyright 2014, Hong Kong University of Science and Technology. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

/*
 * This file defines procedures for time and memory measurements.
 *
 * by Xiao Xiao
 * initial: 2009.9
 */

#ifndef PROFILE_HELPER_H
#define PROFILE_HELPER_H

extern 
void show_res_use( const char*, double = 1.0 );

extern 
void* my_malloc( int );

#endif
