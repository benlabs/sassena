/*
 *  common.hpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef COMMON_HPP_
#define COMMON_HPP_

#ifdef CMAKE_CXX_COMPILER_ENV_VAR
// needed for xdrfile module
#define C_PLUSPLUS
#endif

// standard header
#include <complex>
#include <cstring>
#include <sstream>
#include <fstream>
#include <vector>

// special library headers

// other headers

// define here
#define PRECISION_TYPE_SINGLE

// triggers here

#ifdef PRECISION_TYPE_DOUBLE
typedef double coor_t;
typedef double coor2_t;
#endif
#ifdef PRECISION_TYPE_SINGLE
typedef float coor_t;
typedef double coor2_t;
#endif


#endif
