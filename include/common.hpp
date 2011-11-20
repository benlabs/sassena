/** \file
The content of this file is included by any other file within the project. Use it to apply application wide modifications.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/

#ifndef COMMON_HPP_
#define COMMON_HPP_

//#ifdef CMAKE_CXX_COMPILER_ENV_VAR
//// needed for xdrfile module
//#define C_PLUSPLUS
//#endif

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
