/** \file
This file contains mathemical routines.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/

#ifndef MATH__SMATH_HPP_
#define MATH__SMATH_HPP_

// common header
#include "common.hpp"

// standard header
#include <complex>
#include <vector>

// special library headers
#include <fftw3.h>

// other headers
#include <control.hpp>
#include <log.hpp>

namespace smath {
    
/** 
Computes the sine within the first oscillation. Trades off speed vs. accuracy.
*/
inline double sine(double& x)
{	
    const double B = 4/M_PI;
    const double C = -4/(M_PI*M_PI);

    double y = B * x + C * x * fabs(x);

    //  const float Q = 0.775;
    const double P = 0.225;

    y = P * (y * abs(y) - y) + y;   // Q * y + P * y * abs(y)
	return y;
}

/** 
Element-wise multiplication of an array
*/
template<class T> void multiply_elements(const T& factor,std::vector<std::complex<T> >& data);

/** 
Element-wise multiplication of an array
*/
void multiply_elements(const double factor,fftw_complex* data,size_t NF);

/** 
Element-wise squaring of an array
*/
template<class T> void square_elements(std::vector<std::complex<T> >& data);

/** 
Element-wise addition of two arrays
*/
template<class T> void add_elements(std::vector<std::complex<T> >& target,const std::vector<std::complex<T> >& source);

/** 
Element-wise addition of two arrays
*/
void add_elements(fftw_complex* target,const fftw_complex* source,size_t NF);

/** 
Replaces the data within the array with its auto-correlated value using FFT (scales with N)
*/
template<class T> void auto_correlate_fftw(std::vector<std::complex<T> >& data,fftw_plan p1,fftw_plan p2);

/** 
Replaces the data within the array with its auto-correlated value using the direct (scales with N*N)
*/
template<class T> void auto_correlate_direct(std::vector<std::complex<T> >& data);

/** 
Computes the sum of an array
*/
template<class T> std::complex<T> reduce(const std::vector<std::complex<T> >& data);

/** 
Computes the sum of an array
*/
template<class T> std::complex<T> reduce(const fftw_complex* data,size_t N);
  
  /** 
Replaces the data within the array with its auto-correlated value using FFT (scales with N)
*/
void auto_correlate_fftw(std::vector<std::complex<double> >& data,fftw_plan p1,fftw_plan p2,fftw_complex* fftw_planspace);

/** 
Replaces the data within the array with its auto-correlated value using FFT (scales with N)
*/
void auto_correlate_fftw(fftw_complex* data,fftw_plan p1,fftw_plan p2,size_t NF);

/** 
Replaces the data within the array with its auto-correlated value using FFT (scales with N)
*/
void auto_correlate_direct(fftw_complex* data,size_t N);

/** 
Element-wise squaring of an array
*/
void square_elements(fftw_complex* data,size_t N);

/** 
Element-wise addition of two arrays
*/
template<class T> void add_elements(std::vector<std::complex<T> >& target,const fftw_complex* source,size_t N);
}

#endif

// end of file
