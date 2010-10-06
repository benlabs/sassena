/*
 *  smath.hpp
 *
 *  Created on: May 26, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
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

template<class T> void square_elements(std::vector<std::complex<T> >& data);
template<class T> void add_elements(std::vector<std::complex<T> >& target,const std::vector<std::complex<T> >& source);
template<class T> void auto_correlate_fftw(std::vector<std::complex<T> >& data);
template<class T> void auto_correlate_direct(std::vector<std::complex<T> >& data);
template<class T> std::complex<T> reduce(const std::vector<std::complex<T> >& data);

}

#endif

// end of file
