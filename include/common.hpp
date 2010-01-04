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
#include <sstream>
#include <vector>

// special library headers

// other headers

inline std::string to_s(double value) {
	std::stringstream ss;
	ss << value;
	return ss.str();
}

inline std::string to_s(long value) {
	std::stringstream ss;
	ss << value;
	return ss.str();	
}

inline std::string to_s(size_t value) {
	std::stringstream ss;
	ss << value;
	return ss.str();	
}
inline std::string to_s(int value) {
	std::stringstream ss;
	ss << value;
	return ss.str();	
}

std::vector<double> flatten(std::vector<std::complex<double> >& cvalues);

std::vector<std::complex<double> > compress(std::vector<double>& rvalues);

#endif

// end of file
