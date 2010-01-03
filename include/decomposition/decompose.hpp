/*
 *  decompose.hpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef DECOMPOSE_HPP_
#define DECOMPOSE_HPP_

// common header
#include "common.hpp"

// standard header
#include <fstream>
#include <string>
#include <vector>

// special library headers

// other headers

// this class is used to evenly distribute the total number evenly among the slices
class EvenDecompose : public std::vector<size_t> {
		
public:
	EvenDecompose(size_t total,size_t slices); 
	
	std::vector<size_t> indexes_for(size_t pos);
	
	size_t max();
	size_t min();
	
	size_t rank_of(size_t index);
};

// this class evenly distributes with a modulo flavor. R stands for reversing, ie. sequence is left->right,right->left,....
class RModuloDecompose : public std::vector<size_t> {
	std::vector<std::vector<size_t> > indexes_for_slice;
public:
	RModuloDecompose(size_t total,size_t slices); 
	
	std::vector<size_t> indexes_for(size_t pos);
};

#endif

// end of file