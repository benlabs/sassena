/*
 *  decompose.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

// direct header
#include "decompose.hpp"

// standard header
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>

// special library headers

// other headers

using namespace std;

EvenDecompose::EvenDecompose(size_t total, size_t slices) {
	if (slices<=0) throw;

	size_t minperslice = total/slices;
	size_t leftover = total % minperslice;
	
	for(size_t i = 0; i < slices; ++i)
	{
		size_t thisslice = minperslice;
		if (leftover>i) thisslice +=1;
		this->push_back(thisslice);
	}
}

vector<size_t> EvenDecompose::indexes_for(size_t pos) {
	vector<size_t> result(this->at(pos));
	size_t offset=0;
	for(size_t i = 0; i < pos; ++i)
	{
		offset += this->at(i);
	}
	size_t maxi = this->at(pos);
	for(size_t i = 0; i < maxi; ++i)
	{
		result[i]=offset+i;
	}
	return result;
}

size_t EvenDecompose::max() {
	// EvenDecompose is deterministic in that the last entry has always the most assigned
	if (this->size()>0) {
		return this->at(0);
	} else {
		return 0;
	}
}

size_t EvenDecompose::min() {
	// EvenDecompose is deterministic in that the last entry has always the least assigned
	if (this->size()>0) {
		return this->at(this->size()-1);
	} else {
		return 0;
	}
}

size_t EvenDecompose::rank_of(size_t index) {
	size_t total = 0;
	size_t rank = 0;
	for(size_t i = 0; i < this->size(); ++i)
	{
		total += this->at(i);
		if (index>total) break;
		rank++;
	}
	
	if ((index<0) || (index>=total)) throw;

	return rank;
}

RModuloDecompose::RModuloDecompose(size_t total, size_t slices) {
	if (slices<=0) throw;

	indexes_for_slice.resize(slices);

	// we have to save indexes for this Decompose
	for(size_t i = 0; i < total; ++i)
	{
		size_t row = i / slices;
		if ((row%2)==0) { 
			size_t thisslice = i % slices;	
			indexes_for_slice[thisslice].push_back(i);		
		} 
		else {
			size_t thisslice = slices - (i % slices) - 1;
			indexes_for_slice[thisslice].push_back(i);
		}
		
	}
}

vector<size_t> RModuloDecompose::indexes_for(size_t pos) {
	return indexes_for_slice[pos];
}
