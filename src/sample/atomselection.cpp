/*
 *  atomselection.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

// direct header
#include "sample/atomselection.hpp"

// standard header
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>

// special library headers
#include <boost/regex.hpp>

// other headers
#include "sample/atom.hpp"
#include "sample/atoms.hpp"
#include "decomposition/decompose.hpp"
#include "control.hpp"
#include "log.hpp"


using namespace std;

//Atomselection& Atomselection::operator+=(const Atomselection& that) {
//    
//    for(size_t i = 0; i < that.indexes.size(); ++i)
//    {
//        add(that.indexes[i]);
//    }
//}
//
//Atomselection& Atomselection::operator-=(const Atomselection& that) {
//    
//    for(size_t i = 0; i < that.indexes.size(); ++i)
//    {
//        remove(that.indexes[i]);
//    }
//}

void Atomselection::add(size_t index) {
	indexes.push_back(index);
}

void Atomselection::remove(size_t index) {
	vector<size_t>::iterator ai = find(indexes.begin(),indexes.end(),index);
	indexes.erase(ai);
}

std::vector<Atomselection> Atomselection::slice(size_t parts) {

    std::vector<Atomselection> result(parts);
    
    size_t total = indexes.size();
        
    for(size_t i = 0; i < total; ++i)
    {
        size_t rank = (i*parts)/(total);
        result[rank].add(indexes[i]);
    }
	
	return result;
}

// end of file
