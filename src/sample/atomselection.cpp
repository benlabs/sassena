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
#include <boost/lexical_cast.hpp>

// other headers
#include "sample/atom.hpp"
#include "sample/atoms.hpp"
#include "decomposition/decompose.hpp"
#include "control.hpp"
#include "log.hpp"


using namespace std;

IndexAtomselection::IndexAtomselection(std::vector<size_t> ids) :
    ids_(ids)
{
    
}

size_t IndexAtomselection::operator[](size_t index) {
    return ids_[index];
}

size_t IndexAtomselection::size() {
    return ids_.size();
}


RangeAtomselection::RangeAtomselection(size_t from,size_t to) : 
    from_(from), 
    to_(to)
{
    if (to_<from_) {
        Err::Inst()->write("Range Atomselection specified is negative range!");
        Err::Inst()->write(string("from = ")+boost::lexical_cast<string>(from_));
        Err::Inst()->write(string("to = ")+boost::lexical_cast<string>(to_));
        throw;
    }
}

size_t RangeAtomselection::operator[](size_t index) {
    return from_+index;
}

size_t RangeAtomselection::size() {
    return (1+to_)-from_;
}

// end of file
