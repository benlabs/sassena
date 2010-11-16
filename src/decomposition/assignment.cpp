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
#include "decomposition/assignment.hpp"

// standard header
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>

// special library headers

// other headers
#include "log/log.hpp"

using namespace std;

DivAssignment::DivAssignment(size_t NN,size_t rank, size_t NAF) :
 Assignment(NN,rank,NAF)
{
    size_ = 0;
    offset_=0;
    bool init=true;
    size_t start = (rank_*NAF_)/NN_;    
    for(size_t i = start; i < NAF_; ++i)
    {
        if ( rank_ == (i*NN_ / NAF_) ) {
            size_++;
            if (init) {
                init=false;
                offset_=i;
            }
        } else {
            if (!init) break;
        }
    }
}

size_t DivAssignment::operator[](size_t index) {

    if (index<size_) {
        return offset_+index;
    } else {
        Err::Inst()->write("Assignment out of bounds, operator[]");
        throw;
    }
    
}

size_t DivAssignment::size() {
    return size_;
}

size_t DivAssignment::offset() {
    return offset_;
}

bool DivAssignment::contains(size_t i) {
    if (i>=offset_) {
        if (i<offset_+size_) {
            return true;
        }
    }
    return false;
}

size_t DivAssignment::index(size_t i) {
    if (!contains(i)) throw;
    
    return (i-offset_);
}

ModAssignment::ModAssignment(size_t NN,size_t rank, size_t NAF) :
 Assignment(NN,rank,NAF)
{
    size_ = 0;
    offset_=0;
    bool init=true;
    for(size_t i = 0; i < NAF_; ++i)
    {
        if ( rank_ == (i % NN_) ) {
            size_++;
            if (init) {
                init=false;
                offset_=i;
            }
        }
    }
}

size_t ModAssignment::operator[](size_t index) {
    if (index<size_) {
        return offset_+index*NN_;        
    } else {
        Err::Inst()->write("Assignment out of bounds, operator[]");
        throw;
    }
}

size_t ModAssignment::size() {
    return size_;
}

size_t ModAssignment::offset() {
    return offset_;
}

bool ModAssignment::contains(size_t i) {
    if (i<offset_) return false;
    
    if (((i-offset_) % NN_) == 0) return true;
    
    return false;
}

size_t ModAssignment::index(size_t i) {
    if (!contains(i)) throw;
    
    return (i-offset_)/NN_;
}