/*
 *  assignment.hpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef DECOMPOSITION__ASSIGNMENT_HPP_
#define DECOMPOSITION__ASSIGNMENT_HPP_

// common header
#include "common.hpp"

// standard header
#include <fstream>
#include <string>
#include <vector>

// special library headers

// other headers

class Assignment {
protected:
    size_t NN_;
    size_t NAF_;
    size_t rank_;

    size_t size_;
    size_t offset_;

public:
    Assignment(size_t NN, size_t rank,size_t NAF) : NN_(NN), NAF_(NAF), rank_(rank), size_(0), offset_(0) {}
    
    
    virtual size_t operator[](size_t index) = 0 ;
    virtual size_t size() = 0 ;
    virtual size_t offset() = 0 ;
    virtual bool contains(size_t i) = 0;
    virtual size_t index(size_t i) = 0;
};

class ModAssignment : Assignment {

public:
    
    ModAssignment(size_t NN, size_t rank,size_t NAF); 
    
    size_t operator[](size_t index);
    size_t size();
    size_t offset();
    bool contains(size_t i);
    size_t index(size_t i);
};

class DivAssignment : Assignment {

public:
    
    DivAssignment(size_t NN, size_t rank,size_t NAF); 
    
    size_t operator[](size_t index);
    size_t size();
    size_t offset();
    bool contains(size_t i);
    size_t index(size_t i);
};



#endif

// end of file
