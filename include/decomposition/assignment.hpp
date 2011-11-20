/** \file
This file contains operational classes which mimic data types and encapsulate assignment information. 

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
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

/** 
Interface which represents a generic assignment of a number of tasks towards a number of nodes
*/
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

/** 
Specific assignment class which places tasks according to modulo logic. 
*/
class ModAssignment : Assignment {

public:
    
    ModAssignment(size_t NN, size_t rank,size_t NAF); 
    
    size_t operator[](size_t index);
    size_t size();
    size_t offset();
	size_t max(); // retreives the maximum size of any possible assignment
    bool contains(size_t i);
    size_t index(size_t i);
};

/** 
Specific assignment class which places tasks according to div logic. 
*/
class DivAssignment : Assignment {

public:
    
    DivAssignment(size_t NN, size_t rank,size_t NAF); 
    
    size_t operator[](size_t index);
    size_t size(); // retrieve the scoped size of the assignment
    size_t offset(); // retrieves the absolute offset
	size_t max(); // retreives the maximum size of any possible assignment
    bool contains(size_t i);
    size_t index(size_t i);
};



#endif

// end of file
