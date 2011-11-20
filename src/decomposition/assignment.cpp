/** \file
This file contains operational classes which mimic data types and encapsulate assignment information. 

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
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
#include <boost/lexical_cast.hpp>

// other headers
#include "log/log.hpp"

using namespace std;

DivAssignment::DivAssignment(size_t NN,size_t rank, size_t NAF) :
 Assignment(NN,rank,NAF)
{
    size_t first = (rank_*NAF_)/NN_; 
    size_t next = ((rank_+1)*NAF_)/NN_;
	offset_ = (rank_*NAF_)/NN_; 
	size_ =  next - first; 
}

size_t DivAssignment::operator[](size_t index) {

    if (index<size_) {
        return offset_+index;
    } else {
        Err::Inst()->write("Assignment out of bounds, operator[]");
        Err::Inst()->write(string("offset_=")+boost::lexical_cast<string>(offset_));
        Err::Inst()->write(string("size_=")+boost::lexical_cast<string>(size_));
        Err::Inst()->write(string("index=")+boost::lexical_cast<string>(index));
        Err::Inst()->write(string("NN_=")+boost::lexical_cast<string>(NN_));
        Err::Inst()->write(string("NAF_=")+boost::lexical_cast<string>(NAF_));
        throw;
    }
    
}

size_t DivAssignment::size() {
    return size_;
}

size_t DivAssignment::offset() {
    return offset_;
}

size_t DivAssignment::max() {
	size_t max = NAF_/NN_;
	if ( (NAF_%NN_)!=0 ) max+=1;
	return max;
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
    size_ = NAF_/NN_;
	if ( (NAF_%NN_)!=0 ) {
		if ( rank < ( NAF_ - NN_*(NAF_/NN_) ) ) size_+=1;
	}
    offset_=rank;
}

size_t ModAssignment::operator[](size_t index) {
    if (index<size_) {
        return offset_+index*NN_;        
    } else {
        Err::Inst()->write("Assignment out of bounds, operator[]");
        Err::Inst()->write(string("offset_=")+boost::lexical_cast<string>(offset_));
        Err::Inst()->write(string("size_=")+boost::lexical_cast<string>(size_));
        Err::Inst()->write(string("index=")+boost::lexical_cast<string>(index));
        Err::Inst()->write(string("NN_=")+boost::lexical_cast<string>(NN_));
        Err::Inst()->write(string("NAF_=")+boost::lexical_cast<string>(NAF_));
        throw;
    }
}

size_t ModAssignment::size() {
    return size_;
}

size_t ModAssignment::offset() {
    return offset_;
}

size_t ModAssignment::max() {
	size_t max = NAF_/NN_;
	if ( (NAF_%NN_)!=0 ) max+=1;
	return max;
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

// end of file
