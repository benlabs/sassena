/*
 *  This file is part of the software sassena
 *
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008-2010 Benjamin Lindner
 *
 */

#ifndef SAMPLE__ATOMSELECTION_HPP_
#define SAMPLE__ATOMSELECTION_HPP_

// common header
#include "common.hpp"

// standard header
#include <string>
#include <vector>

// special library headers
#include <boost/serialization/access.hpp>
#include <boost/serialization/list.hpp>

class IAtomselection {
protected:
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {

    }
    
public:
    virtual size_t operator[](size_t index) = 0;    
    virtual size_t size() = 0;
};

class IndexAtomselection : public IAtomselection {
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<IAtomselection, IndexAtomselection>(*this);
        ar & ids_;
    }
	///////////////////
    std::vector<size_t> ids_;    

public:	
    IndexAtomselection() {}
    IndexAtomselection(std::vector<size_t> ids);

    size_t operator[](size_t index);    
    size_t size();	
};

class RangeAtomselection  : public IAtomselection {
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<IAtomselection, RangeAtomselection>(*this);
        ar & from_;
        ar & to_;
    }
	///////////////////
	
    size_t from_;
    size_t to_;
public:	
    RangeAtomselection() : from_(0),to_(0) {}
    RangeAtomselection(size_t from, size_t to);
    
    size_t operator[](size_t index);    
    size_t size();
};

#endif

// end of file
