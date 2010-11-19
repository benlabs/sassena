/*
 *  This file is part of the software sassena
 *
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008-2010 Benjamin Lindner
 *
 */

#ifndef SAMPLE__FRAMESET_INDEX_HPP_
#define SAMPLE__FRAMESET_INDEX_HPP_

// common header
#include "common.hpp"

// standard header
#include <fstream>
#include <string>
#include <vector>

// special library headers
#include <boost/regex.hpp> 
#include <boost/serialization/access.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/filesystem.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

// other headers
#include "sample/atoms.hpp"
#include "sample/atomselection.hpp"
#include "sample/frame.hpp"

class FramesetIndex : public std::vector<std::streamoff> {
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<std::vector<std::streamoff> >(*this);
		ar & signature;
    }

    std::vector<char> signature;

public:
    static std::vector<char> generate_signature(std::string filename);
    void set_signature(std::vector<char> sig) { signature = sig; }
    std::vector<char> get_signature() { return signature; }
    
    void load(std::string filename);
    void save(std::string filename);
};

#endif
