/*
 *  This file is part of the software sassena
 *
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008-2010 Benjamin Lindner
 *
 */

// direct header
#include "sample/frameset_index.hpp"

// standard header
#include <fstream>


// special library headers
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include "io/xdrfile/xdrfile.h"
#include "io/xdrfile/xdrfile_xtc.h"
#include "io/xdrfile/xdrfile_trr.h"

// other headers
#include "sample/atomselection.hpp"
#include "sample/center_of_mass.hpp"
#include "control.hpp"
#include "log.hpp"


using namespace std;

std::vector<char> FramesetIndex::generate_signature(std::string filename) {
    std::ifstream binin(filename.c_str(),ios::binary);
    std::streamoff buffersizehalf = 1*1024*1024;
    std::vector<char> buffer(2*buffersizehalf);

    binin.seekg(0,ios::end);
    std::streamoff endpos = binin.tellg();
    std::streamoff filesize = endpos;
    binin.seekg(0,ios::beg);
    
    if (filesize>2*buffersizehalf) {
        binin.read(&(buffer[0]),sizeof(char)*buffersizehalf);
        binin.seekg(ios_base::end-buffersizehalf);
        binin.read(&(buffer[buffersizehalf]),sizeof(char)*buffersizehalf);        
    } else {
        buffer.resize(filesize);
        binin.read(&(buffer[0]),sizeof(char)*filesize);
    } 
    return buffer;
}

void FramesetIndex::load(std::string filename) { 
    std::ifstream in(filename.c_str()); 
    boost::archive::text_iarchive ar(in); 
    ar >> *this; 
}

void FramesetIndex::save(std::string filename) { 
    std::ofstream out(filename.c_str()); 
    boost::archive::text_oarchive ar(out); 
    ar << *this; 
}

void FramesetIndex::select(std::vector<size_t>& selection) {
	std::vector<std::streamoff> copy;
	
	for(size_t i = 0; i < selection.size(); ++i) copy.push_back(selection[i]);
	this->clear();
	for(size_t i = 0; i < copy.size(); ++i)	this->push_back(copy[i]);	
	
}

// end of file
