/*
 *  frames.hpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef FRAMES_HPP_
#define FRAMES_HPP_

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

// other headers
#include "atoms.hpp"
#include "atomselection.hpp"
#include "frame.hpp"

class FrameFilePosLocator {
public:
	// context
	long frame_number_offset;

	// file dependent
	long number_of_frames;
	std::string filename;	
	// seek position where data starts within file
	std::streamoff init_byte_pos;
};

// byte_pos of frame N:
// init_byte_pos + block_size_byte*(N-1)
class DCDFrameFilePosLocator : public FrameFilePosLocator {
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
    	ar &  frame_number_offset;
     	ar &  filename;
     	ar &  number_of_frames;
     	ar &  number_of_atoms;
     	ar &  flag_ext_block1;
     	ar &  flag_ext_block2;
     	ar &  init_byte_pos;
     	ar &  block_size_byte;
     	ar &  x_byte_offset;
     	ar &  y_byte_offset;
     	ar &  z_byte_offset;
     	ar &  block1_byte_offset;
     	ar &  block2_byte_offset;
    }
	///////////////////
public:

	long number_of_atoms;
	long flag_ext_block1;
	long flag_ext_block2;
	std::streamoff block_size_byte;
	std::streamoff x_byte_offset;
	std::streamoff y_byte_offset;
	std::streamoff z_byte_offset;
	std::streamoff block1_byte_offset;
	std::streamoff block2_byte_offset;
};

class Frames {
	friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		// define this
		throw;
    }
	///////////////////
	
	std::vector<FrameFilePosLocator> framefileposlocators;	
	
	std::map<size_t,Frame> framecache;

	size_t framecache_max;
	size_t currentframe_i;

	// find framefileposlocator for given framenumber
	FrameFilePosLocator& find_locator(int framenumber);
public:
	// unit cell behaviour:
	bool wrapping;
	std::string centergroup;


	Frames() : framecache_max(2);
	
	size_t size();
	
	virtual void add_file(const std::string filename,Atoms& atoms);
	void load(int framenumber,std::vector<Atomselection>& atomselections = {});
	Frame& current();
	
//	void write(std::string filename, std::string af = "pdb") { atoms.write(filename,currentframe(),af); }	
	
	void clear_cache();
};

// specializations:
// DCD
// ...

class DCDFrames : public Frames {
	void read_data(int framenumber,Frame& cf);
public:
	void add_file(const std::string filename,Atoms& atoms);
};


#endif
