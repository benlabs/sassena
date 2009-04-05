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
#include "frame.hpp"

class DcdHeader {	
public:
	int32_t headsize;
	int32_t fingerprint;
	int32_t number_of_frames;
	int32_t dummy1;
	int32_t timesteps_between_frames;
	char buf1[24];
	float size_of_timestep;
	int32_t flag_ext_block1;
	int32_t flag_ext_block2;
	// size1 == size2
};

// byte_pos of frame N:
// init_byte_pos + block_size_byte*(N-1)
class DcdFramesetDescriptor {
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
	long frame_number_offset;

	std::string filename;
	long number_of_frames;
	long number_of_atoms;
	long flag_ext_block1;
	long flag_ext_block2;
	std::streamoff init_byte_pos;
	std::streamoff block_size_byte;
	std::streamoff x_byte_offset;
	std::streamoff y_byte_offset;
	std::streamoff z_byte_offset;
	std::streamoff block1_byte_offset;
	std::streamoff block2_byte_offset;
};

// need to install serialize for std::streampos and std::streamoff
namespace boost {
namespace serialization {

//template<class Archive> 
//void serialize(Archive & ar, std::fstream::off_type & g, const unsigned int version)
//{
//
//}

} // namespace serialization
} // namespace boost

//


class DcdFrames {
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        ar & total_frames;
        ar & framesets;
    }
	/////////////////// 
	
private:
	public:
	long total_frames;
	// this list helps mapping the absolute frame numbers to the file index
	std::vector<DcdFramesetDescriptor> framesets;

	// this is a helper function. It delivers the right descriptor for a given framenumber
	bool find_FramesetDescriptor(int framenumber,DcdFramesetDescriptor& dcd_desc);


	DcdFrames() :  total_frames(0) {}
	DcdFrames(const std::string filename,Atoms& atoms) :  total_frames(0) { add_file(filename,atoms); }
	~DcdFrames() {}

	void add_file(const std::string filename,Atoms& atoms,int recursion_trigger = 5);
	int number_of_framesets() { return framesets.size(); }
	int number_of_frames() { return total_frames; }

	void read(int frame_number,DcdFrame& blank_frame);
};


#endif
