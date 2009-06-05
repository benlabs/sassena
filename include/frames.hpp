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
#include <boost/serialization/map.hpp>
#include <boost/serialization/serialization.hpp>

// other headers
#include "atoms.hpp"
#include "atomselection.hpp"
#include "frame.hpp"

// Forward:
class Frameset;
class DCDFrameset;
class PDBFrameset;
class TRRFrameset;
class XTCFrameset;

////////////////////////////////////////////////////////////////////////////////
// wrapper class to 'store' frames in
////////////////////////////////////////////////////////////////////////////////

class Frames {
	friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & framecache;
		ar & framecache_max;
		ar & currentframe_i;
		ar & number_of_frames;
		ar & wrapping;
		ar & centergroup_selection;
		// list all possible derived classes for Frameset
        ar.register_type(static_cast<DCDFrameset*>(NULL));
        ar.register_type(static_cast<PDBFrameset*>(NULL));
        ar.register_type(static_cast<XTCFrameset*>(NULL));
        ar.register_type(static_cast<TRRFrameset*>(NULL));
		ar & framesets;
    }
	///////////////////
	
//	std::vector<FrameFilePosLocator> framefileposlocators;	
	std::vector<Frameset*> framesets;
	
	std::map<size_t,Frame> framecache;

	size_t framecache_max;
	size_t currentframe_i;

	size_t number_of_frames;

	// maps a global framenumber to a specific framseet
	Frameset& find_frameset(size_t framenumber);

	// translate global framenumber into frameset specific one
	size_t scope_framenumber(size_t framenumber);
	void test_framenumber(size_t framenumber);
	
public:
	// unit cell behaviour:
	bool wrapping;
	Atomselection centergroup_selection;

	Frames() : framecache_max(2), currentframe_i(-1), number_of_frames(0) {}
	~Frames();
	
	size_t size();

	// push a frameset, frameset is specialized
	size_t add_frameset(const std::string filename,const std::string filetype,Atoms& atoms);

	// load a frame into the framecache, set as current
	void load(size_t framenumber,Atoms& atoms,std::map<std::string,Atomselection>& atomselections);
	void load(size_t framenumber,Atoms& atoms,Atomselection& atomselection);
	
	// retrieve the current frame
	Frame& current();
	
//	void write(std::string filename, std::string af = "pdb") { atoms.write(filename,currentframe(),af); }	
	
	void clear_cache();
};


////////////////////////////////////////////////////////////////////////////////
// abstraction class. real frameset has to be format dependent
////////////////////////////////////////////////////////////////////////////////

class Frameset {
	friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & frame_number_offset;
		ar & number_of_frames;
    }
public:
	virtual ~Frameset() {}
	
	size_t frame_number_offset;
	size_t number_of_frames;
	
	// derived classes have to support this!
	virtual void read_frame(size_t internalframenumber,Frame& cf) {}
};

////////////////////////////////////////////////////////////////////////////////
// specialized framesets, format dependent
////////////////////////////////////////////////////////////////////////////////

// DCD Frameset and dependents

class DCDHeader {	
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

class DCDFrameset : public Frameset {
	friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & boost::serialization::base_object<Frameset>(*this);
		ar & filename;
		ar & init_byte_pos;
		ar & number_of_atoms;
		ar & flag_ext_block1;
		ar & flag_ext_block2;
		ar & block_size_byte;
		ar & x_byte_offset;
		ar & y_byte_offset;
		ar & z_byte_offset;
		ar & block1_byte_offset;
		ar & block2_byte_offset;
    }

	std::string filename;	
	// seek position where data starts within file
	std::streamoff init_byte_pos;
	
	long number_of_atoms;
	unsigned long flag_ext_block1;
	unsigned long flag_ext_block2;
	std::streamoff block_size_byte;
	std::streamoff x_byte_offset;
	std::streamoff y_byte_offset;
	std::streamoff z_byte_offset;
	std::streamoff block1_byte_offset;
	std::streamoff block2_byte_offset;

	bool detect(const std::string filename);	
public:
	
	// allow construction w/o reading file -> call init manually
	DCDFrameset() {}
	void init(std::string filename,size_t framenumber_offset);

	// this constructor should be called by default
	DCDFrameset(std::string filename,size_t frame_number_offset) { init(filename,frame_number_offset); }
	
	// internalframenumber used for positioning file pointer, data loaded into Frame argument
	void read_frame(size_t internalframenumber,Frame& cf);
	
};


// PDB Frameset and dependents

class PDBFrameset : public Frameset {
	friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & boost::serialization::base_object<Frameset>(*this);
		ar & filename;
		ar & number_of_atoms;
    }

	std::string filename;
	long number_of_atoms;
	
	bool detect(const std::string filename);	
public:
	
	// allow construction w/o reading file -> call init manually
	PDBFrameset() {}
	void init(std::string filename,size_t framenumber_offset);

	// this constructor should be called by default
	PDBFrameset(std::string filename,size_t frame_number_offset) { init(filename,frame_number_offset); }
	
	// internalframenumber used for positioning file pointer, data loaded into Frame argument
	void read_frame(size_t internalframenumber,Frame& cf);
	
};


// XTC Frameset and dependents

class XTCFrameset : public Frameset {
	friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & boost::serialization::base_object<Frameset>(*this);
		ar & filename;
		ar & number_of_atoms;
		ar & frame_byte_offsets;		
    }

	std::string filename;
	long number_of_atoms;
	
	std::vector<std::ios::streamoff> frame_byte_offsets;	
	
	bool detect(const std::string filename);	
public:
	
	// allow construction w/o reading file -> call init manually
	XTCFrameset() {}
	void init(std::string filename,size_t framenumber_offset);

	// this constructor should be called by default
	XTCFrameset(std::string filename,size_t frame_number_offset) { init(filename,frame_number_offset); }
	
	// internalframenumber used for positioning file pointer, data loaded into Frame argument
	void read_frame(size_t internalframenumber,Frame& cf);
	
};


// TRR Frameset and dependents

class TRRFrameset : public Frameset {
	friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & boost::serialization::base_object<Frameset>(*this);
		ar & filename;
		ar & number_of_atoms;
		ar & frame_byte_offsets;
    }

	std::string filename;
	long number_of_atoms;
	
	std::vector<std::streamoff> frame_byte_offsets;		
	
	bool detect(const std::string filename);	
public:
	
	// allow construction w/o reading file -> call init manually
	TRRFrameset() {}
	void init(std::string filename,size_t framenumber_offset);

	// this constructor should be called by default
	TRRFrameset(std::string filename,size_t frame_number_offset) { init(filename,frame_number_offset); }
	
	// internalframenumber used for positioning file pointer, data loaded into Frame argument
	void read_frame(size_t internalframenumber,Frame& cf);
	
};

#endif
