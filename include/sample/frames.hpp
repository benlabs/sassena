/** \file
This file contains a management class for framesets and specifies framesets for different trajectory types.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/

#ifndef SAMPLE__FRAMES_HPP_
#define SAMPLE__FRAMES_HPP_

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
#include "../vendor/xdrfile-1.1.1/include/xdrfile.h"

// other headers
#include "sample/atoms.hpp"
#include "sample/atomselection.hpp"
#include "sample/frame.hpp"
#include "sample/frameset_index.hpp"

// Forward:
class Frameset;
class FileFrameset;
class DCDFrameset;
class PDBFrameset;
class TRRFrameset;
class XTCFrameset;
class CloneFrameset;

/** 
Management class for framesets. 
*/
class Frames {
	friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & number_of_frames;
		// list all possible derived classes for Frameset
        ar.register_type(static_cast<CloneFrameset*>(NULL));
        ar.register_type(static_cast<DCDFrameset*>(NULL));
        ar.register_type(static_cast<PDBFrameset*>(NULL));
        ar.register_type(static_cast<XTCFrameset*>(NULL));
        ar.register_type(static_cast<TRRFrameset*>(NULL));
		ar & framesets;
    }
	///////////////////
	
//	std::vector<FrameFilePosLocator> framefileposlocators;	
	std::vector<Frameset*> framesets;

	size_t number_of_frames;

	// maps a global framenumber to a specific framseet
	Frameset& find_frameset(size_t framenumber);

	void test_framenumber(size_t framenumber);
	
public:
	// unit cell behaviour:

	Frames() : number_of_frames(0) {}
	
	size_t size();
	
	// push a frameset, frameset is specialized
	size_t add_frameset(const std::string filename,const std::string filetype,size_t first, size_t last, bool last_set, size_t stride,const std::string index_filename,bool index_default,size_t clones);

	// load a Frame
	Frame load(size_t framenumber);	
};


/** 
Interface for framesets. Real framesets are based on a particular format.
*/
class Frameset {
	friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & frame_number_offset;
		ar & number_of_frames;
		ar & number_of_atoms;
    }

public:
	size_t frame_number_offset;
	size_t number_of_frames;
	size_t number_of_atoms;

    // derived classes have to support this!
    virtual void read_frame(size_t framenumber,Frame& cf) = 0;
	
};

/** 
Derived frameset which implements file access functionality.
*/
class FileFrameset : public Frameset {
	friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
 		ar & boost::serialization::base_object<Frameset>(*this);
		ar & filename;	
//		ar & first;
//		ar & last;
//		ar & last_set;
//		ar & stride;		
		ar & frameset_index_;
    }

public:
	virtual ~FileFrameset() {}
	
    FramesetIndex frameset_index_;
	
	
	std::string filename;
	
//	size_t first;
//	size_t last;
//	bool last_set;
//	size_t stride;
	
    void trim_index(size_t first,size_t last,bool last_set,size_t stride);
    void load_index(std::string cache_filename);
    void save_index(std::string cache_filename);

    virtual void generate_index() = 0;
};

////////////////////////////////////////////////////////////////////////////////
// specialized framesets, format dependent
////////////////////////////////////////////////////////////////////////////////

/** 
Type class representing the DCD file header structure
*/
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

/** 
Frameset class for the DCD trajectory file format
*/
class DCDFrameset : public FileFrameset {
	friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & boost::serialization::base_object<FileFrameset>(*this);
		ar & init_byte_pos;
		ar & flag_ext_block1;
		ar & flag_ext_block2;
		ar & block_size_byte;
		ar & x_byte_offset;
		ar & y_byte_offset;
		ar & z_byte_offset;
		ar & block1_byte_offset;
		ar & block2_byte_offset;
    }

	// seek position where data starts within file
	std::streamoff init_byte_pos;
	
	unsigned long flag_ext_block1;
	unsigned long flag_ext_block2;
	std::streamoff block_size_byte;
	std::streamoff x_byte_offset;
	std::streamoff y_byte_offset;
	std::streamoff z_byte_offset;
	std::streamoff block1_byte_offset;
	std::streamoff block2_byte_offset;
	
	std::ifstream dcdfile;

	bool detect(const std::string filename);	
public:
	
	// allow construction w/o reading file -> call init manually
	DCDFrameset() {}
	~DCDFrameset();
	void init(std::string filename,size_t framenumber_offset);

	// this constructor should be called by default
	DCDFrameset(std::string filename, size_t frame_number_offset) { init(filename,frame_number_offset); }
	
	// internalframenumber used for positioning file pointer, data loaded into Frame argument
	void read_frame(size_t framenumber,Frame& cf);
	
	// fill frame_offsets. non-seekable files have to be scanned completely!
	void generate_index();
	void close(); // closes the associated filestream
};


/** 
Frameset class for the PDB trajectory file format, based on VMD animate write
*/
class PDBFrameset : public FileFrameset {
	friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & boost::serialization::base_object<FileFrameset>(*this);
        ar & default_uc;
    }

	std::vector<CartesianCoor3D> default_uc; 
	
	bool detect(const std::string filename);	
public:
	
	
	// allow construction w/o reading file -> call init manually
	PDBFrameset() {}
	void init(std::string filename,size_t framenumber_offset);

	// this constructor should be called by default
	PDBFrameset(std::string filename, size_t frame_number_offset) { init(filename,frame_number_offset); }
	
	// internalframenumber used for positioning file pointer, data loaded into Frame argument
	void read_frame(size_t framenumber,Frame& cf);

	// fill frame_offsets. non-seekable files have to be scanned completely!
	void generate_index();
};


/** 
Frameset class for the XTC trajectory file format
*/
class XTCFrameset : public FileFrameset {
	friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & boost::serialization::base_object<FileFrameset>(*this);
    }
	
	bool detect(const std::string filename);
	XDRFILE* p_xdrfile;
public:
	
	// allow construction w/o reading file -> call init manually
	XTCFrameset() : p_xdrfile(NULL) {}
	~XTCFrameset();
	void init(std::string filename,size_t framenumber_offset);	

	// this constructor should be called by default
	XTCFrameset(std::string filename, size_t frame_number_offset) { init(filename,frame_number_offset); }
	
	// internalframenumber used for positioning file pointer, data loaded into Frame argument
	void read_frame(size_t framenumber,Frame& cf);

	// fill frame_offsets. non-seekable files have to be scanned completely!
	void generate_index();
	void close(); // closes the associated filepointer
};

/** 
Frameset class for the TRR trajectory file format
*/
class TRRFrameset : public FileFrameset {
	friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & boost::serialization::base_object<FileFrameset>(*this);
    }

	bool detect(const std::string filename);
	XDRFILE* p_xdrfile;	
public:
	
	// allow construction w/o reading file -> call init manually
	TRRFrameset() : p_xdrfile(NULL) {}
	~TRRFrameset();
	void init(std::string filename,size_t framenumber_offset);

	// this constructor should be called by default
	TRRFrameset(std::string filename, size_t frame_number_offset) { init(filename,frame_number_offset); }
	
	// internalframenumber used for positioning file pointer, data loaded into Frame argument
	void read_frame(size_t framenumber,Frame& cf);
	
	// fill frame_offsets. non-seekable files have to be scanned completely!
	void generate_index();
	void close(); // closes the associated filepointer
};

/** 
Pseudo frameset which is a reference to another one and can be used to efficiently multiply coordinates.
*/
class CloneFrameset : public Frameset {
	friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & boost::serialization::base_object<Frameset>(*this);
        ar & p_originalframeset;
    }
    
    Frameset* p_originalframeset;
    
public:
    CloneFrameset() {}
    CloneFrameset(Frameset* original,size_t nof);
	// internalframenumber used for positioning file pointer, data loaded into Frame argument
	void read_frame(size_t framenumber,Frame& cf);

};

#endif

// end of file
