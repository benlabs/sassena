/** \file 
This file contains a class which contains routines to write the in-memory-stored trajectory to a file. This feature is mainly used for consistency check and visualization of the artifical motions.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/


#ifndef STAGER__COORDINATE_WRITER_HPP_
#define STAGER__COORDINATE_WRITER_HPP_

// common header
#include "common.hpp"

// standard header
#include <string>
#include <vector>
#include <fstream>

// other headers

/** 
Interface for writing coordinates staged within the distributed memory into a file.
*/
class ICoordinateWriter {
public:
    // set layout
    virtual void init() = 0;    
    virtual void prepare() = 0;    
    virtual void write(coor_t* data,size_t blockoffset, size_t myblocks) =0 ;
};

/** 
Writes the coordinates staged within the distributed memory as a DCD file
*/
class DCDCoordinateWriter : public ICoordinateWriter {
    // filehandle:
    
    size_t blocks_;
    size_t entries_;
    
    std::string file_;
    std::streamoff data_offset_;
    
public:
    DCDCoordinateWriter(std::string file,size_t blocks,size_t entries); // initializes a file with the trajectory format

    // set layout
    void init();
    void prepare();
    void write(coor_t* data,size_t blockoffset, size_t myblocks);
    
};

#endif

//end of file
