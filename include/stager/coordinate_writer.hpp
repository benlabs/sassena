/*
 *  This file is part of the software sassena
 *
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008-2010 Benjamin Lindner
 *
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


////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////

class ICoordinateWriter {
public:
    // set layout
    virtual void init(size_t blocks,size_t entries) = 0;    
    virtual void prepare() = 0;    
    virtual void write(coor_t* data,size_t blockoffset, size_t myblocks) =0 ;
};

class DCDCoordinateWriter : public ICoordinateWriter {
    // filehandle:
    
    size_t blocks_;
    size_t entries_;
    
    std::string file_;
    std::streamoff data_offset_;
    
public:
    DCDCoordinateWriter(std::string file); // initializes a file with the trajectory format

    // set layout
    void init(size_t blocks,size_t entries);
    void prepare();
    void write(coor_t* data,size_t blockoffset, size_t myblocks);
    
};

#endif

//end of file
