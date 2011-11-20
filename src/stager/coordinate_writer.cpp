/** \file
This file contains a class which contains routines to write the in-memory-stored trajectory to a file. This feature is mainly used for consistency check and visualization of the artifical motions.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/
 
// direct header
#include "stager/coordinate_writer.hpp"

// standard header
#include <fstream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>

// special library headers
#include <boost/lexical_cast.hpp>

// other headers
#include "control.hpp"
#include "log.hpp"


using namespace std;

DCDCoordinateWriter::DCDCoordinateWriter(std::string file,size_t blocks,size_t entries) :
   	blocks_(blocks),
	entries_(entries),
	file_(file)
{
    
}

void DCDCoordinateWriter::init() {    
    std::ofstream out(file_.c_str(),ios_base::binary | ios_base::ate);
    
    int32_t marker = 0; // fortran 4 byte marker		
    char* p_marker = (char*) (&marker);

	// headsize:
    marker= 4*21;
	out.write(p_marker,sizeof(int32_t)); // FORTRAN Marker / head size 1
    std::string ident("CORD");
    out.write(ident.c_str(),sizeof(int32_t)); // MAGIC IDENTIFIER 2
    int32_t nof = blocks_;
    out.write( (char*) (&nof),sizeof(int32_t)); // number of frames 3
    marker= 0; out.write(p_marker,sizeof(int32_t)); // dummy 4
    int32_t one = 1;
	out.write((char*)(&one),sizeof(int32_t)); // timestep between frames 5

	out.write(p_marker,sizeof(int32_t)); // timesteps in simulation
	out.write(p_marker,sizeof(int32_t)); // nstep or istart, nsavc
	out.write(p_marker,sizeof(int32_t)); // unassigned 8
	out.write(p_marker,sizeof(int32_t)); // unassigned 9
	out.write(p_marker,sizeof(int32_t)); // unassigned 10
	out.write(p_marker,sizeof(int32_t)); // unassigned 11

    float dt = 1;
	out.write((char*) (&dt),sizeof(float)); // size of timestep 12
    int32_t marker_false = 0;
    int32_t marker_true = 1;
	out.write((char*)(&marker_true),sizeof(int32_t)); // optional block 1, deactivate 14
	
	out.write((char*)(&marker_false),sizeof(int32_t)); // optional block 2, deactivate 15
	out.write(p_marker,sizeof(int32_t)); // unassigned 16
	out.write(p_marker,sizeof(int32_t)); // unassigned 17
	out.write(p_marker,sizeof(int32_t)); // unassigned 18
	out.write(p_marker,sizeof(int32_t)); // unassigned 19
	out.write(p_marker,sizeof(int32_t)); // unassigned 20
	out.write(p_marker,sizeof(int32_t)); // unassigned 20
	out.write(p_marker,sizeof(int32_t)); // unassigned 20
	
    marker=24; // !=0 indicates CHARMM format
	out.write(p_marker,sizeof(int32_t)); // unassigned 21

    // final header marker
    marker= 4*21;
	out.write(p_marker,sizeof(int32_t)); // unassigned
    marker= 4;
	out.write(p_marker,sizeof(int32_t)); // unassigned
    marker= 0;
	out.write(p_marker,sizeof(int32_t)); // unassigned

    marker= 4;
	out.write(p_marker,sizeof(int32_t)); // unassigned
    marker= 4;
	out.write(p_marker,sizeof(int32_t)); // unassigned
    marker= entries_; 
	out.write(p_marker,sizeof(int32_t)); // number of atoms
    marker= 4;
	out.write(p_marker,sizeof(int32_t)); // unassigned
    
    out.close();
}

void DCDCoordinateWriter::prepare() {
    std::ifstream in(file_.c_str());
    in.seekg(0,ios_base::end);
    data_offset_= in.tellg();
}

void DCDCoordinateWriter::write(coor_t* data,size_t blockoffset, size_t myblocks) {
    // the data comes in in XYZ XYZ format
    // we need to restructure it for DCD
    
    std::fstream out;
    out.open(file_.c_str(),ios_base::in | ios_base::out | ios_base::binary);
    size_t blockbytesize = 2*sizeof(int32_t)+ 6*sizeof(double) + 3*2*sizeof(int32_t)+ 3*entries_*sizeof(float);
        
    std::streamoff byteoffset = data_offset_ + blockoffset*blockbytesize;
    out.seekp(byteoffset);
    
    int32_t marker = 0; // fortran 4 byte marker		
    char* p_marker = (char*)(&marker);
    
    // need some workspace:
    float* p_entrydatabuffer = (float*) malloc(entries_*sizeof(float));
    // write by blocks
    for(size_t i = 0; i < myblocks; ++i)
    {
        marker = 6*sizeof(double);
		double emptycell[6]={0.0,0.0,0.0,0.0,0.0,0.0};
		char* p_emptycell = (char*) &(emptycell[0]);
        out.write(p_marker,sizeof(int32_t));
		out.write(p_emptycell,6*sizeof(double));
        out.write(p_marker,sizeof(int32_t));	
        // k iterates coordinates X,Y,Z
        for(size_t k = 0; k < 3; ++k)
        {
            // prepare data
            for(size_t j = 0; j < entries_; ++j) p_entrydatabuffer[j]=data[i*entries_*3 + j*3 + k];
            // write data
            marker = entries_*sizeof(float);
            out.write(p_marker,sizeof(int32_t));
            out.write((char*)p_entrydatabuffer, entries_*sizeof(float));
            out.write(p_marker,sizeof(int32_t));
        }        
    }
    delete p_entrydatabuffer;
    out.close();
}