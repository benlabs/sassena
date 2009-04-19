/*
 *  frames.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

// direct header
#include "frames.hpp"

// standard header
#include <fstream>

// special library headers

// other headers
#include "atomselection.hpp"
#include "settings.hpp"

using namespace std;


void Frames::test_framenumber() {
	if ((currentframe_i>=size()) || (currentframe_i<0)) {
		cerr << "ERROR>> " << " currentframe integer out of bound" << endl;
	}
}

Frames::~Frames() {
	for (size_t i=0;i<framesets.size();i++) {
		delete framesets[i];
	}	
}

void Frames::add_frameset(const std::string filename,const std::string filetype,Atoms& atoms) {
	// get current total number of frames -> offset
	// set offset
	
	// Detect frame type
	if (filetype=="dcd") {
		// create an empty frameset and then initialize w/ filename
		Frameset* fsp = new DCDFrameset(filename,number_of_frames);
		framesets.push_back(fsp);
		// increase frame counter
		number_of_frames += fsp->number_of_frames;
	}
	else {
		cerr << "ERROR>> filetype '" << filetype << "' not supported." << endl;
		throw;
	}
}

size_t Frames::size() {
	return number_of_frames;
}

Frameset& Frames::find_frameset(size_t framenumber) {
	test_framenumber();
	for (size_t i=0;i<framesets.size();i++) {
		if (framesets[i]->frame_number_offset<=framenumber) {
			if ((framenumber-framesets[i]->frame_number_offset)<framesets[i]->number_of_frames) return *framesets[i];
		}
	} 

	// if we miss, frame not found!
	cerr << "ERROR>> " << "Frame not found: " << framenumber << endl;
	throw;
}

size_t Frames::scope_framenumber(size_t framenumber) {
	Frameset& fs = find_frameset(framenumber);
	return framenumber-fs.frame_number_offset;
}

void Frames::load(size_t framenumber,Atoms& atoms,std::map<std::string,Atomselection>& atomselections) {
	// get frameset
	// call 
	if (framecache.find(framenumber)!=framecache.end()) {
		currentframe_i = framenumber;
	}
	else {
		if (framecache.size()>framecache_max) {
			framecache.clear();
		}
		// locate frameset 
		Frameset& fs = find_frameset(framenumber);
		// prepare empty frame
		Frame& cf = framecache[framenumber];
		
		// fill frame w/ data
		fs.read_frame(scope_framenumber(framenumber),cf);

		// postprocessing frame data
		// do any wrapping "before" creating coordinate sets, otherwise the coordinates will run out of sync
		if (wrapping) {
			cf.origin = cf.cofm(atoms,centergroup_selection);
			cf.wrap(); 					
		}	
		
		// for each group in atomselections, block data for performance
		for (std::map<std::string,Atomselection>::iterator asi=atomselections.begin();asi!=atomselections.end();asi++) {
			cf.push_selection(asi->second);
		}
	}
	currentframe_i = framenumber;
	
}


void Frames::load(size_t framenumber,Atoms& atoms,Atomselection& atomselection) {
	// get frameset
	// call 
	if (framecache.find(framenumber)!=framecache.end()) {
		currentframe_i = framenumber;
	}
	else {
		if (framecache.size()>framecache_max) {
			framecache.clear();
		}
		// locate frameset 
		Frameset& fs = find_frameset(framenumber);
		// prepare empty frame
		Frame& cf = framecache[framenumber];
		
		// fill frame w/ data
		fs.read_frame(scope_framenumber(framenumber),cf);

		// postprocessing frame data
		// do any wrapping "before" creating coordinate sets, otherwise the coordinates will run out of sync
		if (wrapping) {
			cf.origin = cf.cofm(atoms,centergroup_selection);
			cf.wrap(); 					
		}	
		
		// for each group in atomselections, block data for performance
		cf.push_selection(atomselection);
	}
	currentframe_i = framenumber;
	
}


Frame& Frames::current() {
	test_framenumber();
	return framecache[currentframe_i];
}



void Frames::clear_cache() {
	framecache.clear();
}



// constructor = analyze file and store frame locator information
void DCDFrameset::init(std::string fn,size_t fno) {
	
	// general information
	frame_number_offset = fno;
	filename = fn;
	
	// test if dcd file, if not throw!
	if (!detect(fn)) {
		cerr << "ERROR>> " << "file '" << fn << "' appears to not be a DCD file. INT Magic failed" << endl;
		throw;
	}
	
	// locator specific information
	
	int32_t marker; // fortran 4 byte marker		

	// if this is a dcd file just go ahead, otherwise assume a text file
	// where each line represents the path of a dcdfile.
	// this concept can be nested...

	ifstream dcdfile(Settings::get_filepath(fn).c_str());
	DCDHeader dcdheader;
	dcdfile.read((char*) &dcdheader,sizeof(dcdheader));
	dcdfile.seekg(23*sizeof(int32_t),ios_base::beg);

	dcdfile.read((char*) &marker,sizeof(int32_t));
	dcdfile.seekg(marker,ios_base::cur);
	dcdfile.read((char*) &marker,sizeof(int32_t));

	long noa;
	dcdfile.read((char*) &marker,sizeof(int32_t));
	dcdfile.read((char*) &noa,sizeof(long));
	dcdfile.read((char*) &marker,sizeof(int32_t));

//	if (number_of_atoms!=atoms.size()) { cerr << "Atom number mismatch (dcd)" << number_of_atoms << " vs. (pdb)" << atoms.size() << endl; throw; }
	// determine byte positions, i.e. read first frame

	number_of_frames = dcdheader.number_of_frames;
	number_of_atoms = noa;
	init_byte_pos = dcdfile.tellg();
	flag_ext_block1 = dcdheader.flag_ext_block1;
	flag_ext_block2 = dcdheader.flag_ext_block2;

	if (flag_ext_block1) {
		dcdfile.read((char*) &marker,sizeof(int32_t));
		block1_byte_offset = dcdfile.tellg()-init_byte_pos;
		dcdfile.seekg(marker,ios_base::cur);
		dcdfile.read((char*) &marker,sizeof(int32_t));
	} else {
		block1_byte_offset = dcdfile.tellg()-init_byte_pos;
	}

	dcdfile.read((char*) &marker,sizeof(int32_t));
	x_byte_offset = dcdfile.tellg()-init_byte_pos;
	dcdfile.seekg(number_of_atoms*sizeof(float),ios_base::cur);
	dcdfile.read((char*) &marker,sizeof(int32_t));

	dcdfile.read((char*) &marker,sizeof(int32_t));
	y_byte_offset = dcdfile.tellg()-init_byte_pos;
	dcdfile.seekg(number_of_atoms*sizeof(float),ios_base::cur);
	dcdfile.read((char*) &marker,sizeof(int32_t));

	dcdfile.read((char*) &marker,sizeof(int32_t));
	z_byte_offset = dcdfile.tellg()-init_byte_pos;
	dcdfile.seekg(number_of_atoms*sizeof(float),ios_base::cur);
	dcdfile.read((char*) &marker,sizeof(int32_t));

	if (flag_ext_block2) {
		dcdfile.read((char*) &marker,sizeof(int32_t));
		block2_byte_offset = dcdfile.tellg()-init_byte_pos;
		dcdfile.seekg(marker,ios_base::cur);
		dcdfile.read((char*) &marker,sizeof(int32_t));
	} else {
		block2_byte_offset = dcdfile.tellg()-init_byte_pos;
	}

	block_size_byte = dcdfile.tellg()-init_byte_pos;

	 // end dcd file read out
}

bool DCDFrameset::detect(const string filename) {
	ifstream dcdfile(Settings::get_filepath(filename).c_str(),ios::binary);
	char fp1[4];fp1[0]=0x54;fp1[1]=0x00;fp1[2]=0x00;fp1[3]=0x00;
	char fp2[4];fp2[0]=0x43;fp2[1]=0x4f;fp2[2]=0x52;fp2[3]=0x44;
	char buf[92]; dcdfile.read(buf,92);
	if ( (memcmp(&buf[0],&fp1[0],4)==0) && (memcmp(&buf[88],&fp1[0],4)==0) && (memcmp(&buf[4],&fp2[0],4)==0) ) return true;	else return false;
}

void DCDFrameset::read_frame(size_t internalframenumber,Frame& cf) {

	ifstream dcdfile(Settings::get_filepath(filename).c_str());
	
	double unit_cell_block[6];
	
	// This reads the unit cell
	if (flag_ext_block1) {
		dcdfile.seekg(init_byte_pos+internalframenumber*block_size_byte+block1_byte_offset,ios_base::beg);
		for (int i=0;i<6;i++) {
			dcdfile.read((char*) &(unit_cell_block[i]),sizeof(double)); 
		}
		// convert dcd first block to unit cell
		
		cf.unitcell.push_back(CartesianCoor3D(unit_cell_block[0],0,0));
		cf.unitcell.push_back(CartesianCoor3D(unit_cell_block[1],unit_cell_block[2],0));
		cf.unitcell.push_back(CartesianCoor3D(unit_cell_block[3],unit_cell_block[4],unit_cell_block[5]));
	}
	else {
		// signal that there is no unit cell!
		cerr << "ERROR>> " << " Frames w/o unit cell parameters are not supported" << endl;
		throw;
		// cf. = false;
	}
	
	cf.number_of_atoms = number_of_atoms;

	dcdfile.seekg(init_byte_pos+internalframenumber*block_size_byte+x_byte_offset,ios_base::beg);
	for (int i=0;i<number_of_atoms;i++) {
		float temp; dcdfile.read((char*) &temp,sizeof(float)); cf.x.push_back(temp);
	}
	dcdfile.seekg(init_byte_pos+internalframenumber*block_size_byte+y_byte_offset,ios_base::beg);
	for (int i=0;i<number_of_atoms;i++) {
		float temp; dcdfile.read((char*) &temp,sizeof(float)); cf.y.push_back(temp);
	}
	dcdfile.seekg(init_byte_pos+internalframenumber*block_size_byte+z_byte_offset,ios_base::beg);
	for (int i=0;i<number_of_atoms;i++) {
		float temp; dcdfile.read((char*) &temp,sizeof(float)); cf.z.push_back(temp);
	}

	// This block has unknown feature
	if (flag_ext_block2) {
		dcdfile.seekg(init_byte_pos+internalframenumber*block_size_byte+block2_byte_offset,ios_base::beg);
		for (int i=0;i<number_of_atoms;i++) {
			float temp; dcdfile.read((char*) &temp,sizeof(float));
			// skip it
			// cf.block2.push_back(temp);
		}
	}

}

// end of file
