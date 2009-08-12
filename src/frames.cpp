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
#include <boost/filesystem.hpp>
#include "xdrfile/xdrfile.h"
#include "xdrfile/xdrfile_xtc.h"
#include "xdrfile/xdrfile_trr.h"

// other headers
#include "atomselection.hpp"
#include "log.hpp"
#include "parameters.hpp"
#include "settings.hpp"


using namespace std;


void Frames::test_framenumber(size_t framenumber) {
	if ((framenumber>=size()) || (framenumber<0)) {
		cerr << "ERROR>> " << " framenumber (" << framenumber << ") integer out of bound" << endl;
	}
}

Frames::~Frames() {
	for (size_t i=0;i<framesets.size();i++) {
		delete framesets[i];
	}	
}

size_t Frames::add_frameset(const std::string filename,const std::string filetype,size_t first, size_t last, bool last_set, size_t stride, Atoms& atoms) {
	
	// Detect frame type
	if (filetype=="dcd") {
		// create an empty frameset and then initialize w/ filename
		Frameset* fsp = new DCDFrameset(filename, first,  last,  last_set,  stride, number_of_frames);
		framesets.push_back(fsp);
		// increase frame counter
		number_of_frames += fsp->number_of_frames;
		return fsp->number_of_frames;
	}
	else if (filetype=="dcdlist") {
		// read text file line by line
		ifstream listfile(filename.c_str());
		string line;
		size_t total_number_of_frames = 0;
		while (listfile >> line) {
			// skip empty lines or lines prepended w/ a '#' !
			if ((line.size() >0) && (line.substr(0,1) == "#")) continue;
			total_number_of_frames += add_frameset(line.c_str(), "dcd", first,  last,  last_set,  stride,  atoms);
		}
		return total_number_of_frames;
	}
	else if (filetype=="pdb") {
		Frameset* fsp = new PDBFrameset(filename, first,  last,  last_set,  stride,number_of_frames);
		framesets.push_back(fsp);
		// increase frame counter
		number_of_frames += fsp->number_of_frames;
		return fsp->number_of_frames;		
	}
	else if (filetype=="pdblist") {
		// read text file line by line
		ifstream listfile(filename.c_str());
		string line;
		size_t total_number_of_frames = 0;		
		while (listfile >> line) {
			// skip empty lines or lines prepended w/ a '#' !
			if ((line.size() >0) && (line.substr(0,1) == "#")) continue;			
			total_number_of_frames += add_frameset(line.c_str(), "pdb",first,  last,  last_set,  stride,  atoms);
		}
		return total_number_of_frames;		
	}
	else if (filetype=="xtc") {
		// create an empty frameset and then initialize w/ filename
		Frameset* fsp = new XTCFrameset(filename, first,  last,  last_set,  stride,number_of_frames);
		framesets.push_back(fsp);
		// increase frame counter
		number_of_frames += fsp->number_of_frames;
		return fsp->number_of_frames;
	}	
	else if (filetype=="trr") {
		// create an empty frameset and then initialize w/ filename
		Frameset* fsp = new TRRFrameset(filename, first,  last,  last_set,  stride,number_of_frames);
		framesets.push_back(fsp);
		// increase frame counter
		number_of_frames += fsp->number_of_frames;
		return fsp->number_of_frames;
	}		
	else {
		cerr << "ERROR>> filetype '" << filetype << "' not supported." << endl;
		throw;
	}
	
	return -1; 
}

size_t Frames::size() {
	return number_of_frames;
}

size_t Frames::cache_size() {
	return framecache.size();
}

Frameset& Frames::find_frameset(size_t framenumber) {
	test_framenumber(framenumber);
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
	return (framenumber - fs.frame_number_offset);
}

void Frames::load(size_t framenumber,Atoms& atoms,std::map<std::string,Atomselection>& atomselections) {
	// get frameset
	// call 	
	if (framecache.find(framenumber)!=framecache.end()) {
		currentframe_i = framenumber;
	}
	else {
		if ((Params::Inst()->limits.framecache_max>0) && framecache.size()>Params::Inst()->limits.framecache_max) {
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
		if ((Params::Inst()->limits.framecache_max>0) && framecache.size()>Params::Inst()->limits.framecache_max) {
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
	test_framenumber(currentframe_i);
	return framecache[currentframe_i];
}



void Frames::clear_cache() {
	framecache.clear();
}

size_t Frameset::init_frame_byte_offsets() {

	// test if a cache file exists
	// by default: cache-file name = filename.frame_byte_offsets
	using namespace boost::filesystem;
	
	path fpath(filename);
	string cache_filename = (path(Params::Inst()->runtime.config_rootpath) / fpath.filename()).string() + ".frame_byte_offsets"; // by default the cache files are in the config directory

	if (ifstream(cache_filename.c_str()).good()) {
		Info::Inst()->write(string("Reading Frame Offset Positions from Cache File: ")+cache_filename);
		read_frame_byte_offsets(cache_filename);
	} else {
		Info::Inst()->write("Scanning Trajectory for Frame Offset Position Information");
		scan_frame_byte_offsets();
		Info::Inst()->write(string("Writing Frame Offset Positions to Cache File: ")+cache_filename);		
		write_frame_byte_offsets(cache_filename);
	}

	// apply range modifier (first,last,stride)
	trim_frame_byte_offsets();
	
	return frame_byte_offsets.size();
}

void Frameset::read_frame_byte_offsets(string cache_filename) {
	ifstream cache_file(cache_filename.c_str());

	std::ios::streamoff frame_offset;	
	while (cache_file>>frame_offset) {
		frame_byte_offsets.push_back(frame_offset);
	}

}

void Frameset::write_frame_byte_offsets(string cache_filename) {
	ofstream cache_file(cache_filename.c_str());
	
	for(size_t i = 0; i < frame_byte_offsets.size(); ++i)
	{
		cache_file << frame_byte_offsets[i] << endl;
	}
}

void Frameset::trim_frame_byte_offsets() {
	
	vector<std::ios::streamoff> lfo;
	
	for(size_t i = 0; i < frame_byte_offsets.size(); ++i)
	{
		if (i<first) continue;
		if (last_set && (i>last))	break;
		if ((i % stride)==0) lfo.push_back(frame_byte_offsets[i]);		
	}
	
	if (frame_byte_offsets.size()!=lfo.size()) {
		Info::Inst()->write(string("Applied Range(first,last,stride) reduced number of frames from ")+to_s(frame_byte_offsets.size())+string(" to ")+to_s(lfo.size()));
	}
	
	frame_byte_offsets = lfo;
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

	ifstream dcdfile(fn.c_str());
	DCDHeader dcdheader;
	dcdfile.read((char*) &dcdheader,sizeof(dcdheader));
	dcdfile.seekg(23*sizeof(int32_t),ios_base::beg);

	dcdfile.read((char*) &marker,sizeof(int32_t));
	dcdfile.seekg(marker,ios_base::cur);
	dcdfile.read((char*) &marker,sizeof(int32_t));

	uint32_t noa;
	dcdfile.read((char*) &marker,sizeof(int32_t));
	dcdfile.read((char*) &noa,sizeof(uint32_t));
	dcdfile.read((char*) &marker,sizeof(int32_t));

//	if (number_of_atoms!=atoms.size()) { cerr << "Atom number mismatch (dcd)" << number_of_atoms << " vs. (pdb)" << atoms.size() << endl; throw; }
	// determine byte positions, i.e. read first frame

	number_of_frames = dcdheader.number_of_frames;
	number_of_atoms = static_cast<size_t>(noa);
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
	
	// this overwrites the value we get from the dcd header
	number_of_frames = init_frame_byte_offsets();
	
}

void DCDFrameset::scan_frame_byte_offsets() {
	// With DCD Files frame_byte_offsets are calculatable:
	std::ios::streamoff current_offset = init_byte_pos;
	for(size_t i = 0; i < number_of_frames; ++i)
	{
		frame_byte_offsets.push_back(i*block_size_byte + init_byte_pos);
	}
}

bool DCDFrameset::detect(const string filename) {
	ifstream dcdfile(filename.c_str(),ios::binary);
	char fp1[4];fp1[0]=0x54;fp1[1]=0x00;fp1[2]=0x00;fp1[3]=0x00;
	char fp2[4];fp2[0]=0x43;fp2[1]=0x4f;fp2[2]=0x52;fp2[3]=0x44;
	char buf[92]; dcdfile.read(buf,92);
	if ( (memcmp(&buf[0],&fp1[0],4)==0) && (memcmp(&buf[88],&fp1[0],4)==0) && (memcmp(&buf[4],&fp2[0],4)==0) ) return true;	else return false;
}

void DCDFrameset::read_frame(size_t internalframenumber,Frame& cf) {

	ifstream dcdfile(filename.c_str());
	
	double unit_cell_block[6];
	
	// This reads the unit cell
	if (flag_ext_block1) {
		dcdfile.seekg(frame_byte_offsets[internalframenumber]+block1_byte_offset,ios_base::beg);
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

	dcdfile.seekg(frame_byte_offsets[internalframenumber]+x_byte_offset,ios_base::beg);
	for (size_t i=0;i<number_of_atoms;i++) {
		float temp; dcdfile.read((char*) &temp,sizeof(float)); cf.x.push_back(temp);
	}
	dcdfile.seekg(frame_byte_offsets[internalframenumber]+y_byte_offset,ios_base::beg);
	for (size_t i=0;i<number_of_atoms;i++) {
		float temp; dcdfile.read((char*) &temp,sizeof(float)); cf.y.push_back(temp);
	}
	dcdfile.seekg(frame_byte_offsets[internalframenumber]+z_byte_offset,ios_base::beg);
	for (size_t i=0;i<number_of_atoms;i++) {
		float temp; dcdfile.read((char*) &temp,sizeof(float)); cf.z.push_back(temp);
	}

	// This block has unknown feature
	if (flag_ext_block2) {
		dcdfile.seekg(frame_byte_offsets[internalframenumber]+block2_byte_offset,ios_base::beg);
		for (size_t i=0;i<number_of_atoms;i++) {
			float temp; dcdfile.read((char*) &temp,sizeof(float));
			// skip it
			// cf.block2.push_back(temp);
		}
	}
}

// PDB Frameset implementation:


// constructor = analyze file and store frame locator information
void PDBFrameset::init(std::string fn,size_t fno) {
	
	// general information
	frame_number_offset = fno;
	filename = fn;
	
	// test if dcd file, if not throw!
	if (!detect(fn)) {
		cerr << "ERROR>> " << "file '" << fn << "' appears not to be a PDB file" << endl;
		throw;
	}
	
	ifstream pdbfile(filename.c_str());	
	
	// scan the first frame
	// also load the first CRYST1 as default unit cell
	string line;
	default_uc.resize(3);
	size_t number_of_atom_entries=0;
	while (getline(pdbfile,line)) {
		// if entry == crystal -> unit cell
		if ( line.substr(0,6) != "CRYST1" ) continue;

		double a = atof(line.substr(6,9).c_str());
		double b = atof(line.substr(15,9).c_str());
		double c = atof(line.substr(24,9).c_str());
		double alpha = atof(line.substr(33,7).c_str()) * 2 * M_PI / 360 ;
		double beta =  atof(line.substr(40,7).c_str()) * 2 * M_PI / 360 ;
		double gamma = atof(line.substr(47,7).c_str()) * 2 * M_PI / 360 ;															
		default_uc[0]=CartesianCoor3D(a,0,0);
		default_uc[1]=CartesianCoor3D(b*cos(gamma),b*sin(gamma),0);
		default_uc[2]=CartesianCoor3D(c*cos(beta),c*cos(alpha)*sin(beta),c*sin(alpha)*sin(beta));
		break;
	}
	
	pdbfile.seekg(ios::beg);
	
	while (getline(pdbfile,line)) {					
		if ( line.substr(0,6) == "ATOM  " ) number_of_atom_entries++;
		if ( line.substr(0,3) == "END" ) break;							
	}	
	// locator specific information
	
	pdbfile.close();
	
	number_of_frames = init_frame_byte_offsets();
	
}

void PDBFrameset::scan_frame_byte_offsets() {
	
	ifstream pdbfile(filename.c_str());	
	
	size_t cfn=0;
	size_t lines=0;
	bool lastentryisater = false;
	std::ios::streamoff filepos = pdbfile.tellg();
	while (true) {
		string line;
		if (!getline(pdbfile,line)) break;

		if ( line.substr(0,3) == "END" ) {
			frame_byte_offsets.push_back(filepos);
			filepos = pdbfile.tellg();
			lastentryisater = true;
		}
		else {
			lastentryisater = false;
		}
	}
	
	// save last unterminated frame
	if (!lastentryisater) {
			frame_byte_offsets.push_back(filepos);
	}
	
}

bool PDBFrameset::detect(const string filename) {
	ifstream pdbfile(filename.c_str());
	string line; 
	bool atomentry=false;
	while (getline(pdbfile,line)) {
		if ( line.substr(0,6) == "ATOM  " ) { 
			atomentry = true;
			break;
		}
	}
	if (atomentry) return true; else return false;
}

void PDBFrameset::read_frame(size_t internalframenumber,Frame& cf) {

	ifstream pdbfile(filename.c_str());
	
	vector<CartesianCoor3D> uc = default_uc; 

	pdbfile.seekg(frame_byte_offsets[internalframenumber]);
	
	size_t number_of_atom_entries = 0;
	string line;
	while (getline(pdbfile,line)) {
		// if entry == crystal -> unit cell
		if ( line.substr(0,6) == "CRYST1" ) {
			double a = atof(line.substr(6,9).c_str());
			double b = atof(line.substr(15,9).c_str());
			double c = atof(line.substr(24,9).c_str());
			double alpha = atof(line.substr(33,7).c_str()) * 2 * M_PI / 360 ;
			double beta =  atof(line.substr(40,7).c_str()) * 2 * M_PI / 360 ;
			double gamma = atof(line.substr(47,7).c_str()) * 2 * M_PI / 360 ;															
			uc[0]=CartesianCoor3D(a,0,0);
			uc[1]=CartesianCoor3D(b*cos(gamma),b*sin(gamma),0);
			uc[2]=CartesianCoor3D(c*cos(beta),c*cos(alpha)*sin(beta),c*sin(alpha)*sin(beta));
		}
						
		// if entry == atom -> coordinates
		if ( line.substr(0,6) == "ATOM  " ) {
			number_of_atom_entries++;			
			cf.x.push_back(atof(line.substr(30,8).c_str()) );
			cf.y.push_back(atof(line.substr(38,8).c_str()) );
			cf.z.push_back(atof(line.substr(46,8).c_str()) );
		}
		
		// if entry == END -> break
		if ( line.substr(0,3) == "END" ) break;
	}

	cf.unitcell = uc;
	cf.number_of_atoms = number_of_atom_entries;
}


// XTC Frameset implementation:

// constructor = analyze file and store frame locator information
void XTCFrameset::init(std::string fn,size_t fno) {
	
	// general information
	frame_number_offset = fno;
	filename = fn;
	
	// test if xtc file, if not throw!
	if (!detect(fn)) {
		cerr << "ERROR>> " << "file '" << fn << "' appears not to be a XTC file" << endl;
		throw;
	}
	
  	/* This function returns the number of atoms in the xtc file in *natoms */
	int natoms = 0;
  	int retval = read_xtc_natoms(const_cast<char*>(filename.c_str()),&natoms);
	number_of_atoms = natoms;

	number_of_frames = init_frame_byte_offsets();
	
	return;
}

void XTCFrameset::scan_frame_byte_offsets() {
	
    XDRFILE* p_xdrfile =  xdrfile_open(const_cast<char*>(filename.c_str()),"r");	
	FILE* fp = get_filepointer(p_xdrfile);
	int step =0; float t =0;
	float prec =  1000.0;
	rvec* coords = (rvec*) malloc(sizeof(rvec)*number_of_atoms);
	matrix box; 	

	size_t cfn=0;
	while (true) {
		std::ios::streamoff filepos = ftell( fp );
		int retval = read_xtc(p_xdrfile,number_of_atoms,&step,&t,box,coords,&prec);
		if (retval != exdrOK) break;
		frame_byte_offsets.push_back(filepos);
		cfn++;		
	}
	
	free(coords);
	xdrfile_close(p_xdrfile);
}

bool XTCFrameset::detect(const string filename) {

  	/* This function returns the number of atoms in the xtc file in *natoms */
	int natoms = 0;
  	int retval = read_xtc_natoms(const_cast<char*>(filename.c_str()),&natoms);
	
	if (retval==exdrOK) return true; else return false;
}

void XTCFrameset::read_frame(size_t internalframenumber,Frame& cf) {
	
	// read a specific frame
    XDRFILE* p_xdrfile =  xdrfile_open(const_cast<char*>(filename.c_str()),"r");	
	FILE* fp = get_filepointer(p_xdrfile);
	int step =0; float t =0;
	float prec =  1000.0;
	rvec* coords = (rvec*) malloc(sizeof(rvec)*number_of_atoms);
	matrix box; 	

	fseek( fp, frame_byte_offsets[internalframenumber] , SEEK_SET );
	int retval = read_xtc(p_xdrfile,number_of_atoms,&step,&t,box,coords,&prec);
		
	vector<CartesianCoor3D> uc(3); 

	uc[0] = 10.0* CartesianCoor3D(box[0][0],box[0][1],box[0][2]);
	uc[1] = 10.0* CartesianCoor3D(box[1][0],box[1][1],box[1][2]);
	uc[2] = 10.0* CartesianCoor3D(box[2][0],box[2][1],box[2][2]);

	for(size_t i = 0; i < number_of_atoms; ++i)
	{
		cf.x.push_back(10.0* coords[i][0]);
		cf.y.push_back(10.0* coords[i][1]);
		cf.z.push_back(10.0* coords[i][2]);				
	}
	

	cf.unitcell = uc;
	cf.number_of_atoms = number_of_atoms;
	
	free(coords);
	
	xdrfile_close(p_xdrfile);
	return;
}



// TRR Frameset implementation:

// constructor = analyze file and store frame locator information
void TRRFrameset::init(std::string fn,size_t fno) {
	
	// general information
	frame_number_offset = fno;
	filename = fn;
	
	// test if xtc file, if not throw!
	if (!detect(fn)) {
		cerr << "ERROR>> " << "file '" << fn << "' appears not to be a TRR file" << endl;
		throw;
	}
	
  	/* This function returns the number of atoms in the xtc file in *natoms */
	int natoms = 0;
  	int retval = read_trr_natoms(const_cast<char*>(filename.c_str()),&natoms);
	number_of_atoms = natoms;

	number_of_frames = init_frame_byte_offsets();
	
	return;
}

void TRRFrameset::scan_frame_byte_offsets() {
	
    XDRFILE* p_xdrfile =  xdrfile_open(const_cast<char*>(filename.c_str()),"r");	
	FILE* fp = get_filepointer(p_xdrfile);
	int step =0; float t =0;
	float lambda =  0.0;
	rvec* coords = (rvec*) malloc(sizeof(rvec)*number_of_atoms);
	matrix box; 	

	size_t cfn=0;
	while (true) {
		std::ios::streamoff filepos = ftell( fp );
		int retval = read_trr(p_xdrfile,number_of_atoms,&step,&t,&lambda,box,coords,NULL,NULL);
		if (retval != exdrOK) break;
		frame_byte_offsets.push_back(filepos);
		cfn++;		
	}
	
	free(coords);
	xdrfile_close(p_xdrfile);
}

bool TRRFrameset::detect(const string filename) {

  	/* This function returns the number of atoms in the xtc file in *natoms */
	int natoms = 0;
  	int retval = read_trr_natoms(const_cast<char*>(filename.c_str()),&natoms);
	
	if (retval==exdrOK) return true; else return false;
}

void TRRFrameset::read_frame(size_t internalframenumber,Frame& cf) {
	
	// read a specific frame
    XDRFILE* p_xdrfile =  xdrfile_open(const_cast<char*>(filename.c_str()),"r");	
	FILE* fp = get_filepointer(p_xdrfile);
	int step =0; float t =0;
	float lambda = 0;
	rvec* coords = (rvec*) malloc(sizeof(rvec)*number_of_atoms);
	matrix box; 	

	fseek( fp, frame_byte_offsets[internalframenumber] , SEEK_SET );
	int retval = read_trr(p_xdrfile,number_of_atoms,&step,&t,&lambda,box,coords,NULL,NULL);
		
	vector<CartesianCoor3D> uc(3); 

	uc[0] = 10.0* CartesianCoor3D(box[0][0],box[0][1],box[0][2]);
	uc[1] = 10.0* CartesianCoor3D(box[1][0],box[1][1],box[1][2]);
	uc[2] = 10.0* CartesianCoor3D(box[2][0],box[2][1],box[2][2]);

	for(size_t i = 0; i < number_of_atoms; ++i)
	{
		cf.x.push_back(10.0* coords[i][0]);
		cf.y.push_back(10.0* coords[i][1]);
		cf.z.push_back(10.0* coords[i][2]);				
	}
	

	cf.unitcell = uc;
	cf.number_of_atoms = number_of_atoms;
	
	free(coords);
	
	xdrfile_close(p_xdrfile);
	return;
}

// end of file
