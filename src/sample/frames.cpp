/** \file
This file contains a management class for framesets and specifies framesets for different trajectory types.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/

// direct header
#include "sample/frames.hpp"

// standard header
#include <fstream>

// special library headers
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include "../vendor/xdrfile-1.1.1/include/xdrfile.h"
#include "../vendor/xdrfile-1.1.1/include/xdrfile_xtc.h"
#include "../vendor/xdrfile-1.1.1/include/xdrfile_trr.h"

// other headers
#include "sample/atomselection.hpp"
#include "sample/center_of_mass.hpp"
#include "control.hpp"
#include "log.hpp"



using namespace std;


void Frames::test_framenumber(size_t framenumber) {
	if ((framenumber>=size()) || (framenumber<0)) {
		cerr << "ERROR>> " << " framenumber (" << framenumber << ") integer out of bound" << endl;
	}
}

size_t Frames::add_frameset(const std::string filename,const std::string filetype,size_t first, size_t last, bool last_set, size_t stride,const std::string index_filename,bool index_default,size_t clones) {
	
    if (clones==0) return 0;
	
	// Detect frame type
	if (filetype=="dcd") {
        size_t this_number_of_frames =0;
		// create an empty frameset and then initialize w/ filename
		DCDFrameset* fsp = new DCDFrameset(filename, number_of_frames);
		if (!index_default) {
			fsp->load_index(index_filename);		
		} else {
	        fsp->generate_index();			
		}
        fsp->trim_index(first,  last,  last_set,  stride);
		framesets.push_back(fsp);
		
        Frameset* p_frameset = static_cast<Frameset*>(fsp);
						
		// increase frame counter
		number_of_frames += fsp->number_of_frames;
		this_number_of_frames += fsp->number_of_frames;
        for(size_t i = 1; i < clones; ++i)
        {
    		Frameset* fspc = new CloneFrameset(p_frameset,number_of_frames);
			framesets.push_back(fspc);
    		number_of_frames += fspc->number_of_frames;	
    		this_number_of_frames += fspc->number_of_frames;		    
		}
		
		return this_number_of_frames;
	}
	else if (filetype=="dcdlist") {
		// read text file line by line
		ifstream listfile(filename.c_str());
		string line;
		size_t total_number_of_frames = 0;
		while (listfile >> line) {
			// skip empty lines or lines prepended w/ a '#' !
			if ((line.size() >0) && (line.substr(0,1) == "#")) continue;
            for(size_t i = 0; i < clones; ++i)
            {
    			total_number_of_frames += add_frameset(line.c_str(), "dcd", first,  last,  last_set,  stride,index_filename,index_default,1);                
            }
		}
		return total_number_of_frames;
	}
	else if (filetype=="pdb") {
        size_t this_number_of_frames =0;

		FileFrameset* fsp = new PDBFrameset(filename,number_of_frames);
		if (!boost::filesystem::exists(index_filename)) {
            Err::Inst()->write(string("Frameset ")+filename+string(" requires an Index file. "));
            Err::Inst()->write(string("No Frameset Index found at given location: ")+index_filename);
            throw;
		}
	    fsp->load_index(index_filename);
        fsp->trim_index(first,  last,  last_set,  stride);
		framesets.push_back(fsp);
		
		Frameset* p_frameset = static_cast<Frameset*>(fsp);
						
		// increase frame counter
		number_of_frames += fsp->number_of_frames;
		this_number_of_frames += fsp->number_of_frames;
        for(size_t i = 1; i < clones; ++i)
        {
    		Frameset* fspc = new CloneFrameset(p_frameset,number_of_frames);
			framesets.push_back(fspc);
    		number_of_frames += fspc->number_of_frames;	
    		this_number_of_frames += fspc->number_of_frames;		    
		}
		
		return this_number_of_frames;		
	}
	else if (filetype=="pdblist") {
		// read text file line by line
		ifstream listfile(filename.c_str());
		string line;
		size_t total_number_of_frames = 0;		
		while (listfile >> line) {
			// skip empty lines or lines prepended w/ a '#' !
			if ((line.size() >0) && (line.substr(0,1) == "#")) continue;
			for(size_t i = 0; i < clones; ++i)
            {			
			    total_number_of_frames += add_frameset(line.c_str(), "pdb",first,  last,  last_set,  stride,index_filename,index_default,1);
			}
		}
		return total_number_of_frames;		
	}
	else if (filetype=="xtc") {
        size_t this_number_of_frames =0;
        	    
		// create an empty frameset and then initialize w/ filename
		FileFrameset* fsp = new XTCFrameset(filename, number_of_frames);
		if (!boost::filesystem::exists(index_filename)) {
            Err::Inst()->write(string("Frameset ")+filename+string(" requires an Index file. "));
            Err::Inst()->write(string("No Frameset Index found at given location: ")+index_filename);
            throw;
		}
		
		fsp->load_index(index_filename);
		fsp->trim_index(first,  last,  last_set,  stride);
		framesets.push_back(fsp);
		
		Frameset* p_frameset = static_cast<Frameset*>(fsp);
						
		// increase frame counter
		number_of_frames += fsp->number_of_frames;
		this_number_of_frames += fsp->number_of_frames;
        for(size_t i = 1; i < clones; ++i)
        {
    		Frameset* fspc = new CloneFrameset(p_frameset,number_of_frames);
			framesets.push_back(fspc);
    		number_of_frames += fspc->number_of_frames;	
    		this_number_of_frames += fspc->number_of_frames;		    
		}
		
		return this_number_of_frames;
	}	
	else if (filetype=="trr") {
        size_t this_number_of_frames =0;
        	    
		// create an empty frameset and then initialize w/ filename
		FileFrameset* fsp = new TRRFrameset(filename, number_of_frames);
		if (!boost::filesystem::exists(index_filename)) {
            Err::Inst()->write(string("Frameset ")+filename+string(" requires an Index file. "));
            Err::Inst()->write(string("No Frameset Index found at given location: ")+index_filename);
            throw;
		}
		
		fsp->load_index(index_filename);  
		fsp->trim_index(first,  last,  last_set,  stride);
		framesets.push_back(fsp);

		Frameset* p_frameset = static_cast<Frameset*>(fsp);
						
		// increase frame counter
		number_of_frames += fsp->number_of_frames;
		this_number_of_frames += fsp->number_of_frames;
        for(size_t i = 1; i < clones; ++i)
        {
    		Frameset* fspc = new CloneFrameset(p_frameset,number_of_frames);
			framesets.push_back(fspc);
    		number_of_frames += fspc->number_of_frames;	
    		this_number_of_frames += fspc->number_of_frames;		    
		}
		
		return this_number_of_frames;
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

Frame Frames::load(size_t framenumber) {
	Frameset& fs = find_frameset(framenumber);

	// prepare empty frame
    Frame cf;		
	fs.read_frame(framenumber,cf);
	
    return cf;
}

void FileFrameset::trim_index(size_t first,size_t last,bool last_set,size_t stride) {
	
	vector<std::ios::streamoff> lfo;
	
	for(size_t i = 0; i < frameset_index_.size(); ++i)
	{
		if (i<first) continue;
		if (last_set && (i>last))	break;
		if ((i % stride)==0) lfo.push_back(frameset_index_[i]);		
	}
	
	if (frameset_index_.size()!=lfo.size()) {
		Info::Inst()->write(string("Applied Range(first,last,stride) reduced number of frames from ")+boost::lexical_cast<string>(frameset_index_.size())+string(" to ")+boost::lexical_cast<string>(lfo.size()));
	}

    frameset_index_.clear();
    for(size_t i = 0; i < lfo.size(); ++i)
    {
        frameset_index_.push_back(lfo[i]);
    }
    number_of_frames = frameset_index_.size();
}
void FileFrameset::load_index(std::string index_filename) {
    
    frameset_index_.load(index_filename);
    std::vector<char> sig = FramesetIndex::generate_signature(filename);

    if (sig!=frameset_index_.get_signature()) {
        Err::Inst()->write(string("FrameIndex file (")+index_filename+string(") incompatible with selected trajectory: ")+filename);
        throw;
    }    
}

void FileFrameset::save_index(std::string index_filename) {
    frameset_index_.save(index_filename);
}

void DCDFrameset::generate_index() {
    frameset_index_.clear();
    frameset_index_.set_signature(FramesetIndex::generate_signature(filename));
    
    for(size_t i = 0; i < number_of_frames; ++i)
	{
        frameset_index_.push_back(i*block_size_byte + init_byte_pos);
	}	
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
	if (!dcdfile.is_open()) {
		dcdfile.open(filename.c_str(),ios::binary);
	} else {
		dcdfile.seekg(0,ios_base::beg);
	}
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
//	number_of_frames = init_frame_byte_offsets();
}

bool DCDFrameset::detect(const string filename) {
	if (!dcdfile.is_open()) {
		dcdfile.open(filename.c_str(),ios::binary);
	} else {
		dcdfile.seekg(0,ios_base::beg);
	}
	
	char fp1[4];fp1[0]=0x54;fp1[1]=0x00;fp1[2]=0x00;fp1[3]=0x00;
	char fp2[4];fp2[0]=0x43;fp2[1]=0x4f;fp2[2]=0x52;fp2[3]=0x44;
	char buf[92]; dcdfile.read(buf,92);
	if ( (memcmp(&buf[0],&fp1[0],4)==0) && (memcmp(&buf[88],&fp1[0],4)==0) && (memcmp(&buf[4],&fp2[0],4)==0) ) return true;	else return false;
}

void DCDFrameset::close() {
	dcdfile.close();
}

DCDFrameset::~DCDFrameset() {
	try {
		close();
	} catch (...) {
		Warn::Inst()->write("File close operation failed!");
	}
}

void DCDFrameset::read_frame(size_t framenumber,Frame& cf) {

    size_t internalframenumber = framenumber - frame_number_offset;

	if (!dcdfile.is_open()) {
		dcdfile.open(filename.c_str(),ios::binary);
	}
	
	double unit_cell_block[6];
	
	// This reads the unit cell
	if (flag_ext_block1) {
		dcdfile.seekg(frameset_index_[internalframenumber]+block1_byte_offset,ios_base::beg);
		for (int i=0;i<6;i++) {
			dcdfile.read((char*) &(unit_cell_block[i]),sizeof(double)); 
		}
		// convert dcd first block to unit cell
		
		cf.unitcell.push_back(CartesianCoor3D(unit_cell_block[0],0,0));
		cf.unitcell.push_back(CartesianCoor3D(unit_cell_block[1],unit_cell_block[2],0));
		cf.unitcell.push_back(CartesianCoor3D(unit_cell_block[3],unit_cell_block[4],unit_cell_block[5]));
	}
	else {
		cf.unitcell.push_back(CartesianCoor3D(0,0,0));
		cf.unitcell.push_back(CartesianCoor3D(0,0,0));
		cf.unitcell.push_back(CartesianCoor3D(0,0,0));		
	}
	
	cf.number_of_atoms = number_of_atoms;

	dcdfile.seekg(frameset_index_[internalframenumber]+x_byte_offset,ios_base::beg);
	for (size_t i=0;i<number_of_atoms;i++) {
		float temp; dcdfile.read((char*) &temp,sizeof(float)); cf.x.push_back(temp);
	}
	dcdfile.seekg(frameset_index_[internalframenumber]+y_byte_offset,ios_base::beg);
	for (size_t i=0;i<number_of_atoms;i++) {
		float temp; dcdfile.read((char*) &temp,sizeof(float)); cf.y.push_back(temp);
	}
	dcdfile.seekg(frameset_index_[internalframenumber]+z_byte_offset,ios_base::beg);
	for (size_t i=0;i<number_of_atoms;i++) {
		float temp; dcdfile.read((char*) &temp,sizeof(float)); cf.z.push_back(temp);
	}

	// This block has unknown feature
	if (flag_ext_block2) {
		dcdfile.seekg(frameset_index_[internalframenumber]+block2_byte_offset,ios_base::beg);
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
	bool unitcell=false;
	while (getline(pdbfile,line)) {
		// if entry == crystal -> unit cell
		if ( line.substr(0,6) != "CRYST1" ) continue;
		unitcell=true;
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
	
	if (!unitcell) {
		default_uc[0]=CartesianCoor3D(0,0,0);
		default_uc[1]=CartesianCoor3D(0,0,0);
		default_uc[2]=CartesianCoor3D(0,0,0);		
	}
	
	pdbfile.seekg(ios::beg);
	
	while (getline(pdbfile,line)) {					
		if ( line.substr(0,6) == "ATOM  " ) number_of_atom_entries++;
		if ( line.substr(0,3) == "END" ) break;							
	}	
	// locator specific information
	
	pdbfile.close();
	
}


void PDBFrameset::generate_index() {
	
    std::vector<std::streamoff> frame_byte_offsets;
	ifstream pdbfile(filename.c_str());	
	
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

	frameset_index_.clear();
    frameset_index_.set_signature(FramesetIndex::generate_signature(filename));
		
    for(size_t i = 0; i < frame_byte_offsets.size(); ++i)
	{
        frameset_index_.push_back(frame_byte_offsets[i]);

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

void PDBFrameset::read_frame(size_t framenumber,Frame& cf) {

    size_t internalframenumber = framenumber - frame_number_offset;

	ifstream pdbfile(filename.c_str());
	
	vector<CartesianCoor3D> uc = default_uc; 

	pdbfile.seekg(frameset_index_[internalframenumber]);
	
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
  	if (read_xtc_natoms(const_cast<char*>(filename.c_str()),&natoms) != exdrOK) {
		cerr << "ERROR>> " << " initializing file: " << fn << endl;
        throw;
  	}
	number_of_atoms = natoms;
	
	return;
}

void XTCFrameset::generate_index() {
	
    std::vector<size_t> frame_byte_offsets;
	if (p_xdrfile==NULL) {
	    p_xdrfile =  xdrfile_open(const_cast<char*>(filename.c_str()),"r");	
	    if (p_xdrfile==NULL) {
	        Err::Inst()->write(string("Unable to open file: ")+filename);
	        throw;
	    }		
	}
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
		if (cfn<10) {
            Info::Inst()->write(string("Reading frame ")+boost::lexical_cast<string>(cfn));
		} else if (cfn<100) {
		    if ((cfn%10)==0) {
                Info::Inst()->write(string("Reading frame ")+boost::lexical_cast<string>(cfn));		        
		    }
		} else if ((cfn%100)==0) {
            Info::Inst()->write(string("Reading frame ")+boost::lexical_cast<string>(cfn));		        		    
		}
		
		cfn++;		
	}
	
	free(coords);
	
	frameset_index_.clear();
    frameset_index_.set_signature(FramesetIndex::generate_signature(filename));
    		
    for(size_t i = 0; i < frame_byte_offsets.size(); ++i)
	{
        frameset_index_.push_back(frame_byte_offsets[i]);
	}
}

bool XTCFrameset::detect(const string filename) {

  	/* This function returns the number of atoms in the xtc file in *natoms */
	int natoms = 0;
  	int retval = read_xtc_natoms(const_cast<char*>(filename.c_str()),&natoms);
	
	if (retval==exdrOK) return true; else return false;
}

void XTCFrameset::close() {
	xdrfile_close(p_xdrfile);
	p_xdrfile=NULL;
}

XTCFrameset::~XTCFrameset() {
	try {
		close();
	} catch (...) {
		Warn::Inst()->write("File close operation failed!");
	}
}

void XTCFrameset::read_frame(size_t framenumber,Frame& cf) {
	
	size_t internalframenumber = framenumber - frame_number_offset;
    
	// read a specific frame
	if (p_xdrfile==NULL) {
	    p_xdrfile =  xdrfile_open(const_cast<char*>(filename.c_str()),"r");	
	    if (p_xdrfile==NULL) {
	        Err::Inst()->write(string("Unable to open file: ")+filename);
	        throw;
	    }		
	}
	FILE* fp = get_filepointer(p_xdrfile);
	int step =0; float t =0;
	float prec =  1000.0;
	rvec* coords = (rvec*) malloc(sizeof(rvec)*number_of_atoms);
	matrix box; 	

	fseek( fp, frameset_index_[internalframenumber] , SEEK_SET );
	read_xtc(p_xdrfile,number_of_atoms,&step,&t,box,coords,&prec);
		
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
  	read_trr_natoms(const_cast<char*>(filename.c_str()),&natoms);
	number_of_atoms = natoms;

	return;
}

void TRRFrameset::generate_index() {
	
    std::vector<std::streamoff> frame_byte_offsets;
	if (p_xdrfile==NULL) {
	    p_xdrfile =  xdrfile_open(const_cast<char*>(filename.c_str()),"r");	
	    if (p_xdrfile==NULL) {
	        Err::Inst()->write(string("Unable to open file: ")+filename);
	        throw;
	    }		
	}
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
		if (cfn<10) {
            Info::Inst()->write(string("Reading frame ")+boost::lexical_cast<string>(cfn));
		} else if (cfn<100) {
		    if ((cfn%10)==0) {
                Info::Inst()->write(string("Reading frame ")+boost::lexical_cast<string>(cfn));		        
		    }
		} else if ((cfn%100)==0) {
            Info::Inst()->write(string("Reading frame ")+boost::lexical_cast<string>(cfn));		        		    
		}
		
		cfn++;		
	}
	
	free(coords);
	
	frameset_index_.clear();
    frameset_index_.set_signature(FramesetIndex::generate_signature(filename));
    		
    for(size_t i = 0; i < frame_byte_offsets.size(); ++i)
	{
        frameset_index_.push_back(frame_byte_offsets[i]);
	}    
}

bool TRRFrameset::detect(const string filename) {

  	/* This function returns the number of atoms in the xtc file in *natoms */
	int natoms = 0;
  	int retval = read_trr_natoms(const_cast<char*>(filename.c_str()),&natoms);
	
	if (retval==exdrOK) return true; else return false;
}

void TRRFrameset::close() {
	xdrfile_close(p_xdrfile);
	p_xdrfile=NULL;
}

TRRFrameset::~TRRFrameset() {
	try {
		close();
	} catch (...) {
		Warn::Inst()->write("File close operation failed!");
	}
}

void TRRFrameset::read_frame(size_t framenumber,Frame& cf) {
	
	size_t internalframenumber = framenumber - frame_number_offset;
    
	// read a specific frame
	if (p_xdrfile==NULL) {
	    p_xdrfile =  xdrfile_open(const_cast<char*>(filename.c_str()),"r");	
	    if (p_xdrfile==NULL) {
	        Err::Inst()->write(string("Unable to open file: ")+filename);
	        throw;
	    }		
	}
	FILE* fp = get_filepointer(p_xdrfile);
	int step =0; float t =0;
	float lambda = 0;
	rvec* coords = (rvec*) malloc(sizeof(rvec)*number_of_atoms);
	matrix box; 	

	fseek( fp, frameset_index_[internalframenumber] , SEEK_SET );
	read_trr(p_xdrfile,number_of_atoms,&step,&t,&lambda,box,coords,NULL,NULL);
		
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
	
}

CloneFrameset::CloneFrameset(Frameset* original, size_t nof) {
    p_originalframeset = original;
    frame_number_offset = nof;
    number_of_frames = original->number_of_frames;
    number_of_atoms = original->number_of_atoms;
}

void CloneFrameset::read_frame(size_t framenumber,Frame& cf) {
    // delegate to original frameset
    size_t internalframenumber = framenumber - frame_number_offset;
	
    p_originalframeset->read_frame(internalframenumber+p_originalframeset->frame_number_offset,cf);
}

// end of file
