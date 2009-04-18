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

// special library headers

// other headers
#include "atomselection.hpp"
#include "settings.hpp"

using namespace std;


//bool detect_dcdfile(const string dcdfilename) {
//	ifstream dcdfile(Settings::get_filepath(dcdfilename).c_str(),ios::binary);
//	char fp1[4];fp1[0]=0x54;fp1[1]=0x00;fp1[2]=0x00;fp1[3]=0x00;
//	char fp2[4];fp2[0]=0x43;fp2[1]=0x4f;fp2[2]=0x52;fp2[3]=0x44;
//	char buf[92]; dcdfile.read(buf,92);
//	if ( (memcmp(&buf[0],&fp1[0],4)==0) && (memcmp(&buf[88],&fp1[0],4)==0) && (memcmp(&buf[4],&fp2[0],4)==0) ) return true;	else return false;
//}

void Frames::test_framenumber() {
	if ((currentframe_i>=totalframes) || (currentframe_i<0)) {
		cerr << "ERROR>> " << " currentframe integer out of bound" << endl;
	}
}

FrameFilePosLocator& Frames::find_locator(int framenumber) {
	test_framenumber();
	for (std::vector<FrameFilePosLocator>::iterator ffpli=framefileposlocators.begin();ffpli!=framefileposlocators.end();ffpli++) {
		if ( (framenumber>=ffpli->frame_number_offset) && ((framenumber-ffpli->frame_number_offset)<ffpli->number_of_frames)) {
			return ffpli;
		}
	}
	// else: locator not found!
	throw;
}

Frame& Frames::current() {
	test_framenumber();
	return framecache[currentframe_i];
}

void Frames::load(int framenumber,std::vector<Atomselection> atomselections) { 
	if (framecache.find(framenumber)!=framecache.end()) {
		currentframe_i = framenumber;
	}
	else {
		if (framecache.size()>=framecache_max) framecache.clear();
		Frame& cf = framecache[framenumber]; // create an empty frame
		// load the real data
		read_data(framenumber,cf); // member function overloaded by specialized class
		// do postprocessing
		if (wrapping) {
			cf.origin = cf.cofm(atoms,atomselections[centergroup]);
			cf.wrap(); 					
		}	
		// for each group in atomselections, block data for performance
		for (std::vector<Atomselection>::iterator asi=atomselections.begin();asi!=atomselections.end();asi++) {
			cf.push_selection(asi);
		}
	}
	currentframe_i = framenumber;
}

void Frames::clear_cache() {
	framecache.clear();
}

// specialization: 
// each subclass has to implement: 
// addfile()
// readfile()



void DCDFrames::add_file(const std::string dcdfilename,Atoms& atoms,int recursion_trigger) {

	int32_t marker; // fortran 4 byte marker		

	if (recursion_trigger<1) { cerr << "DcdFiles::add_file. Too much recursion. Escaping." << endl; throw; }

	// if this is a dcd file just go ahead, otherwise assume a text file
	// where each line represents the path of a dcdfile.
	// this concept can be nested...

	ifstream dcdfile(Settings::get_filepath(dcdfilename).c_str());
	DcdHeader dcdheader;
	dcdfile.read((char*) &dcdheader,sizeof(dcdheader));
#ifdef SASSIM_DEBUG
	print_dcdheader(dcdheader);
#endif
	dcdfile.seekg(23*sizeof(int32_t),ios_base::beg);

	dcdfile.read((char*) &marker,sizeof(int32_t));
	dcdfile.seekg(marker,ios_base::cur);
	dcdfile.read((char*) &marker,sizeof(int32_t));

	uint32_t number_of_atoms;
	dcdfile.read((char*) &marker,sizeof(int32_t));
	dcdfile.read((char*) &number_of_atoms,sizeof(uint32_t));
	dcdfile.read((char*) &marker,sizeof(int32_t));

	if (number_of_atoms!=atoms.size()) { cerr << "Atom number mismatch (dcd)" << number_of_atoms << " vs. (pdb)" << atoms.size() << endl; throw; }
#ifdef SASSIM_DEBUG
	clog << "DCD Number of atoms: " << number_of_atoms << endl;
#endif
	// determine byte positions, i.e. read first frame
	DcdFramesetDescriptor dcd_desc;
	dcd_desc.filename = dcdfilename;
	dcd_desc.number_of_frames = dcdheader.number_of_frames;
	dcd_desc.number_of_atoms = number_of_atoms;
	dcd_desc.init_byte_pos = dcdfile.tellg();
	dcd_desc.flag_ext_block1 = dcdheader.flag_ext_block1;
	dcd_desc.flag_ext_block2 = dcdheader.flag_ext_block2;

	if (dcd_desc.flag_ext_block1) {
		dcdfile.read((char*) &marker,sizeof(int32_t));
		dcd_desc.block1_byte_offset = dcdfile.tellg()-dcd_desc.init_byte_pos;
		dcdfile.seekg(marker,ios_base::cur);
		dcdfile.read((char*) &marker,sizeof(int32_t));
	} else {
		dcd_desc.block1_byte_offset = dcdfile.tellg()-dcd_desc.init_byte_pos;
	}

	dcdfile.read((char*) &marker,sizeof(int32_t));
	dcd_desc.x_byte_offset = dcdfile.tellg()-dcd_desc.init_byte_pos;
	dcdfile.seekg(number_of_atoms*sizeof(float),ios_base::cur);
	dcdfile.read((char*) &marker,sizeof(int32_t));

	dcdfile.read((char*) &marker,sizeof(int32_t));
	dcd_desc.y_byte_offset = dcdfile.tellg()-dcd_desc.init_byte_pos;
	dcdfile.seekg(number_of_atoms*sizeof(float),ios_base::cur);
	dcdfile.read((char*) &marker,sizeof(int32_t));

	dcdfile.read((char*) &marker,sizeof(int32_t));
	dcd_desc.z_byte_offset = dcdfile.tellg()-dcd_desc.init_byte_pos;
	dcdfile.seekg(number_of_atoms*sizeof(float),ios_base::cur);
	dcdfile.read((char*) &marker,sizeof(int32_t));

	if (dcd_desc.flag_ext_block2) {
		dcdfile.read((char*) &marker,sizeof(int32_t));
		dcd_desc.block2_byte_offset = dcdfile.tellg()-dcd_desc.init_byte_pos;
		dcdfile.seekg(marker,ios_base::cur);
		dcdfile.read((char*) &marker,sizeof(int32_t));
	} else {
		dcd_desc.block2_byte_offset = dcdfile.tellg()-dcd_desc.init_byte_pos;
	}

	dcd_desc.block_size_byte = dcdfile.tellg()-dcd_desc.init_byte_pos;
	dcd_desc.frame_number_offset = total_frames;
	framesets.push_back(dcd_desc);

	total_frames += dcd_desc.number_of_frames;

	 // end dcd file read out
}


void DCDFrames::read(int frame_number,DcdFrame& blank_frame) {

	// frame_number is absolute, we have to convert it to relative frame_number first
	// and select the right file
//	cout << "read DCD Frame " << frame_number << endl;
	// get the right file descriptor;
	DcdFramesetDescriptor dcd_desc;
	if (!(find_FramesetDescriptor(frame_number,dcd_desc))) {cerr << "no descriptor found for frame number" << frame_number << endl; throw;} ;	
	
	// make frame_number relative;
	frame_number -= dcd_desc.frame_number_offset;

//	cout << "relative frame_number: " << frame_number << endl;	


	// small check: frame_number in range?
	if ( (frame_number < 0) || (frame_number>=dcd_desc.number_of_frames) ) {cerr << "frame_number out of range" << endl; throw;}

	ifstream dcdfile(Settings::get_filepath(dcd_desc.filename).c_str());
	if (dcd_desc.flag_ext_block1) {
		dcdfile.seekg(dcd_desc.init_byte_pos+frame_number*dcd_desc.block_size_byte+dcd_desc.block1_byte_offset,ios_base::beg);
		for (int i=0;i<6;i++) {
			double temp; dcdfile.read((char*) &temp,sizeof(double)); blank_frame.block1.push_back(temp);
//			cout << "dcd unit cell , read: " << temp << endl;
		}
		blank_frame.unit_cell_status = true;
	}
	else {
		blank_frame.unit_cell_status = false;
	}

//	cout << "unit cell read: " << endl;	
	blank_frame.number_of_atoms = dcd_desc.number_of_atoms;

	dcdfile.seekg(dcd_desc.init_byte_pos+frame_number*dcd_desc.block_size_byte+dcd_desc.x_byte_offset,ios_base::beg);
	for (int i=0;i<dcd_desc.number_of_atoms;i++) {
		float temp; dcdfile.read((char*) &temp,sizeof(float)); blank_frame.x.push_back(temp);
	}
	dcdfile.seekg(dcd_desc.init_byte_pos+frame_number*dcd_desc.block_size_byte+dcd_desc.y_byte_offset,ios_base::beg);
	for (int i=0;i<dcd_desc.number_of_atoms;i++) {
		float temp; dcdfile.read((char*) &temp,sizeof(float)); blank_frame.y.push_back(temp);
	}
	dcdfile.seekg(dcd_desc.init_byte_pos+frame_number*dcd_desc.block_size_byte+dcd_desc.z_byte_offset,ios_base::beg);
	for (int i=0;i<dcd_desc.number_of_atoms;i++) {
		float temp; dcdfile.read((char*) &temp,sizeof(float)); blank_frame.z.push_back(temp);
	}

//	cout << "coordinates read: " << endl;	

	if (dcd_desc.flag_ext_block2) {
		dcdfile.seekg(dcd_desc.init_byte_pos+frame_number*dcd_desc.block_size_byte+dcd_desc.block2_byte_offset,ios_base::beg);
		for (int i=0;i<dcd_desc.number_of_atoms;i++) {
			float temp; dcdfile.read((char*) &temp,sizeof(float)); blank_frame.block2.push_back(temp);
		}
	}
	
//	cout << "Dcd frame read: " << endl;	
}


// end of file
