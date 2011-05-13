/*
 *  This file is part of the software sassena
 *
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008-2010 Benjamin Lindner
 *
 */
 
// direct header
#include "common.hpp"

// standard header
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <sys/time.h>
#include <vector>

// special library headers
#include <boost/asio.hpp>
#include <boost/date_time.hpp>
#include <boost/regex.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/exception/all.hpp>
#include <boost/filesystem.hpp>
#include <boost/mpi.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/program_options.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>

// other headers
#include "log.hpp"
#include "control.hpp"
#include "sample/frames.hpp"
#include "sample/frameset_index.hpp"

#include "SassenaConfig.hpp"

using namespace std;

namespace po = boost::program_options;

void print_title() {
    
	Info::Inst()->write("This software is being developed by Benjamin Lindner.                    ");
	Info::Inst()->write("For help, suggestions or correspondense use:                             ");
	Info::Inst()->write("ben@benlabs.net, Benjamin Lindner (Main Developer, Impl. & Maintenance)  ");
	Info::Inst()->write("franc@cmm.ki.si, Franci Merzel (Methodology)                             ");
	Info::Inst()->write("For publications include the following references:                       ");
	Info::Inst()->write(".........................................................................");
	Info::Inst()->write("1. Sassena - Scattering Calculations on Parallel Computers               ");
	Info::Inst()->write("   to be published                                                       ");		
	Info::Inst()->write(".........................................................................");
    Info::Inst()->write(string("Version Information: ") + string(Sassena_VERSIONSTRING));
	Info::Inst()->write("");

}

void print_description() {
        Info::Inst()->write(".................................................................");
    	Info::Inst()->write("......................D.E.S.C.R.I.P.T.I.O.N......................");
    	Info::Inst()->write(".................................................................");	

    	Info::Inst()->write("This binary generates a trajectory index to allow random access");	
    	Info::Inst()->write("to frames within a molecular dynamics trajectory. ");
    	Info::Inst()->write(".................................................................");	
    	
}

bool init_commandline(int argc,char** argv,po::variables_map& vm) {
    
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("trajectory",po::value<string>(),"name of the trajectory file")
        ("format",po::value<string>(),"format of the trajectory file")
        ("index",po::value<string>(),"name of the trajectory index file (defaults to trajectory file name with txn extension)")          
    ;
    
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    
    
    if (vm.find("help")!=vm.end()) {
        cout << desc << endl;
        return false;
    }
    
    if (vm.find("trajectory")==vm.end()) {
        Err::Inst()->write("Require name of trajectory file");
        cout << desc << endl;
        return false;
    }
        
    if (vm.find("format")==vm.end()) {
        Err::Inst()->write("Require format of trajectory file");
        cout << desc << endl;
        return false;
    }
        
    return true;
}

//void read_parameters(po::variables_map vm) {
//    Info::Inst()->write("Checking command line for parameter overwrite");
    // second stage , allow command line parameters to overwrite defaults in the configuration file
//}

int main(int argc,char* argv[]) {

	// make sure any singleton class exists:
	Info::Inst();
	Err::Inst();
	Warn::Inst();
	
	Info::Inst()->set_prefix(string("Info>>"));
	Warn::Inst()->set_prefix(string("Warn>>"));
	Err::Inst()->set_prefix(string("Err>>"));

    bool initstatus = true;
		
    po::variables_map vm;    
    print_title();
    print_description();
    initstatus = init_commandline(argc,argv,vm);
    if (!initstatus) {
        return 0;
    }

//    read_parameters(vm);

    
    using namespace boost::filesystem;
    path trjpath(vm["trajectory"].as<string>());
    if (!exists(trjpath)) {
        Err::Inst()->write(string("No trajectory found with that path: ")+trjpath.filename().string());
        return -1;
    }

    boost::filesystem::path ipath;
    if (vm.find("index")==vm.end()) {
        std::string tf = vm["trajectory"].as<string>();
        boost::filesystem::path fpath = tf;
    	string fdir;
        if (fpath.parent_path().is_complete()) {
    		fdir = fpath.parent_path().string();             
        } else {
    		fdir = (initial_path() / fpath).string();    	    
    	}
		ipath = (path(fdir).parent_path() / fpath.stem() ).string()+ string(".tnx");
        Warn::Inst()->write("Index file name not specified. Defaulting to filename with tnx extension:");
        Warn::Inst()->write(ipath.string());
    } else {
        ipath = vm["index"].as<string>();        
    }

    if (exists(ipath)) {
        int n=0;
        while (boost::filesystem::exists(ipath.filename().string()+".backup-"+boost::lexical_cast<string>(n))) n++;
        path newidxpath = ipath.filename().string()+".backup-"+boost::lexical_cast<string>(n);
        Warn::Inst()->write(string("Moving old index file to ")+newidxpath.filename().string());        
        boost::filesystem::rename(ipath,newidxpath);
    }

    std::string format = vm["format"].as<string>();

    FileFrameset* p_fs;
    if (format=="dcd") {
        p_fs = new DCDFrameset(trjpath.string(),0);
    } else if (format=="pdb") {
        p_fs = new PDBFrameset(trjpath.string(),0);        
    } else if (format=="xtc") {
        p_fs = new XTCFrameset(trjpath.string(),0);                
    } else if (format=="trr") {
        p_fs = new TRRFrameset(trjpath.string(),0);                        
    } else {
        Err::Inst()->write(string("Format not recognized: ")+format);
        throw;
    }
    Info::Inst()->write("generating index ...");
    p_fs->generate_index();
    Info::Inst()->write(string("writing index to: ")+ipath.string());
    p_fs->save_index(ipath.string());
    
    if (p_fs!=NULL) delete p_fs;
    
    Info::Inst()->write(string("Done."));
    
	return 0;
}


// end of file
