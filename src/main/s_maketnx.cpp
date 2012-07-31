/** \file
This executable routine is used to generate a frame index file which may be required for certain file formats. Sassena requires trajectory files to be seekable (i.e. the frame offset position has to be known without reading the trajectory first). Since some file format do not allow to precompute the offset positions, they have to be indexed first. 

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
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
        ("selection",po::value<string>(),"name of file containing a sub-selection of frames")
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
		std::stringstream ss; ss << "No trajectory found with that path: " << trjpath.filename();
		Err::Inst()->write(ss.str());
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
		std::stringstream ssindex;
        while (true) {
			ssindex.str("");
			ssindex << ipath.filename() << ".backup-0"<< boost::lexical_cast<string>(n);
			if (!boost::filesystem::exists(ssindex.str())) break;
			n++;
		}
        path newidxpath = ssindex.str();
        std::stringstream ss; ss << "Moving old index file to " << newidxpath.filename();
        Warn::Inst()->write(ss.str());
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

	if (vm.find("selection")!=vm.end()) {
		std::vector<size_t> frame_selection;
		ifstream selection_file(vm["selection"].as<string>().c_str());
		size_t framenumber;
		std::string line;
		while (std::getline(selection_file,line)) {
			framenumber = boost::lexical_cast<size_t>(line);
			frame_selection.push_back(framenumber);
		}
		
	    Info::Inst()->write(string("Reducing the frame index to a number of ")+boost::lexical_cast<string>(frame_selection.size()));		
		p_fs->frameset_index_.select(frame_selection);		
	}

    Info::Inst()->write(string("writing index to: ")+ipath.string());
    p_fs->save_index(ipath.string());
    
    if (p_fs!=NULL) delete p_fs;
    
    Info::Inst()->write(string("Done."));
    
	return 0;
}


// end of file
