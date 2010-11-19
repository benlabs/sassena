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

    	Info::Inst()->write("This binary generates a default configuration file as input for");	
    	Info::Inst()->write("sassena");
    	Info::Inst()->write(".................................................................");	
    	
}

bool init_commandline(int argc,char** argv,po::variables_map& vm) {
    
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("config",po::value<string>()->default_value("scatter.xml"),"name of the configuration file")
    ;
    
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    
    if (vm.find("help")!=vm.end()) {
        cout << desc << endl;
        return false;
    }
    
    if (vm["config"].defaulted()) {
        Info::Inst()->write("No configuration filename specified. Will write to scatter.xml");
    }
    if (boost::filesystem::exists(vm["config"].as<string>())) {
        size_t n=0;
        boost::filesystem::path fp = vm["config"].as<string>();
        std::string fn = fp.string() + ".backup-" + boost::lexical_cast<string>(n);
        while (boost::filesystem::exists(fn)) {
            fn = fp.string() + ".backup-" + boost::lexical_cast<string>(++n);
        }
        Warn::Inst()->write("Moving old configuration file to "+fn);
        boost::filesystem::rename(fp,fn);
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
    
//    Params::Inst()->write_default(vm["config"].as<string>());
    std::ofstream conf(vm["config"].as<string>().c_str());
    
    conf << "<root>" << endl;
    conf << " <sample>" << endl;
    conf << "  <structure>" << endl;
    conf << "   <file>sample.pdb</file>" << endl;    
    conf << "   <format>pdb</format>" << endl;    
    conf << "  </structure>" << endl;
    conf << "  <framesets>" << endl;
    conf << "   <frameset>" << endl;
    conf << "    <file>sample.dcd</file>" << endl;    
    conf << "    <format>dcd</format>" << endl;                
    conf << "   </frameset>" << endl;
    conf << "  </framesets>" << endl;    
    conf << " </sample>" << endl;
    
    conf << " <scattering>" << endl;
    conf << "  <type>all</type>" << endl;
    conf << "  <target>system</target>" << endl;
    conf << "  <vectors>" << endl;
    conf << "   <type>scan</type>" << endl;
    conf << "   <scan>" << endl;
    conf << "    <from>0.1</from>" << endl;
    conf << "    <to>1.5</to>" << endl;    
    conf << "    <points>10</points>" << endl;    
    conf << "    <basevector>" << endl;    
    conf << "     <x>1.0</x>" << endl;    
    conf << "     <y>0.0</y>" << endl;    
    conf << "     <z>0.0</z>" << endl;    
    conf << "    </basevector>" << endl;    
    conf << "   </scan>" << endl;    
    conf << "  </vectors>" << endl;
    conf << "  <average>" << endl;
    conf << "   <orientation>" << endl;
    conf << "    <type>vectors</type>" << endl;
    conf << "    <vectors>" << endl;
    conf << "     <type>sphere</type>" << endl;
    conf << "     <algorithm>boost_uniform_on_sphere</algorithm>" << endl;
    conf << "     <resolution>50</resolution>" << endl;
    conf << "    </vectors>" << endl;
    conf << "   </orientation>" << endl;        
    conf << "  </average>" << endl;
    conf << " </scattering>" << endl;
    
    conf << "</root>" << endl;
    
    conf.close();
    Info::Inst()->write(string("Configuration file written to ")+vm["config"].as<string>());
    Info::Inst()->write("Have a nice day!");
    
    
	return 0;
}


// end of file
