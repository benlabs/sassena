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
#include <boost/date_time.hpp>
#include <boost/regex.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/exception/all.hpp>
#include <boost/filesystem.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/thread.hpp>

// other headers
#include "control.hpp"
#include "log.hpp"
#include "sample/sample.hpp"

#include "SassenaConfig.hpp"

using namespace std;

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

	Info::Inst()->write("This binary computes the scattering intensities directly from");	
	Info::Inst()->write("a molecular dynamics trajectory. ");	
	Info::Inst()->write(".................................................................");	
}



void print_initialization() {
    Info::Inst()->write(".................................................................");
	Info::Inst()->write("...................I.N.I.T.I.A.L.I.Z.A.T.I.O.N...................");
	Info::Inst()->write(".................................................................");	
}


int main(int argc,char* argv[]) {

    // The rank 0 node is responsible for the progress output and to inform the user
    // compute nodes should be silent all the time, except when errors occur.
    // In that case the size and the rank should be included into the error message in the following way:

	// make sure any singleton class exists:
	Info::Inst();
	Err::Inst();
	Warn::Inst();
	
	Info::Inst()->set_prefix(boost::lexical_cast<string>(string(".Info>>")));
	Warn::Inst()->set_prefix(boost::lexical_cast<string>(string(".Warn>>")));
	Err::Inst()->set_prefix(boost::lexical_cast<string>(string(".Err>>")));
	
	Params* params = Params::Inst();
	Database* database = Database::Inst();

	Sample sample;
	
    print_title();
    print_description();
    
    print_initialization();
	    
    params->init(argc,argv);
    database->init();	    
    sample.init();

		
	//------------------------------------------//
	//
	// Communication of the sample
	// At this point it is ILLEGAL to change anything within the sample.
	//
	//------------------------------------------//

Info::Inst()->write("Exchanging sample, database & params information with compute nodes... ");
Info::Inst()->write("params... ");
Info::Inst()->write("database... ");
Info::Inst()->write("sample... ");




// params


// database

std::stringstream paramsstream;
{
    boost::archive::text_oarchive ar(paramsstream);
    ar << *params;
}
{
    boost::archive::text_iarchive ar(paramsstream);
    ar >> *params;    
}

std::stringstream databasestream;
{
    boost::archive::text_oarchive ar(databasestream);
    ar << *database;
}
{
    boost::archive::text_iarchive ar(databasestream);
    ar >> *database;    
}

std::stringstream samplestream;
{
    boost::archive::text_oarchive ar(samplestream);
    ar << sample;
}
{
    boost::archive::text_iarchive ar(samplestream);
    ar >> sample;    
}	

Info::Inst()->write("done");

}