/*
 *  sample.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

// direct header
#include "sample/sample.hpp"

// standard header
#include <fstream>
#include <sstream>

// special library headers
#include <boost/lexical_cast.hpp>

// other headers
#include "control.hpp"
#include "log.hpp"

using namespace std;

Sample::Sample() {
    coordinate_sets.set_atoms(atoms);
}

void Sample::init() {
    	
    Info::Inst()->write("Initializing sample...");
    	
    Params* params = Params::Inst();
    // create the sample via structure file	
    Info::Inst()->write(string("Reading structure from file: ")+params->sample.structure.file);
    add_atoms(params->sample.structure.file,params->sample.structure.format);
    Info::Inst()->write(string("Done. Atoms read: ")+boost::lexical_cast<string>(atoms.size()));
    
    // add selections / groups
    for(map<string,SampleGroupParameters>::iterator sgpi = params->sample.groups.begin();sgpi!=params->sample.groups.end();sgpi++)
    {
    	SampleGroupParameters& sp = sgpi->second;
    	Info::Inst()->write(string("Reading selection file: ")+sp.file);
    	if (sp.format=="ndx") {
        	vector<pair<string,Atomselection> > sels = atoms.select_ndx(sp.file);    	    
        	for(size_t i = 0; i < sels.size(); ++i)
        	{
                atoms.push_selection(sels[i].first,sels[i].second);
        	}
    	} else if (sp.format=="pdb") {
        	Atomselection sel = atoms.select_pdb(sp.file,sp.select,sp.select_value);    	    
            atoms.push_selection(sp.name,sel);
    	} else {
            Err::Inst()->write("Selection file format not understood");
            throw;
    	}
    }
    
    // add custom selections here ( if not already set! )
    if (atoms.selections.find("system")!=atoms.selections.end()) {
        Warn::Inst()->write("Renaming selection system to system_RENAMED_BY_SASSENA");
        Warn::Inst()->write("system is a reserved selection word and includes all atoms");
        map<string,Atomselection>::iterator si = atoms.selections.find("system");
        atoms.selections["system_RENAMED_BY_SASSENA"]=si->second;
        atoms.selections.erase(si);
    }
	// this shortcut creates a full selection
	Atomselection sel = atoms.select();		
    atoms.push_selection("system",sel); // provide two alternative ways to access a full selection
    atoms.system_selection = sel;

    coordinate_sets.init();
    Info::Inst()->write(string("Total number of coordinate sets found: ")+boost::lexical_cast<string>(coordinate_sets.size()));
	coordinate_sets.set_atoms(atoms);
	coordinate_sets.set_selection(atoms.selections["system"]);
}

// end of file
