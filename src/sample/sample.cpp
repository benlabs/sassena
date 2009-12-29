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

// other headers
#include "control.hpp"

using namespace std;

Sample::Sample() {
    coordinate_sets.set_atoms(atoms);
}

void Sample::init() {
    	
    Params* params = Params::Inst();
    // create the sample via structure file	
    Info::Inst()->write(string("Reading structure information from file: ")+params->sample.structure.file);
    add_atoms(params->sample.structure.file,params->sample.structure.format);
    Info::Inst()->write(string("Done. Atoms read: ")+to_s(atoms.size()));
    
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

    Params::Inst()->runtime.limits.cache.coordinate_sets = 1;    	
    coordinate_sets.init();
    Info::Inst()->write(string("Total number of coordinate sets found: ")+to_s(coordinate_sets.size()));
	coordinate_sets.set_atoms(atoms);
	coordinate_sets.set_selection(atoms.selections["system"]);
	
	// adjust the coordinate sets cache to the size we need	
    size_t memusage_per_cs = 3*sizeof(double)*atoms.size();
    Params::Inst()->runtime.limits.cache.coordinate_sets = long(Params::Inst()->limits.memory.coordinate_sets/memusage_per_cs);
    Info::Inst()->write(string("Number of cacheable coordinate sets per node: ")+to_s(Params::Inst()->runtime.limits.cache.coordinate_sets));
    
    // fault if not enough memory for one coordinate set
    if (memusage_per_cs>Params::Inst()->limits.memory.coordinate_sets) {
		Err::Inst()->write(string("Insufficient Buffer size for coordinate set."));            
		Err::Inst()->write(string("Size required:")+to_s(memusage_per_cs)+string(" bytes"));            
		Err::Inst()->write(string("Configuration Parameter: limits.memory.coordinate_sets"));
        throw;
    }
}

void Sample::set_kappa(std::string group, double kappa) {
	atoms.assert_selection(group);
	Atomselection& as = atoms.selections[group];

	for (size_t i=0;i<as.indexes.size();i++) {
		atoms[as.indexes[i]].kappa=kappa;		
	}
}

// end of file
