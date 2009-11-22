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
    	atoms.add_selection(sp.name,sp.file,sp.format,sp.select,sp.select_value);
    }
    
    // add custom selections here ( if not already set! )
    if (atoms.selections.find("system")==atoms.selections.end()) {
    	// this shortcut creates a full selection
    	atoms.add_selection("system",true);		
    }
    
    // apply deuteration
    for(size_t i = 0; i < params->sample.deuter.size(); ++i)
    {
    	deuter(params->sample.deuter[i]);
    }

    Params::Inst()->runtime.limits.cache.coordinate_sets = 1;    	
    coordinate_sets.init();
    Info::Inst()->write(string("Total number of coordinate sets found: ")+to_s(coordinate_sets.size()));
	
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

void Sample::deuter(std::string group) {
	atoms.assert_selection(group);
	
	for (Atomselection::iterator asi=atoms.selections[group].begin();asi!=atoms.selections[group].end();asi++) {
		atoms[*asi].name="deuterium";		
	}	
}

// does random substitution
void Sample::deuter(std::string group,double coverage,int seed = -1) {
	
	if (seed==-1) seed = time(NULL);
	srand(seed);

	atoms.assert_selection(group);
	Atomselection& as = atoms.selections[group];
	
	// in order to exactly get the "percentage" given by coverage 100% = 1.0
	// we have to create a temporary list and test afterwards 
	size_t total_size = as.size();
	size_t target_size = int(coverage)*total_size;

	vector<bool> deuterarray(total_size);
	size_t real_size = 0;
	for(size_t i = 0; i < total_size; ++i)
	{
		if (coverage==0.0) {
			deuterarray[i]=false;
		}
		else if (coverage==1.0) {
			deuterarray[i]=true;			
		}
		else if ( coverage > (1.0*rand()/RAND_MAX) ) {
			deuterarray[i]=true;
			real_size++;
		}
		else {
			deuterarray[i]=false;			
		}
	}
	
	Info::Inst()->write(string(" Deuterated with random substitution, seed = ") +to_s(seed) +string(", totalsize/targetsize/realsize=") + to_s(total_size) +string("/") + to_s(target_size) + string("/") + to_s(real_size));
	
	for(size_t i = 0; i < total_size; ++i)
	{
		if (deuterarray[i]) atoms[as[i]].name = "deuterium";
	}
}

void Sample::set_kappa(std::string group, double kappa) {
	atoms.assert_selection(group);
	Atomselection& as = atoms.selections[group];

	for (Atomselection::iterator asi=as.begin();asi!=as.end();asi++) {
		atoms[*asi].kappa=kappa;		
	}
}

// end of file
