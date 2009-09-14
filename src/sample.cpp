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
#include "sample.hpp"

// standard header
#include <fstream>
#include <sstream>

// special library headers

// other headers
#include "log.hpp"

using namespace std;

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
	
	clog << "INFO>> " << " Deuterated with random substitution, seed = " << seed << ", totalsize/targetsize/realsize=" << total_size << "/" << target_size << "/" << real_size << endl;
	
	for(size_t i = 0; i < total_size; ++i)
	{
		if (deuterarray[i]) atoms[as[i]].name = "deuterium";
	}
}

// end of file
