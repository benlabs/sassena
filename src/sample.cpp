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
#include "settings.hpp"

using namespace std;

void Sample::deuter(std::string group) {
	Atomselection as = atomselections[group];
	for (Atomselection::iterator asi=as.begin();asi!=as.end();asi++) {
		atoms[*asi].name="deuterium";		
	}
}

// does random substitution
void Sample::deuter(std::string group,double coverage,int seed = -1) {
	
	if (seed==-1) seed = time(NULL);
	srand(seed);

	Atomselection as = atomselections[group];
	
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

void Sample::add_selection(std::string name, std::string filename, std::string format,std::string select,double select_value) {
	if (format=="ndx") {
			ifstream ndxfile(Settings::get_filepath(filename).c_str());
			int linecounter = 0;
			string line; 
			map<string,vector<size_t> > indexes;
			int bracketcounter =0;
			string name = "";
			while (getline(ndxfile,line)) {
				size_t pos = line.find("[");
				if ( pos !=string::npos )  {
					size_t pos2 = line.find("]");
					if (pos2==string::npos) {
						cerr << "ERROR>> " << "ndx file is missing closing bracket" << endl;
						throw;
					}
					stringstream cleannamestream(line.substr(pos+1,pos2-pos));
					string cleanname; cleannamestream >> cleanname;
					name = cleanname;
				}
				else if (name!="") {
					stringstream thisline(line);
					size_t index=0;
					while (thisline>>index) {
						indexes[name].push_back(index);
					}
				}					
			}
			if (select!="") {
				if (indexes.find(select)!=indexes.end()) {
					atomselections[select] = Atomselection(atoms,indexes[select],select);									
					clog << "INFO>> " << "Adding selection: " << select << endl;						
				}
			}
			else {
				for (map<string,vector<size_t> >::iterator ii=indexes.begin();ii!=indexes.end();ii++) {
					atomselections[ii->first] = Atomselection(atoms,ii->second,ii->first);									
					clog << "INFO>> " << "Adding selection: " << ii->first << endl;
				}					
			}
	}
	else if (format=="pdb") {
		if (name=="") {
			cerr << "ERROR>> " << "pdb selection entries must specify a name for the selection" << endl;
			throw;
		}
		atomselections[name] = Atomselection(atoms,filename,format,select,select_value,name);			
		clog << "INFO>> " << "Adding selection: " << name << endl;	
		
	}

}

void Sample::add_selection(std::string name,bool select) {
	atomselections[name] = Atomselection(atoms,select,name);	
}

// end of file
