/*
 *  atomselection.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

// direct header
#include "atomselection.hpp"

// standard header
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>

// special library headers
#include <boost/regex.hpp>

// other headers
#include "atom.hpp"
#include "atoms.hpp"
#include "settings.hpp"

using namespace std;

Atomselection::Atomselection(Atoms& atoms,bool select,std::string aname) {
	if (select) {
		for(Atoms::iterator ai=atoms.begin();ai!=atoms.end();ai++) {
			push_back(ai->index);
		}
		booleanarray.resize(atoms.size(),true);		
	}
	else {
		booleanarray.resize(atoms.size(),false);		
	}
	name = aname;
}

Atomselection::Atomselection(Atoms& atoms,std::string filename,std::string format,std::string select,double select_value,std::string aname) {
	if (format=="pdb") {
		ifstream pdbfile(Settings::get_filepath(filename).c_str());
		int linecounter = 0;
		if (select=="beta") {
			string line; double beta;
			while (getline(pdbfile,line)) {
				if (line.substr(0,6)=="ATOM  ") {
					try {
						beta = atof(line.substr(60,65).c_str());
						if (beta==select_value) {
							push_back(atoms[linecounter].index);
						}
						
					} catch (...) { cerr << "Error at reading atomselection line:" << endl << line << endl; throw; }
					linecounter++;				
				}
			}
		}
		else {
			throw("atomselection only supports beta at the moment");
		}
	}
	
	booleanarray.resize(atoms.size(),false);
	for (Atoms::iterator ai=atoms.begin();ai!=atoms.end();ai++) {
		booleanarray[ ai->index ] = true;
	}
	name = aname;
}


Atomselection::Atomselection(Atoms& atoms,std::string filename,std::string format,std::string aname) {
	if (format=="pdb") {
		ifstream pdbfile(Settings::get_filepath(filename).c_str());
		int linecounter = 0;
		string line; 
		while (getline(pdbfile,line)) {
			if (line.substr(0,6)=="ATOM  ") {
			try {
				push_back(atoms[linecounter].index);
			} catch (...) { cerr << "Error at reading atomselection line:" << endl << line << endl; throw; }
			linecounter++;	
			}			
		}
	}
	
	booleanarray.resize(atoms.size(),false);
	for (Atoms::iterator ai=atoms.begin();ai!=atoms.end();ai++) {
		booleanarray[ ai->index ] = true;
	}
	name = aname;
}

void Atomselection::add(const Atom& atom) {
	push_back(atom.index);
	booleanarray[atom.index]=true;
}

void Atomselection::remove(const Atom& atom) {
	iterator ai = find(begin(),end(),atom.index);
	erase(ai);
	booleanarray[atom.index]=false;
}

// end of file
