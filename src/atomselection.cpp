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
#include "sample.hpp"
#include "settings.hpp"

using namespace std;

	// selects by range
Atomselection::Atomselection(Atoms::iterator f,Atoms::iterator l) {
		for(Atoms::iterator ai=f;ai!=l;ai++) {
			push_back(ai->index);
		}
	}
	
	// selects by index+length
Atomselection::Atomselection(Atoms::iterator f,int length) {
	int i=0;
		for (Atoms::iterator ai=f;i<length;ai++,i++) {
			push_back(ai->index);
		}
	}
	// select my "atomname" match
Atomselection::Atomselection(Atoms::iterator f,Atoms::iterator l, std::string match) {
		boost::regex e(match);
		for(Atoms::iterator ai=f;ai!=l;ai++) {
			if (regex_match(ai->name,e)) {
				push_back(ai->index);
			}
		}
	}
	// select whole set of atoms:
//Atomselection::Atomselection( std::string match) : sample(&s) {
//		boost::regex e(match);
//		for(Atoms::iterator ai=s.atoms.begin();ai!=s.atoms.end();ai++) {
//			if (regex_match(ai->name,e)) {
//				push_back(&(*ai));
//			}
//		}
//	}
	
Atomselection::Atomselection(Atoms& atoms) {
		for(Atoms::iterator ai=atoms.begin();ai!=atoms.end();ai++) {
			push_back(ai->index);
		}
	}
	
void Atomselection::add(Atom& atom) {
	push_back(atom.index);
}

Atomselection::Atomselection(Atoms& atoms,std::string filename,std::string format,std::string select,double select_value) {
	if (format=="pdb") {
		ifstream pdbfile(Settings::get_filepath(filename).c_str());
		int linecounter = 0;
		if (select=="beta") {
			string line; double beta;
			while (getline(pdbfile,line)) {
				if (line.substr(0,6)=="ATOM  ") {
					try {
						beta = atof(line.substr(60,65).c_str());
						if (beta==select_value) push_back(atoms[linecounter].index);
					} catch (...) { cerr << "Error at reading atomselection line:" << endl << line << endl; throw; }
					linecounter++;				
				}
			}
		}
		else {
			throw("atomselection only supports beta at the moment");
		}
	}
}


Atomselection::Atomselection(Atoms& atoms,std::string filename,std::string format) {
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
}


void Atomselection::readpdb(Atoms& atoms,std::string pdbname,PDBSELECT select) {
	ifstream pdbfile(Settings::get_filepath(pdbname).c_str());
	int linecounter = 0;
	if (select==BETA) {
		string line; double beta;
		while (getline(pdbfile,line)) {
			if (line.substr(0,6)=="ATOM  ") {
				try {
					beta = atof(line.substr(60,65).c_str());
					if (beta==1.0) push_back(atoms[linecounter].index);
				} catch (...) { cerr << "Error at reading atomselection line:" << endl << line << endl; throw; }
				linecounter++;				
			}
		}
	}
	else {
		throw("atomselection only supports beta at the moment");
	}
}	

//void Atomselection::truncate_cylinder(Sample& sample, CartesianCoor3D origin,CartesianCoor3D zaxis,double radius) {
//	for (Atomselection::iterator asi=begin();asi!=end();asi++) {
//		CylinderCoor3D c = sample->currentframe.coord3D((*asi)->index);
//		if ( c.r > radius ) {
//			erase(asi);
//		}
//	}
//}
//
//void Atomselection::truncate_cshell(Sample& sample, CartesianCoor3D origin,CartesianCoor3D zaxis,double radiusinner,double radiusouter) {
//	for(Atomselection::iterator asi=begin();asi!=end();asi++) {
//		CylinderCoor3D c = sample->currentframe.coord3D((*asi)->index);
//		if ((c.r< radiusinner) || (c.r > radiusouter) ) {
//			erase(asi);
//		}
//	}
//}
//

// end of file
