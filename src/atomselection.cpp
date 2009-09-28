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
#include "decompose.hpp"
#include "log.hpp"
#include "parameters.hpp"
#include "database.hpp"

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

Atomselection::Atomselection(Atoms& atoms,std::vector<size_t> indexes,std::string aname) {
	booleanarray.resize(atoms.size(),false);					
			
	for(std::vector<size_t>::iterator ii=indexes.begin();ii!=indexes.end();ii++) {
		if ((*ii<0) || (*ii>=booleanarray.size())) {
			Err::Inst()->write(string("Index out of bound for Atomselection ")+aname);
			Err::Inst()->write(string("Index was ")+to_s(*ii)+string(", size of atom table: ")+to_s(atoms.size()));			
			throw;
		}
			push_back(*ii);		
			booleanarray.at(*ii) = true;
	}
	name = aname;
}


Atomselection::Atomselection(Atoms& atoms,std::string filename,std::string format,std::string select,double select_value,std::string aname) {
	if (format=="pdb") {
		ifstream pdbfile(filename.c_str());
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
	for (iterator ai=this->begin();ai!=this->end();ai++) {
		booleanarray[ *ai ] = true;
	}
	name = aname;
}


Atomselection::Atomselection(Atoms& atoms,std::string filename,std::string format,std::string aname) {
	if (format=="pdb") {
		ifstream pdbfile(filename.c_str());
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
	for (iterator ai=this->begin();ai!=this->end();ai++) {
		booleanarray[ *ai ] = true;
	}
	name = aname;
}


Atomselection& Atomselection::operator+=(const Atomselection& that) {
	if (that.size()==0) return *this;
	// append all elements from 2nd selection
	insert(this->end(),that.begin(),that.end());
	// enforce a sort...
	sort(this->begin(),this->end());
	// to allow removal of duplicates
	erase(unique(this->begin(),this->end()),this->end());
	// finally set missing boolean values
	for(size_t i = 0; i < that.size(); ++i)
	{
		booleanarray[that[i]]=true;
	}
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

std::vector<Atomselection> Atomselection::slice(size_t number) {
	std::vector<Atomselection> result;
	
	EvenDecompose ed(this->size(),number);

	size_t offset =0;
	for(EvenDecompose::iterator edi = ed.begin(); edi != ed.end(); ++edi)
	{
		result.push_back(subset(offset,*edi));
		offset += *edi;
	}
	return result;
}

Atomselection Atomselection::subset(size_t offset,size_t maxcount) {
	Atomselection result;

	result.booleanarray.resize(this->booleanarray.size(),false); // conserve the original mapping
	
	for(size_t i = offset; i < (offset+maxcount); ++i)
	{
		if ((i>=this->size())) break;
		
		result.push_back(this->at(i));
		result.booleanarray[this->at(i)]=true;
	}
	
	return result;
}

bool Atomselection::is_subset_of(Atomselection& other) {
	// use the boolean-array to test:
	// if this->booleanarray has a 1 where other doesn't , then return false
	std::vector<bool>& otherbool = other.booleanarray;
	for(size_t i = 0; i < booleanarray.size(); ++i)
	{
		if ( (booleanarray[i]==true) && (otherbool[i]==false) ) return false;
	}
	
	return true;
};

// end of file
