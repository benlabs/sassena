/*
 *  atoms.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */
// direct header
#include "sample/atoms.hpp"

// standard header
#include <fstream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>

// other headers
#include "sample/atom.hpp"
#include "sample/atomselection.hpp"
#include "math/coor3d.hpp"
#include "control.hpp"
#include "log.hpp"


using namespace std;


// for the pdb format (ATOM entry look below)
Atoms::Atoms(string filename, string fileformat) {
	add(filename,fileformat);
}
 
// for the pdb format (ATOM entry look below)
void Atoms::add(string filename, string fileformat) {

	ifstream input(filename.c_str());
	if (input.fail()) {
		Err::Inst()->write(string("Couldn't open structure file: ")+filename);
		Err::Inst()->write("Typo in filename?");
		throw;
	}
	int atom_index=this->size();

	if (fileformat == "pdb") {
		string line;
		Atom temp_atom;
		int line_counter=0; int atom_counter =0;
		// lookup table for speedup:
		while (getline(input,line)) {
			if (line.substr(0,6)=="ATOM  ") {
				try {
					temp_atom.name = Database::Inst()->names.pdb.get(line.substr(12,4));
					temp_atom.ID = Database::Inst()->atomIDs.get(temp_atom.name); // get a unique ID for the name
					
					// independent atomic constants:
					temp_atom.mass = Database::Inst()->masses.get(temp_atom.ID);
					
//					temp_atom.original_name = line.substr(12,4);
//					temp_atom.residue_name = line.substr(17,4);
//					temp_atom.chainid = line.substr(21,1);
//					temp_atom.resseq = line.substr(22,4);
//					temp_atom.occupancy = line.substr(54,6);
//					temp_atom.tempFactor = line.substr(60,6);
//					temp_atom.element = line.substr(76,2);
//					temp_atom.charge = line.substr(78,2);																															
//					temp_atom.x = atof(line.substr(30,8).c_str());
//					temp_atom.y = atof(line.substr(38,8).c_str());
//					temp_atom.z = atof(line.substr(46,8).c_str());
//					temp_atom.beta = atof(line.substr(60,65).c_str());
					// make sure kappa is initialized:
				} catch (...) { 
					Err::Inst()->write(string("Error at reading pdb line:"));
					Err::Inst()->write(line);
					throw; 
				}
				temp_atom.index = atom_index++; // assign and THEN increment...
				push_back(temp_atom);
				atom_counter++;
			}
			line_counter++;
		}		
	}	
	
    system_selection = select();
}

/* pdb file format , atom section:
Record Format

COLUMNS        DATA  TYPE    FIELD        DEFINITION
-------------------------------------------------------------------------------------
 1 -  6        Record name   "ATOM  "
 7 - 11        Integer       serial       Atom  serial number.
13 - 16        Atom          name         Atom name.
17             Character     altLoc       Alternate location indicator.
18 - 20        Residue name  resName      Residue name.
22             Character     chainID      Chain identifier.
23 - 26        Integer       resSeq       Residue sequence number.
27             AChar         iCode        Code for insertion of residues.
31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
55 - 60        Real(6.2)     occupancy    Occupancy.
61 - 66        Real(6.2)     tempFactor   Temperature  factor.
77 - 78        LString(2)    element      Element symbol, right-justified.
79 - 80        LString(2)    charge       Charge  on the atom.
*/

//void Atoms::write(string filename,Frame& frame, string fileformat) {
//
//	ofstream ofile(filename.c_str());
//
//	if (fileformat == "pdb") {
//		string line;
//		Atom temp_atom;
//		int line_counter=0; int atom_counter =0;
//		// lookup table for speedup:
//		map<string,string> quicklookup;
//		
//		for(Atoms::iterator ati=begin();ati!=end();ati++) {
//			ofile << setw(6) << "ATOM  ";           // table entry identifier
//			if (ati->index<100000)
//				ofile << setw(5) << ati->index;         // Atom serial number
//			else
//				ofile << "*****";         // Atom serial number			
//			ofile << setw(1) << " ";                // --not used--
//			ofile << setw(4) << ati->original_name; // Atom name.
//			ofile << setw(1) << " ";                // chain identifier
//			ofile << setw(4) << ati->residue_name;  // residue name
//			ofile << setw(1) << ati->chainid;       // Chain identifier
//			ofile << setw(4) << ati->resseq;        // Residue sequence number
//			ofile << setw(1) << " ";                // Code for insertion of residues
//			ofile << setw(3) << "   ";              // --not used--
//			ofile << setw(8) << setiosflags(ios::fixed) << setprecision(3) << frame.x[ati->index];             // Orthogonal coordinates for X in Angstroms
//			ofile << setw(8) << setiosflags(ios::fixed) << setprecision(3) << frame.y[ati->index];             // Orthogonal coordinates for Y in Angstroms
//			ofile << setw(8) << setiosflags(ios::fixed) << setprecision(3) << frame.z[ati->index];             // Orthogonal coordinates for Z in Angstroms
//			ofile << setw(6) << ati->occupancy;     // Occupancy
//			ofile << setw(6) << ati->tempFactor;    // Temperature  factor
//			ofile << setw(10)<< "          ";       // --not used--
//			ofile << setw(2) << ati->element;       // Element symbol, right-justified
//			ofile << setw(2) << ati->charge;        // Charge  on the atom
//			ofile << endl;
//			atom_counter++;
//			line_counter++;
//		}
//	}			
//}


Atomselection Atoms::select(string label) {
    Atomselection selection;

    for(size_t i = 0; i < this->size(); ++i)
    {
        if (label=="") { // default: add
            selection.add(this->at(i).index);
        } else if ( this->at(i).name ==label) {
            selection.add(this->at(i).index);
        }
    }
    
    return selection;
}

Atomselection Atoms::select_pdb(std::string filename, std::string select, double select_value) {
    Atomselection selection;

	ifstream pdbfile(filename.c_str());
	size_t linecounter = 0;
	if (select=="beta") {
		string line; double beta;
		while (getline(pdbfile,line)) {
			if (line.substr(0,6)=="ATOM  ") {
				try {
					beta = atof(line.substr(60,65).c_str());
					if (beta==select_value) {
						selection.add(linecounter);
					}
					
				} catch (...) { cerr << "Error at reading atomselection line:" << endl << line << endl; throw; }
				linecounter++;				
			}
		}
	} else {
        Err::Inst()->write("atomselection only supports beta at the moment");
        throw;
	}

    return selection;
}

vector<pair<string,Atomselection> > Atoms::select_ndx(std::string filename) {
    vector<pair<string,Atomselection> > result;
    
   	ifstream ndxfile(filename.c_str());
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
   				Err::Inst()->write("ndx file is missing closing bracket");					
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
   				indexes[name].push_back(index-1); // ndx files start with 1 as an index
   			}
   		}
   	}
   	
   	for (map<string,vector<size_t> >::iterator ii=indexes.begin();ii!=indexes.end();ii++) {
   		Atomselection selection;
   		for(size_t i = 0; i < ii->second.size(); ++i)
   		{
               selection.add(ii->second[i]);
   		}
        result.push_back(make_pair<string,Atomselection>(ii->first,selection));
   		Info::Inst()->write(string("Adding selection: ") + ii->first);
   	}

    return result;
}

void Atoms::push_selection(std::string name, Atomselection& selection) {
    selections[name] = selection;
}


void Atoms::assert_selection(std::string groupname) {
	if (selections.find(groupname)==selections.end()) {
		Err::Inst()->write(string("Assertion of selection: '")+groupname+string("' failed"));
		throw;
	}
}

void Atoms::clear_selections() {
    selections.clear();
}
// end of file
