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
#include "atoms.hpp"

// standard header
#include <fstream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>

// special library headers
#include <boost/regex.hpp>

// other headers
#include "atom.hpp"
#include "coor3d.hpp"
#include "settings.hpp"

using namespace std;


// for the pdb format (ATOM entry look below)
Atoms::Atoms(string filename, string fileformat) {
	add(filename,fileformat);
}

// for the pdb format (ATOM entry look below)
void Atoms::add(string filename, string fileformat) {

	ifstream input(Settings::get_filepath(filename).c_str());

	int atom_index=this->size();

	if (fileformat == "pdb") {
		string line;
		Atom temp_atom;
		int line_counter=0; int atom_counter =0;
		// lookup table for speedup:
		map<string,string> quicklookup;
		while (getline(input,line)) {
			if (line.substr(0,6)=="ATOM  ") {
				try {
					temp_atom.name = guess_atomname(line.substr(12,4),fileformat,quicklookup);
					temp_atom.original_name = line.substr(12,4);
					temp_atom.residue_name = line.substr(17,4);
					temp_atom.chainid = line.substr(21,1);
					temp_atom.resseq = line.substr(22,4);
					temp_atom.occupancy = line.substr(54,6);
					temp_atom.tempFactor = line.substr(60,6);
					temp_atom.element = line.substr(76,2);
					temp_atom.charge = line.substr(78,2);	
																														
					temp_atom.x = atof(line.substr(30,8).c_str());
					temp_atom.y = atof(line.substr(38,8).c_str());
					temp_atom.z = atof(line.substr(46,8).c_str());
					temp_atom.beta = atof(line.substr(60,65).c_str());
					temp_atom.mass = Settings::get("atom_masses")[temp_atom.name];	
					// make sure kappa is initialized:
					temp_atom.kappa = 1.0;	
					temp_atom.volume = Settings::get_volume(temp_atom.name);
				} catch (...) { cerr << "Error at reading pdb line:" << endl << line << endl; throw; }
				temp_atom.index = atom_index++; // assign and THEN increment...
				push_back(temp_atom);
				atom_counter++;
			}
			line_counter++;
		}		
	
	#ifdef DEBUG
		clog << "Finished reading pdb file" << endl;
		clog << "Lines read: " << line_counter << ", Atoms read: " << atom_counter << endl;
	#endif
	}	
}

void Atoms::read_deuter(ifstream& input,string fileformat) {
	
	string line; double beta;
	int atom_counter =0;
	while (getline(input,line)) {
		if (line.substr(0,6)=="ATOM  ") {
			try {
				beta = atof(line.substr(60,65).c_str());
				if (beta==1.0) {
					(*this)[atom_counter].name="deuterium";
				}

				atom_counter++;
			} catch (...) { cerr << "Error at reading deuterpdb line:" << endl << line << endl; throw; }
		}
	}
}

void Atoms::read_particle(ifstream& input, string fileformat) {
	
	string line; double beta;
	int atom_counter =0;
	while (getline(input,line)) {
		if (line.substr(0,6)=="ATOM  ") {
			try {
				beta = atof(line.substr(60,65).c_str());
				if (beta==1.0) (*this)[atom_counter].particle=true; else (*this)[atom_counter].particle=false;
				atom_counter++;
			} catch (...) { cerr << "Error at reading particlepdb line:" << endl << line << endl; throw; }
		}
	}
}


void Atoms::read_solvent(ifstream& input, string fileformat) {
	
	string line; double beta;
	int atom_counter =0;
	while (getline(input,line)) {
		if (line.substr(0,6)=="ATOM  ") {
			try {
				beta = atof(line.substr(60,65).c_str());
				if (beta==1.0) (*this)[atom_counter].solvent=true; else (*this)[atom_counter].solvent=false;
				atom_counter++;
			} catch (...) { cerr << "Error at reading solventpdb line:" << endl << line << endl; throw; }
		}
	}
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

void Atoms::print_statistics() {
	clog << "Size of Atom-Vector: " << (*this).size() << endl;
	map<string,int> counts;
	for (vector<Atom>::iterator atom_i=(*this).begin();atom_i!=(*this).end();atom_i++) {
		if (counts.find(atom_i->name)!=counts.end()) counts[atom_i->name]++; else counts[atom_i->name]=1;
	}
	for (map<string,int>::iterator counts_i=counts.begin();counts_i!=counts.end();counts_i++) {
		clog << counts_i->first << " counted " << counts_i->second << " times." << endl;
	}
	int particle=0; int solvent=0;
	for (Atoms::iterator ai=begin();ai!=end();ai++) {
		if (ai->particle) particle++;
		if (ai->solvent) solvent++;
	}
	cout << "Stats for the particle/solvent system:" << endl;
	cout << "Particle atoms: " << particle<< endl;
	cout << "Solvent atoms: " << solvent<< endl;
}


/* private helper function: */
string Atoms::guess_atomname(const string atomname,string fileformat,map<string,string>& quicklookup) {	
	// if pdbatomname exactly matches an entry in quicklookup, take it
	// otherwise do a regular expression analysis and add result to quicklookup
	
	libconfig::Setting& s = Settings::get("atom_names")[fileformat];
	if (quicklookup.find(atomname)!=quicklookup.end()) {
		return quicklookup[atomname];
	}
	else
	{
		int len = s.getLength();
		int matches = 0;
		string match = "<Unknown>";
		for (int i=0;i<len;i++) {
			boost::regex e(s[i]);
			if (regex_match(atomname,e)) {
				if (matches>0) {
					cerr <<  "ERROR>> " <<"PDB atom name matches " << match << " and " << s[i].getName() << endl;
					cerr <<  "ERROR>> " <<"expression to match: '" << atomname << "'" << endl;
					throw;
					} else { match = s[i].getName(); matches++; }
			}
		}
		if (matches==0) {
			cerr << "ERROR>> " << "PDB atom name not recognized: <" << atomname  << "> " << endl;
			throw;
		}
		quicklookup[atomname]=match;

		return match;
	}
}

void Atoms::write(string filename,Frame& frame, string fileformat) {

	ofstream ofile(filename.c_str());

	if (fileformat == "pdb") {
		string line;
		Atom temp_atom;
		int line_counter=0; int atom_counter =0;
		// lookup table for speedup:
		map<string,string> quicklookup;
		
		for(Atoms::iterator ati=begin();ati!=end();ati++) {
			ofile << setw(6) << "ATOM  ";           // table entry identifier
			if (ati->index<100000)
				ofile << setw(5) << ati->index;         // Atom serial number
			else
				ofile << "*****";         // Atom serial number			
			ofile << setw(1) << " ";                // --not used--
			ofile << setw(4) << ati->original_name; // Atom name.
			ofile << setw(1) << " ";                // chain identifier
			ofile << setw(4) << ati->residue_name;  // residue name
			ofile << setw(1) << ati->chainid;       // Chain identifier
			ofile << setw(4) << ati->resseq;        // Residue sequence number
			ofile << setw(1) << " ";                // Code for insertion of residues
			ofile << setw(3) << "   ";              // --not used--
			ofile << setw(8) << setiosflags(ios::fixed) << setprecision(3) << frame.x[ati->index];             // Orthogonal coordinates for X in Angstroms
			ofile << setw(8) << setiosflags(ios::fixed) << setprecision(3) << frame.y[ati->index];             // Orthogonal coordinates for Y in Angstroms
			ofile << setw(8) << setiosflags(ios::fixed) << setprecision(3) << frame.z[ati->index];             // Orthogonal coordinates for Z in Angstroms
			ofile << setw(6) << ati->occupancy;     // Occupancy
			ofile << setw(6) << ati->tempFactor;    // Temperature  factor
			ofile << setw(10)<< "          ";       // --not used--
			ofile << setw(2) << ati->element;       // Element symbol, right-justified
			ofile << setw(2) << ati->charge;        // Charge  on the atom
			ofile << endl;
			atom_counter++;
			line_counter++;
		}		
	
	#ifdef DEBUG
		clog << "Finished writing pdb file" << endl;
		clog << "Lines written: " << line_counter << ", Atoms written: " << atom_counter << endl;
	#endif
	}	
}

// end of file
