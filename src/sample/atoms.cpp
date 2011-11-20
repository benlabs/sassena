/** \file
This file contains a class which defines an atom group, i.e. a structure.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/
 
// direct header
#include "sample/atoms.hpp"

// standard header
#include <fstream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>

// special library headers
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>

// other headers
#include "sample/atom.hpp"
#include "sample/atomselection.hpp"
#include "math/coor3d.hpp"
#include "control.hpp"
#include "log.hpp"


using namespace std;

Atoms::Atoms(string filename, string fileformat) {
	add(filename,fileformat);
}

Atoms::~Atoms() {}
 

/** Implemenation details:
Only PDB supported at the moment. Assumes the following PDB layout:
Record Format

\verbatim
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
void Atoms::add(string filename, string fileformat) {

	ifstream input(filename.c_str());
	if (input.fail()) {
		Err::Inst()->write(string("Couldn't open structure file: ")+filename);
		Err::Inst()->write("Typo in filename?");
		throw;
	}

	if (fileformat == "pdb") {
		string line;
		while (getline(input,line)) {
			if (line.substr(0,6)=="ATOM  ") {
                string pdbname = line.substr(12,4);
                boost::trim(pdbname);
                string name = Database::Inst()->names.pdb.get(pdbname);
                size_t ID = Database::Inst()->atomIDs.get(name);

                ids_.push_back(ID);
			}
		}		
	}	
}

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


IAtomselection* Atoms::select(string expression) {

    std::vector<size_t> ids;
    
    boost::regex expr(expression);
    
    for(size_t i = 0; i < ids_.size(); ++i)
    {
        string name = Database::Inst()->atomIDs.rget(ids_[i]);
        if ( regex_match(name,expr) ) {
            ids.push_back(ids_[i]);
        }
    }
    
    return (new IndexAtomselection(ids));
}

// end of file
