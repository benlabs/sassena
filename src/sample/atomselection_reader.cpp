/** \file
This file contains a class which constructs an atomselection from an input file.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/
 
// direct header
#include "sample/atomselection_reader.hpp"

// standard header
#include <fstream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>

// special library headers
#include "boost/regex.hpp"
#include <boost/algorithm/string.hpp>
 
// other headers
#include "sample/atoms.hpp"
#include "sample/atomselection.hpp"
#include "math/coor3d.hpp"
#include "control.hpp"
#include "log.hpp"


using namespace std;

std::map<std::string,IAtomselection*> AtomselectionReader::read_ndx(std::string filename,std::string selector, std::string expression) {

   	ifstream ndxfile(filename.c_str());
   	string line; 
   	map<string,vector<size_t> > indexes;
   	string name = "";
   	
   	boost::regex expr(expression);
   	
   	if (selector=="name") {
   	    // do nothing
    } else {
            Err::Inst()->write("NDX Atomselection only supports selection by 'name' at the moment");
            throw;
    }
   	
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
            boost::trim(name);
   		}
   		else if (name!="") {
   		    if (regex_match(name,expr)) {
       			stringstream thisline(line);
       			size_t index=0;
       			while (thisline>>index) {
       				indexes[name].push_back(index-1); // ndx files start with 1 as an index
       			}   
   		    }
   		}
   	}
   	
    map<string,IAtomselection*> result;
   	   	
   	for (map<string,vector<size_t> >::iterator ii=indexes.begin();ii!=indexes.end();ii++) {
        result[ii->first] = new IndexAtomselection(ii->second);
   	}

    return result;
}

IAtomselection* AtomselectionReader::read_pdb(std::string filename,std::string selector, std::string expression) {

    std::vector<size_t> ids;

    boost::regex expr(expression);

   	if (selector=="beta") {
   	    // do nothing
    } else if (selector=="segid") {
   	    // do nothing        
    } else {
            Err::Inst()->write("PDB Atomselection only supports beta and segid at the moment");
            throw;
    }

	ifstream pdbfile(filename.c_str());
	size_t linecounter = 0;
	string line; 
	
	while (getline(pdbfile,line)) {
		if (line.substr(0,6)=="ATOM  ") {
		    if (selector=="beta") {
    		    std::string value = line.substr(60,6);	
                boost::trim(value);		        
    		    if (regex_match(value,expr)) ids.push_back(linecounter);		        
		    }
		    if (selector=="segid") {
    		    std::string value = line.substr(72,4);			        
                boost::trim(value);
    		    if (regex_match(value,expr)) ids.push_back(linecounter);		        		        
		    }
		    linecounter++;				
		}
    }

    return (new IndexAtomselection(ids) );
}

// end of file
