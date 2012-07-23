/** \file
This file contains a the main sample class which defines structure, coordinates and selections.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/

// direct header
#include "sample/sample.hpp"

// standard header
#include <fstream>
#include <sstream>

// special library headers
#include <boost/lexical_cast.hpp>

// other headers
#include "control.hpp"
#include "log.hpp"
#include "sample/atomselection_reader.hpp"

using namespace std;

Sample::Sample() {
    coordinate_sets.set_atoms(atoms);
}

void Sample::init() {
    	
    Info::Inst()->write("Initializing sample...");
    	
    Params* params = Params::Inst();
    // create the sample via structure file	
    Info::Inst()->write(string("Reading structure from file: ")+params->sample.structure.filepath);
    add_atoms(params->sample.structure.filepath,params->sample.structure.format);
    Info::Inst()->write(string("Done. Atoms read: ")+boost::lexical_cast<string>(atoms.size()));
    
    // add selections / groups
    for(map<string,SampleSelectionParameters*>::iterator sgpi = params->sample.selections.begin();sgpi!=params->sample.selections.end();sgpi++)
    {
    	Info::Inst()->write(string("Initializing selection: ")+sgpi->first);

        std::string type = sgpi->second->type();
        
        if (type=="index") {
            SampleIndexSelectionParameters* p_sampleselection = static_cast<SampleIndexSelectionParameters*>(sgpi->second);
            IndexAtomselection* p_sel = new IndexAtomselection(p_sampleselection->ids_);
            atoms.selections.set(sgpi->first,p_sel);
            Info::Inst()->write(string("Created selection ")+sgpi->first+string(" , elements: ")+boost::lexical_cast<string>(p_sel->size()));
        } else if (type=="range") {
            SampleRangeSelectionParameters* p_sampleselection = static_cast<SampleRangeSelectionParameters*>(sgpi->second);
            RangeAtomselection* p_sel = new RangeAtomselection(p_sampleselection->from_,p_sampleselection->to_);
            atoms.selections.set(sgpi->first,p_sel);            
            Info::Inst()->write(string("Created selection ")+sgpi->first+string(" , elements: ")+boost::lexical_cast<string>(p_sel->size()));
        } else if (type=="lexical") {
            SampleLexicalSelectionParameters* p_sampleselection = static_cast<SampleLexicalSelectionParameters*>(sgpi->second);
            IAtomselection* p_sel = atoms.select(p_sampleselection->expression_);
            atoms.selections.set(sgpi->first,p_sel);                   
            Info::Inst()->write(string("Created selection ")+sgpi->first+string(" , elements: ")+boost::lexical_cast<string>(p_sel->size()));
        } else if (type=="file") {
            SampleFileSelectionParameters* p_sampleselection = static_cast<SampleFileSelectionParameters*>(sgpi->second);

        	if (p_sampleselection->format_=="ndx") {
                std::map<std::string,IAtomselection*> sels = AtomselectionReader::read_ndx(p_sampleselection->filepath_,p_sampleselection->selector_,p_sampleselection->expression_);
            	for(std::map<std::string,IAtomselection*>::iterator i = sels.begin(); i != sels.end(); ++i)
            	{
                    atoms.selections.set(i->first,i->second);
                    Info::Inst()->write(string("Created selection ")+i->first+string(" , elements: ")+boost::lexical_cast<string>(i->second->size()));
            	}
        	} else if (p_sampleselection->format_=="pdb") {
            	IAtomselection* p_sel = AtomselectionReader::read_pdb(p_sampleselection->filepath_,p_sampleselection->selector_,p_sampleselection->expression_);  	    
                atoms.selections.set(sgpi->first,p_sel);
                Info::Inst()->write(string("Created selection ")+sgpi->first+string(" , elements: ")+boost::lexical_cast<string>(p_sel->size()));
        	} else {
                Err::Inst()->write("Selection file format not understood");
                throw;
        	}
        } else {
            Err::Inst()->write(string("Selection type not understood: ")+type);
            throw;
        }
    	
    }
    
    // add custom selections here ( if not already set! )
    if (atoms.selections.exists("system")) {
        Warn::Inst()->write("Renaming selection system to system_RENAMED_BY_SASSENA");
        Warn::Inst()->write("system is a reserved selection word and includes all atoms");
        atoms.selections.rename("system","system_RENAMED_BY_SASSENA");
    }
	// this shortcut creates a full selection
    RangeAtomselection* p_sel = new RangeAtomselection(0,atoms.size()-1);
    atoms.selections.set("system",p_sel); 
    
    coordinate_sets.init();
    Info::Inst()->write(string("Total number of coordinate sets found: ")+boost::lexical_cast<string>(coordinate_sets.size()));
	coordinate_sets.set_atoms(atoms);
	coordinate_sets.set_selection(atoms.selections["system"]);
}

// end of file
