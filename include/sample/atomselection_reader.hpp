/** \file
This file contains a class which constructs an atomselection from an input file.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/

#ifndef SAMPLE__ATOMSELECTIONREADER_HPP_
#define SAMPLE__ATOMSELECTIONREADER_HPP_

// common header
#include "common.hpp"

// standard header
#include <string>
#include <vector>

// special library headers
#include <boost/serialization/access.hpp>
#include <boost/serialization/list.hpp>

#include "atomselection.hpp"
#include "atoms.hpp"

/** 
Creates an atomselection from a file.
*/
class AtomselectionReader {

public:
    static std::map<std::string,IAtomselection*> read_ndx(std::string filename,std::string selector, std::string expression);
    static IAtomselection* read_pdb(std::string filename,std::string selector, std::string expression);
};


#endif

// end of file
