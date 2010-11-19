/*
 *  This file is part of the software sassena
 *
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008-2010 Benjamin Lindner
 *
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

class AtomselectionReader {

public:
    static std::map<std::string,IAtomselection*> read_ndx(std::string filename,std::string selector, std::string expression);
    static IAtomselection* read_pdb(std::string filename,std::string selector, std::string expression);
};


#endif

// end of file
