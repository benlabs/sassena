/** \file
This file contains a class which provides a management interface for atomselections.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/

#ifndef SAMPLE__ATOMSELECTIONS_HPP_
#define SAMPLE__ATOMSELECTIONS_HPP_

// common header
#include "atomselection.hpp"

// standard header
#include <string>
#include <map>

// special library headers
#include <boost/serialization/access.hpp>
#include <boost/serialization/list.hpp>

/** 
Manages a set of atom selections and their labels.
*/
class Atomselections  {
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        ar.register_type(static_cast<IndexAtomselection*>(NULL));
        ar.register_type(static_cast<RangeAtomselection*>(NULL));
        
        ar & _selections;
    }
	///////////////////

    std::map<std::string,IAtomselection*> _selections;
    
    typedef std::map<std::string,IAtomselection*>::iterator iterator;
public:	
    ~Atomselections();
    
    IAtomselection* operator[](std::string name);
    void set(std::string name,IAtomselection* selection);
    bool exists(std::string name);
    void rename(std::string name,std::string newname);
    
};

#endif

// end of file