/** \file
This file contains a class which provides a management interface for atomselections.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/

// direct header
#include "sample/atomselections.hpp"

// standard header
#include <string>
#include <map>

// special library headers
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

// other headers
#include "log.hpp"


using namespace std;
IAtomselection* Atomselections::operator[](std::string name) {
    Atomselections::iterator it = _selections.find(name);
    if (it==_selections.end()) {
        Err::Inst()->write(string("Selection with that name has not been defined: ")+name);
            Err::Inst()->write(string("Add a selection with the above name to your configuration file"));
        throw;
    }
    return it->second;
}

void Atomselections::set(std::string name,IAtomselection* selection) {
    if (_selections.find(name)!=_selections.end()) {
        Err::Inst()->write(string("Selection with that name already exists: ")+name);
        Err::Inst()->write(string("Check your configuration file"));
        throw;
    }
    _selections[name]=selection;
}

bool Atomselections::exists(std::string name) {
    if (_selections.find(name)!=_selections.end()) {
        return true;
    } else {
        return false;
    }
}

void Atomselections::rename(std::string name,std::string newname) {
    Atomselections::iterator it = _selections.find(name);
    if (it==_selections.end()) {
        Err::Inst()->write(string("Cannot rename selection with that name which not been defined: ")+name);
        Err::Inst()->write(string("You shouldn't see this message. File a bug report."));
        throw;
    }
    if (_selections.find(newname)!=_selections.end()) {
        Err::Inst()->write(string("Cannot rename selection to a name which has already been defined: ")+newname);
        Err::Inst()->write(string("You shouldn't see this message. File a bug report."));
        throw;
    }

    _selections[newname]=it->second;
    _selections.erase(name);
}

Atomselections::~Atomselections() {
    for(Atomselections::iterator i = _selections.begin(); i != _selections.end(); ++i)
    {
        delete(i->second);
    }
}

// end of file

