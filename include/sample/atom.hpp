/** \file
This file contains a class which defines an atom.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/
 
#ifndef SAMPLE__ATOM_HPP_
#define SAMPLE__ATOM_HPP_

// common header
#include "common.hpp"

// standard header
#include <string>

// special library headers
#include <boost/serialization/access.hpp>
#include <boost/serialization/string.hpp>

// other headers

/** 
Type class which represents an atom.
*/
class Atom {
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & ID;
		ar & name;   

//		ar & original_name;
//		ar & residue_name;
//		ar & chainid;
//		ar & resseq;
//		ar & occupancy;
//		ar & tempFactor;
//		ar & element;
//		ar & charge;

		ar & index;  
		ar & mass;

//		ar & beta;
//		ar & x;
//		ar & y;
//		ar & z;		
    }
	///////////////////
public:

	size_t ID;
	std::string name;	

	// additional info needed for "rewriting":
//	std::string original_name;
//	std::string residue_name;
//	std::string chainid;
//	std::string resseq;
//	std::string occupancy;
//	std::string tempFactor;
//	std::string element;
//	std::string charge;
//			
	size_t index;
	double mass;
	// x-ray scattering factor q-dependent.
	// ignore for now
	//	double scattering_factor_xray;
//	double beta;
//	double x,y,z;
//	// indicators:
};

#endif

// end of file
