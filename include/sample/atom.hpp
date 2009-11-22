/*
 *  atom.hpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
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


class Atom {
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & ID;
		ar & original_name;
		ar & residue_name;
		ar & chainid;
		ar & resseq;
		ar & occupancy;
		ar & tempFactor;
		ar & element;
		ar & charge;
		ar & name;   
		ar & index;  
		ar & mass;
		ar & excluded_volume;
		ar & scatteramp;
		ar & beta;
		ar & x;
		ar & y;
		ar & z;		
		ar & particle;
		ar & solvent;
		ar & kappa;
		ar & volume;
    }
	///////////////////
public:

	Atom() : particle(false), solvent(false) {}

	size_t ID;
	// additional info needed for "rewriting":
	std::string original_name;
	std::string residue_name;
	std::string chainid;
	std::string resseq;
	std::string occupancy;
	std::string tempFactor;
	std::string element;
	std::string charge;
			
	std::string name;	
	size_t index;
	double mass;
	double excluded_volume;
	double scatteramp;
	// x-ray scattering factor q-dependent.
	// ignore for now
	//	double scattering_factor_xray;
	double beta;
	double x,y,z;
	// indicators:
	bool particle;
	bool solvent;
	double kappa;
	double volume;
};

#endif
