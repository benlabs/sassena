/*
 *  scatterspectrum.hpp
 *
 *  Created on: May 26, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef SCATTERSPECTRUM_HPP_
#define SCATTERSPECTRUM_HPP_

// common header
#include "common.hpp"

// standard header
#include <complex>
#include <map>
#include <string>
#include <vector>

// special library headers
#include <boost/serialization/access.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>

// other headers
#include "math/coor3d.hpp"


class ScatterSpectrum : public std::vector<std::pair<CartesianCoor3D,std::vector<std::complex<double> > > > {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & boost::serialization::base_object<std::vector<std::pair<CartesianCoor3D,std::vector<std::complex<double> > > > >(*this);
    }
	///////////////////
public:
	ScatterSpectrum() {}	
	ScatterSpectrum(std::vector<ScatterSpectrum> scatspecs);
	
	void write_average(std::string fname,std::string format);
	void write_plain(std::string fname,std::string format);
	
};

#endif

// end of file
