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

#ifndef SCATTER_DEVICES__SCATTERSPECTRUM_HPP_
#define SCATTER_DEVICES__SCATTERSPECTRUM_HPP_

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


class ScatterSpectrum {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        ar & data;
    }
	///////////////////
	
    std::vector<std::pair<CartesianCoor3D,std::vector<std::complex<double> > > > data;
    
public:
	ScatterSpectrum() {}	
	ScatterSpectrum(std::vector<ScatterSpectrum> scatspecs);
	
    void add(CartesianCoor3D q, std::vector<std::complex<double> > spectrum);
    size_t size() { return data.size(); }
    
    std::pair<CartesianCoor3D,std::vector<std::complex<double> > >& operator[](size_t index) { return data[index]; }
	
	void write_average(std::string fname,std::string format);
	void write_plain(std::string fname,std::string format);

	void read_plain(std::string fname,std::string format);
		
    void transform();
};

#endif

// end of file
