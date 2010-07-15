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
        ar & qvectors;
        ar & amplitudes;
    }
	///////////////////
public:
	// 2d character: Nqvectors x 3
    std::vector<CartesianCoor3D> qvectors;
    
    // 3d character: Nqvectors x Nframes x 2
    std::vector< std::vector<std::complex<double> > > amplitudes;
    

	ScatterSpectrum() {}	
	ScatterSpectrum(std::vector<ScatterSpectrum> scatspecs);
	
    void add(CartesianCoor3D q, std::vector<std::complex<double> > spectrum);
    size_t size() { return qvectors.size(); }
    
    void clear();
    
    //std::pair<CartesianCoor3D,std::vector<std::complex<double> > >& operator[](size_t index) { return make_pair(qvectors[index],amplitudes[index]); }
	
	void write_average(std::string fname,std::string format);
	void write_plain(std::string fname,std::string format);

	void read_plain(std::string fname,std::string format);
		
    void transform();
};

#endif

// end of file
