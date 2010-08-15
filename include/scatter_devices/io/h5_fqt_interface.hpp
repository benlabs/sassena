/*
 *  h5_fqt_interface.hpp
 *
 *  Created on: May 26, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef SCATTER_DEVICES__IO__H5_FQT_INTERFACE_HPP_
#define SCATTER_DEVICES__IO__H5_FQT_INTERFACE_HPP_

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
#include <hdf5.h>


// other headers
#include "math/coor3d.hpp"


namespace H5FQTInterface {
    // return an array of indexes
    std::vector<size_t> init(const std::string filename,const std::vector<CartesianCoor3D>& qvectors,size_t nf);
    std::vector<CartesianCoor3D> get_qvectors(const std::string filename,const std::vector<size_t>& qindexes);
    void store(const std::string filename,const  size_t qindex, const std::vector<std::complex<double> >& fqt);
    
};

#endif

// end of file
