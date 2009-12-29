/*
 *  scatter_devices/scatter_devices.hpp
 *
 *  Created on: May 26, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef SCATTER_DEVICES__SCATTER_DEVICE_HPP_
#define SCATTER_DEVICES__SCATTER_DEVICE_HPP_

// common header
#include "common.hpp"

// standard header
#include <complex>
#include <map>
#include <string>
#include <sys/time.h>
#include <vector>

// special library headers
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/mpi.hpp>

// other headers
#include "coor3d.hpp"
#include "particle_trajectory.hpp"
#include "sample.hpp"
#include "scatter_factors.hpp"
#include "timer.hpp"

class ScatterDevice {
public: 
	Timer timer;	
	
	virtual void execute(CartesianCoor3D q) = 0;
	virtual std::vector<std::complex<double> >& get_spectrum() = 0;
	
};


#endif

//end of file
