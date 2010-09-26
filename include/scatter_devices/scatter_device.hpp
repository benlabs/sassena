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
#include "math/coor3d.hpp"
#include "report/timer.hpp"

class ScatterDevice {
public: 
	Timer timer;	
	
    virtual void compute() = 0;
	virtual void next() = 0;
	virtual void write()  = 0;
	
    virtual size_t status() = 0;
    
    virtual double progress() = 0;

};





#endif

//end of file
