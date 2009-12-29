/*
 *  scatterdevices.hpp
 *
 *  Created on: May 26, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef SCATTER_DEVICES__MULTIPOLE_SPHERE_HPP_
#define SCATTER_DEVICES__MULTIPOLE_SPHERE_HPP_

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
#include "sample.hpp"

#include "coor3d.hpp"
#include "particle_trajectory.hpp"
#include "scatter_factors.hpp"
#include "timer.hpp"

#include "scatter_devices/all.hpp"

class AllMSScatterDevice : public AllScatterDevice {
private:
    std::vector<std::pair<long,long> > moments;
    CartesianCoor3D q;
    
	// have to be implemented by concrete classes:
    void init(CartesianCoor3D& q);
    size_t get_numberofmoments();	
	void scatter(size_t moffset,size_t mcount);
    void norm();	

public:
	AllMSScatterDevice(boost::mpi::communicator& thisworld, Sample& sample);

};


#endif

//end of file