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

#ifndef SCATTER_DEVICES__ALL_VECTORS_HPP_
#define SCATTER_DEVICES__ALL_VECTORS_HPP_

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

#include "math/coor3d.hpp"
#include "report/timer.hpp"

#include "scatter_devices/all.hpp"

class AllVectorsScatterDevice : public AllScatterDevice {
private:
    std::vector<CartesianCoor3D> qvectors;
    
	// have to be implemented by concrete classes:
    void init(CartesianCoor3D& q);
    size_t get_numberofmoments();	
	void scatter(size_t moffset,size_t mcount);
    void norm();	

public:
	AllVectorsScatterDevice(
			boost::mpi::communicator scatter_comm,
			boost::mpi::communicator fqt_comm,
			Sample& sample,
			std::vector<std::pair<size_t,CartesianCoor3D> > QVI,
			std::string fqt_filename
	);

};


#endif

//end of file
