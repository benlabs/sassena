/*
 *  background_analyzer.hpp
 *
 *  Created on: May 26, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef BACKGROUND_ANALYZER_HPP_
#define BACKGROUND_ANALYZER_HPP_

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
#include "atomselection.hpp"
#include "coor3d.hpp"
#include "coordinate_sets.hpp"
#include "frame.hpp"
#include "frames.hpp"
#include "particle_trajectory.hpp"
#include "sample.hpp"
#include "scatter_factors.hpp"
#include "timer.hpp"

class BackgroundAnalyzer {
	
};

#endif 

// end of file