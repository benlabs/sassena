/*
 *  background_analyzer.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

// direct header
#include "background_analyzer.hpp"

// standard header
#include <complex>
#include <fstream>

// special library headers
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

// other headers
#include "analysis.hpp"
#include "coor3d.hpp"
#include "decompose.hpp"
#include "control.hpp"
#include "sample/sample.hpp"
#include "smath.hpp"
#include "particle_trajectory.hpp"


using namespace std;

BackgroundAnalyzer::BackgroundAnalyzer(Sample& sample) {
	
}

void BackgroundAnalyzer::execute(Sample& sample) {
	
}