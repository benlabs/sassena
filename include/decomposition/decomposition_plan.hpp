/*
 *  decompositionstrategy2d.hpp
 *
 *  Created on: May 26, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef DECOMPOSITION__DECOMPOSITIONPLAN_HPP_
#define DECOMPOSITION__DECOMPOSITIONPLAN_HPP_

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

class DecompositionPlan {

	size_t m_bestcolwidth;
	size_t m_penalty;
	size_t m_bestworldsplit;
	
	size_t m_nn;
	size_t m_nq;
	size_t m_nf;

	size_t compute_penalty(size_t nq, size_t nf, size_t nn,size_t worldsplit);	
	
	std::vector<std::pair<double,size_t> > scan_imbalance_spectrum(size_t nq, size_t nf, size_t maxnn);
	
public:
	DecompositionPlan(size_t nn,size_t nq,size_t nf);
	
	std::vector<size_t> colors();

	double static_imbalance();
	size_t penalty();
	size_t partitions();
	size_t partitionsize();

};


#endif

// end of file
