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

	boost::mpi::communicator m_thisworldcomm;
	
	size_t m_bestcolwidth;
	size_t m_penalty;
	size_t m_bestworldsplit;
	
	std::vector<size_t> m_qindexes;
	std::vector<size_t> m_frames;

	size_t compute_penalty(size_t nq, size_t nf, size_t nn,size_t worldsplit);	
	
	std::vector<std::pair<double,size_t> > scan_imbalance_spectrum(size_t nq, size_t nf, size_t maxnn);
	
public:
	DecompositionPlan(boost::mpi::communicator thisworld,std::vector<size_t>& qindexes,std::vector<size_t>& frames);
	
	boost::mpi::communicator split();
	std::vector<size_t> frames();
	std::vector<size_t> qindexes();
	
	double static_imbalance();
	size_t penalty();
	size_t worlds();
};


#endif

// end of file