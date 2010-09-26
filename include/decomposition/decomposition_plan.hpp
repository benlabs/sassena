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


 class DecompositionParameters {
     // initial

     size_t m_NN;
     size_t m_NQ;
     size_t m_NAF;
     size_t m_NP;    

     size_t m_penalty;
     size_t m_NAFcycles;
     size_t m_NQcycles;
     size_t m_NNpP;
     
     size_t m_elbytesize;
     size_t m_nbytesize;
  public:
     DecompositionParameters(size_t nn,size_t nq,size_t naf,size_t np,size_t elbytesize);

     size_t penalty();
     double penalty_percent();
     size_t PID(size_t rank); // maps a parition ID to a node rank

     size_t get_NN() { return m_NN; }     
     size_t get_NQ() { return m_NQ; }     
     size_t get_NAF() { return m_NAF; }     
     size_t get_NP() { return m_NP; }     

     size_t get_NQcycles() { return m_NQcycles; }     
     size_t get_NAFcycles() { return m_NAFcycles; }     

     size_t get_NNpP() { return m_NNpP; }     
     size_t nbytesize() { return m_nbytesize; }
 };
 
class DecompositionPlan {

    DecompositionParameters* p_dp_best;
	
public:
	DecompositionPlan(size_t nn,size_t nq,size_t naf,size_t elbytesize,size_t nmaxbytesize);
    ~DecompositionPlan();
    
	std::vector<size_t> colors();

	double static_imbalance();
	size_t penalty();
	size_t partitions();
	size_t partitionsize();

    size_t nbytesize();
};


#endif

// end of file
