/*
 *  DecompositionPlan.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

// direct header
#include "decomposition_plan.hpp"

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
#include "log.hpp"
#include "parameters.hpp"
#include "database.hpp"
#include "sample.hpp"
#include "particle_trajectory.hpp"


using namespace std;

DecompositionPlan::DecompositionPlan(boost::mpi::communicator thisworld,vector<CartesianCoor3D>& qvectors,vector<size_t>& frames) {
	
	size_t nn = thisworld.size();
	size_t nq = qvectors.size();
	size_t nf = frames.size();
	
	// calculate some partioning schemes, init penalty w/ qsplit==1
	size_t bestworldsplit = 1;
	size_t bestcolwidth = nn;
	size_t penalty = compute_penalty(nq,nf,nn,bestworldsplit);

	if (Params::Inst()->debug.decomposition_split==0) {
		for(size_t worldsplit = 2; worldsplit <= ((nn<nq) ? nn : nq ); ++worldsplit)
		{
			size_t newpenalty = compute_penalty(nq,nf,nn,worldsplit);
			size_t colwidth = nn / worldsplit;
			size_t colheight = nf / colwidth + ( (nf % colwidth)==0 ? 0 : 1 );

			// hard limit: frames per node
			if (colheight>1000) break;

			// if equal, prefer higher split
			if (newpenalty<=penalty) {
				penalty = newpenalty;
				bestworldsplit = worldsplit;
				bestcolwidth = colwidth;
			}	
		}		
	} else {
		// calculate some partioning schemes, init penalty w/ qsplit==1
		bestworldsplit = Params::Inst()->debug.decomposition_split;
		bestcolwidth = nn / bestworldsplit;
		penalty = compute_penalty(nq,nf,nn,bestworldsplit);
	}

		
	if (nf<nn) bestcolwidth = nf;	
	
	
	// only allow successful construction if load imbalance is sufficiently small

	double fpenalty = penalty * 1.0/ (nf*nq);
	if (fpenalty > Params::Inst()->limits.static_load_imbalance_max) {
		if (thisworld.rank()==0) {
			Err::Inst()->write(string("Total static load balance too high. Try a different number of nodes."));

			Info::Inst()->write(string("Static imbalance limit: ")+to_s(Params::Inst()->limits.static_load_imbalance_max));
			Info::Inst()->write(string("Computed imbalance: ")+to_s(fpenalty));
			Info::Inst()->write(string("Decomposition Plan: split=")+to_s(bestworldsplit));
			
			Info::Inst()->write("For your layout, the following number of nodes have best partitioning:");

			vector<pair<double,size_t> > lis = scan_imbalance_spectrum(nq,nf,nq*nf);
			stringstream lss;
			size_t liscounter = 0;
			for(size_t i = 0; i < lis.size(); ++i)
			{
				if (lis[i].first!=0.0) liscounter++;
				lss << lis[i].second;
				if (lis[i].first==0.0) lss << "*";
				lss << ", ";
				if (liscounter>10) break;
			}
			Info::Inst()->write(string("Best numbers: ")+lss.str());
			
		}

		throw;
	}	
	
	// store this variables:
	m_bestcolwidth = bestcolwidth;
	m_penalty = penalty;
	m_bestworldsplit = bestworldsplit;
	
	m_qvectors = qvectors;
	m_frames = frames;
	
	
}

boost::mpi::communicator DecompositionPlan::split() {
		
	size_t mycolor = m_thisworldcomm.rank() / m_bestcolwidth;
		
	boost::mpi::communicator local = m_thisworldcomm.split(mycolor);
	return local;
}

vector<CartesianCoor3D> DecompositionPlan::qvectors() {

	vector<CartesianCoor3D> result;
	
	EvenDecompose e(m_qvectors.size(),m_bestworldsplit);
	size_t mycolor = m_thisworldcomm.rank() / m_bestcolwidth;
	
	if (mycolor<e.size()) {
		vector<size_t> qindexes = e.indexes_for(mycolor);
		for(size_t i = 0; i < qindexes.size(); ++i)
		{
			result.push_back( m_qvectors[ qindexes[i] ] );
		}			
	}		
	
	return result;
}

vector<size_t> DecompositionPlan::frames() {
	throw; // not implemented yet.
}

double DecompositionPlan::static_imbalance() {
	return ( 1.0 - m_penalty/(1.0*m_frames.size()*m_qvectors.size()) ); 
}

size_t DecompositionPlan::penalty() {
	return m_penalty;
}

vector<pair<double,size_t> > DecompositionPlan::scan_imbalance_spectrum(size_t nq, size_t nf, size_t maxnn = 10000) {
	vector<pair<double,size_t> > result;
	for(size_t i = 1; i <= maxnn; ++i)
	{
		// computate_penalty requires split being smaller than NN
		size_t maxj = (i<nq ? i : nq);
		for(size_t j = 1; j <= maxj; ++j)
		{
			result.push_back(make_pair(compute_penalty(nq,nf,i,j),i));			
		}
	}
	sort(result.begin(),result.end());
	return result;
}

size_t DecompositionPlan::compute_penalty(size_t nq, size_t nf, size_t nn,size_t worldsplit) {

	if (worldsplit>nn) {
		Err::Inst()->write("worldsplit > nn. That doesn't make sense. This seems to be a coding problem...");
		throw;
	}

	size_t colcount = worldsplit;
	size_t colwidth = nn / worldsplit;
	size_t colheight = nf / colwidth + ( (nf % colwidth)==0 ? 0 : 1 );
	size_t colpenalty = (colheight*colwidth)-nf;
		
	size_t worldpenalty = nn - colwidth * worldsplit;
	
	size_t rowcount = nq / worldsplit + ( (nq % worldsplit)==0 ? 0 : 1 );
	size_t rowpenalty = (nq % colcount) * (colwidth*colheight); // penalty due to insufficient usage of parallel worlds
		
	size_t penalty = nq * colpenalty  + rowcount * worldpenalty + rowpenalty;
	return penalty;
}

size_t DecompositionPlan::worlds() { 
	return m_bestworldsplit;
}
