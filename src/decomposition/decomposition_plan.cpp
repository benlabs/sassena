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
#include "decomposition/decomposition_plan.hpp"

// standard header
#include <complex>
#include <fstream>

// special library headers
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

// other headers
#include "math/coor3d.hpp"
#include "decomposition/decompose.hpp"
#include "control.hpp"
#include "log.hpp"


using namespace std;

DecompositionPlan::DecompositionPlan(size_t nn,size_t nq,size_t nf) {
    
    m_nn = nn;
    m_nq = nq;
    m_nf = nf;
    
	// calculate some partioning schemes, init penalty w/ qsplit==1
	size_t bestworldsplit = 1;
	size_t bestcolwidth = nn;
	size_t penalty = compute_penalty(nq,nf,nn,bestworldsplit);
	
	if (Params::Inst()->limits.decomposition.partitions.automatic) {
        size_t maxworldsplit = ((nn<nq) ? nn : nq ); // natural limit: number of q vectors
	    if (Params::Inst()->limits.decomposition.partitions.max < maxworldsplit) maxworldsplit = Params::Inst()->limits.decomposition.partitions.max;
		for(size_t worldsplit = 2; worldsplit <= maxworldsplit; ++worldsplit)
		{
			size_t newpenalty = compute_penalty(nq,nf,nn,worldsplit);
			size_t colwidth = nn / worldsplit;
			size_t colheight = nf / colwidth + ( (nf % colwidth)==0 ? 0 : 1 );

			// only do a partition split if we can afford it memory-wise!
			if (worldsplit >= ( Params::Inst()->runtime.limits.cache.coordinate_sets / colheight ) ) break;

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
		bestworldsplit = Params::Inst()->limits.decomposition.partitions.count;
		bestcolwidth = nn / bestworldsplit;
		
		if (bestcolwidth<1) {
			Err::Inst()->write("Not allowed: Number of partitions higher than number of nodes");
			Err::Inst()->write("adjust: limits.decomposition.partitions.count");
			throw;
		}
		
		penalty = compute_penalty(nq,nf,nn,bestworldsplit);
	}
		
	if (nf<nn) bestcolwidth = nf;	
	
	// store this variables:
	m_bestcolwidth = bestcolwidth;
	m_penalty = penalty;
	m_bestworldsplit = bestworldsplit;
	
	// only allow successful construction if load imbalance is sufficiently small

	if (static_imbalance() > Params::Inst()->limits.decomposition.static_imbalance) {
			Err::Inst()->write(string("Total static load balance too high. Try a different number of nodes."));

			Info::Inst()->write(string("Static imbalance limit: ")+to_s(Params::Inst()->limits.decomposition.static_imbalance));
			Info::Inst()->write(string("Computed imbalance: ")+to_s(static_imbalance()));
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


		throw;
	}		
}

std::vector<size_t> DecompositionPlan::colors() {
		
	std::vector<size_t> colors;
	for(size_t i=0;i<m_nn;i++) {
		colors.push_back( i / m_bestcolwidth);
	}
	return colors;
}

double DecompositionPlan::static_imbalance() {
	return ( m_penalty/(1.0*m_nf*m_nq) );
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

		size_t leftoverworlds = 0;
		if (worldsplit>nn)  {
			leftoverworlds = worldsplit-nn;
			worldsplit = nn;
		}

		size_t colcount = worldsplit;
		size_t colwidth = nn / worldsplit;
		size_t colheight = nf / colwidth + ( (nf % colwidth)==0 ? 0 : 1 );
		size_t colpenalty = (colheight*colwidth)-nf;
		
		size_t worldpenalty = nn - colwidth * worldsplit;
	
		size_t rowcount = nq / worldsplit + ( (nq % worldsplit)==0 ? 0 : 1 );
		size_t rowpenalty = 0;
		if ((nq%colcount)!=0) {
			rowpenalty = (colcount - (nq % colcount) ) * (colwidth*colheight);
		} 
		
		// penalty due to insufficient usage of parallel worlds
		
		size_t penalty = nq * colpenalty  + rowcount * worldpenalty + rowpenalty + leftoverworlds * (colwidth*colheight) * rowcount;
		return penalty;
}

size_t DecompositionPlan::partitions() {
	return m_bestworldsplit;
}

size_t DecompositionPlan::partitionsize() {
	return m_nn / m_bestworldsplit;
}


// end of file
