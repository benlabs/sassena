/*
 *  analysis.hpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef ANALYSIS_HPP_
#define ANALYSIS_HPP_

// common header
#include "common.hpp"

// standard header
#include <complex>
#include <map>
#include <string>
#include <vector>

// special library headers

// other headers
#include "atomselection.hpp"
#include "coor3d.hpp"
#include "sample.hpp"

extern double b0_cached;

// a wrapper class for all the useful calculations
namespace Analysis {
	std::complex<double> background(Sample& sample,int r, double h, CartesianCoor3D q,std::map<std::string,std::vector<double> >& kappas);
	std::pair<double,double> background_avg(Sample& sample,int r, double h, CartesianCoor3D q);
	
	// use atomselection as argument to allow for selection flexibility 
	void set_scatteramp(Sample& sample,Atomselection as,CartesianCoor3D q,bool background);	

	// this is a wrapper which selects the right type of averaging!
	void scatter(Sample& sample,Atomselection as,CartesianCoor3D q,uint32_t qseed, std::vector<std::complex<double> >& scattering_amplitudes);

	std::complex<double> scatter_none(Sample& sample,Atomselection as,CartesianCoor3D& q);		
	
	void qvectors_unfold_sphere                (std::string avvectors, CartesianCoor3D q, uint32_t qseed,double resolution, std::vector<CartesianCoor3D>& qvectors);
	void qvectors_unfold_cylinder              (std::string avvectors, CartesianCoor3D q, uint32_t qseed,double resolution, std::vector<CartesianCoor3D>& qvectors);
	
	void scatter_vectors             (Sample& sample,Atomselection as,std::vector<CartesianCoor3D>& qvectors,std::vector<std::complex<double> >& scattering_amplitudes);
	void scatter_sphere_multipole    (Sample& sample,Atomselection as,CartesianCoor3D q,double resolution,std::vector<std::complex<double> >& scattering_amplitudes);
	void scatter_cylinder_multipole  (Sample& sample,Atomselection as,CartesianCoor3D q,double resolution,std::vector<std::complex<double> >& scattering_amplitudes);

};


#endif
