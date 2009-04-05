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
	double scatter(Sample& sample,Atomselection as,CartesianCoor3D q);

	std::complex<double> scatter_none(Sample& sample,Atomselection as,CartesianCoor3D q);		
		
	double scatter_sphere_bf_rasterlinear      (Sample& sample,Atomselection as,CartesianCoor3D q,double resolution);	
	double scatter_sphere_bf_rastersurface     (Sample& sample,Atomselection as,CartesianCoor3D q,double resolution);	
	double scatter_sphere_bf_triangulation     (Sample& sample,Atomselection as,CartesianCoor3D q,double resolution);	
	double scatter_sphere_bf_mc_phiacos        (Sample& sample,Atomselection as,CartesianCoor3D q,double resolution);	
	double scatter_sphere_bf_mc_doublesqrt     (Sample& sample,Atomselection as,CartesianCoor3D q,double resolution);	
	double scatter_sphere_bf_mc_quaternion     (Sample& sample,Atomselection as,CartesianCoor3D q,double resolution);	
	double scatter_sphere_bf_mc_boostunisphere (Sample& sample,Atomselection as,CartesianCoor3D q,double resolution);	
	double scatter_sphere_multipole            (Sample& sample,Atomselection as,CartesianCoor3D q,double resolution);
	double scatter_cylinder_bf                 (Sample& sample,Atomselection as,CartesianCoor3D q,double resolution);		
	double scatter_cylinder_multipole          (Sample& sample,Atomselection as,CartesianCoor3D q,double resolution);
	
};


#endif
