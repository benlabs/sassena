/** \file 
This file contains a class which manages the generation of the proper scattering factors which may depend on atom and q vector lengths. 

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/


#ifndef SCATTER_DEVICES__SCATTER_FACTORS_HPP_
#define SCATTER_DEVICES__SCATTER_FACTORS_HPP_

// common header
#include "common.hpp"

// standard header
#include <string>
#include <vector>

// special library headers

// other headers
#include "sample.hpp"

//forward declaration...


/** 
Efficient management class for scattering factors. 
*/
class ScatterFactors {
	Sample* p_sample;
	IAtomselection* p_selection;
	
	bool m_background;
	
	std::vector<double> factors;
    std::vector<double> m_kappas;
public:
	ScatterFactors();
	
	// use these to initialize the scatterfactors set:
	void set_selection(IAtomselection* selection);
	void set_sample(Sample& sample);
	
	IAtomselection* get_selection();
	
	double get(size_t atomselectionindex);	
	std::vector<double>& get_all();
	
	void update(CartesianCoor3D q);
	void update_kappas();
    double compute_background(CartesianCoor3D q);
    
	void set_background(bool status);
};

#endif

// end of file
