/*
 *  scatter_factors.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */
// direct header
#include "scatter_factors.hpp"

// standard header
#include <fstream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>

// other headers
#include "sample/coordinate_sets.hpp"
#include "sample/frame.hpp"
#include "coor3d.hpp"
#include "control.hpp"

using namespace std;

ScatterFactors::ScatterFactors() {
	m_background=true;
}

void ScatterFactors::update(CartesianCoor3D q) {
	
	double background_sl = Params::Inst()->scattering.background.factor;

	for(size_t i = 0; i < p_selection->size(); ++i)
	{
		size_t atomID = p_sample->atoms[p_selection->at(i)].ID;
		double ql = q.length();
		double sf = Database::Inst()->sfactors.get(atomID,ql);
	
		// calculate effective scattering length:
		if (m_background) {
			double k  = p_sample->atoms[p_selection->at(i)].kappa;
			double v = Database::Inst()->volumes.get(atomID);
			double efactor = Database::Inst()->exclusionfactors.get(atomID,k*v,ql);

			sf = sf - background_sl*efactor;
		}

		factors[i] = sf;
	}
	
}

void ScatterFactors::set_selection(Atomselection& selection) {
	factors.resize(selection.size());
	p_selection = &selection;
}

Atomselection& ScatterFactors::get_selection() {
	return *p_selection;
}

void ScatterFactors::set_sample(Sample& sample) {
	p_sample = &sample;
}

double ScatterFactors::get(size_t atomselectionindex) {
	return factors[atomselectionindex];
}	

vector<double>& ScatterFactors::get_all() {
	return factors;
}

void ScatterFactors::set_background(bool status) {
	m_background = status;
}