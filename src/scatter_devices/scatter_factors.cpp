/** \file
This file contains a class which manages the generation of the proper scattering factors which may depend on atom and q vector lengths. 

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/
 
// direct header
#include "scatter_devices/scatter_factors.hpp"

// standard header
#include <fstream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>

// other headers
#include "sample/coordinate_sets.hpp"
#include "sample/frame.hpp"
#include "math/coor3d.hpp"
#include "control.hpp"
#include "log.hpp"

using namespace std;

ScatterFactors::ScatterFactors() {
	m_background=true;
	
	
}

void ScatterFactors::update_kappas() {
    // kappas fully map atoms
    // this way indexes can be used as addresses

    std::vector<ScatteringBackgroundKappaParameters>& kappas = Params::Inst()->scattering.background.kappas;
    
    if (p_sample==NULL) {
        Err::Inst()->write("ScatterFactors::update_kappas: sample not set");
        throw;
    }
    
    for(size_t i = 0; i < kappas.size(); ++i)
    {
        IAtomselection* selection = p_sample->atoms.selections[kappas[i].selection];
        for(size_t j = 0; j < selection->size(); ++j)
        {
            m_kappas[(*selection)[j]] = kappas[i].value;
        }
    }
    
}

void ScatterFactors::update(CartesianCoor3D q) {
	
	double background_sl = Params::Inst()->scattering.background.factor;

	for(size_t i = 0; i < p_selection->size(); ++i)
	{
		size_t atomID = p_sample->atoms[(*p_selection)[i]];
		double ql = q.length();
		double sf = Database::Inst()->sfactors.get(atomID,ql);
	
		// calculate effective scattering length:
		if (m_background) {
			double k  = m_kappas[(*p_selection)[i]];
			double v = Database::Inst()->volumes.get(atomID);
			double efactor = Database::Inst()->exclusionfactors.get(atomID,k*v,ql);

			sf = sf - background_sl*efactor;
		}

		factors[i] = sf;
	}
	
}

void ScatterFactors::set_selection(IAtomselection* selection) {
	factors.resize(selection->size());
	p_selection = selection;
}

IAtomselection* ScatterFactors::get_selection() {
	return p_selection;
}

void ScatterFactors::set_sample(Sample& sample) {
	p_sample = &sample;
	
    m_kappas.assign(p_sample->atoms.size(),1.0); // default to factor of 1.0
    update_kappas();
}

double ScatterFactors::get(size_t atomselectionindex) {
	return factors[atomselectionindex];
}	

vector<double>& ScatterFactors::get_all() {
	return factors;
}

double ScatterFactors::compute_background(CartesianCoor3D q) {

    double efactor_sum = 0;
    double sf_sum=0;
	double ql = q.length();

	for(size_t i = 0; i < p_selection->size(); ++i)
	{
		size_t atomID = p_sample->atoms[(*p_selection)[i]];
		double sf = Database::Inst()->sfactors.get(atomID,ql);

		double k  = m_kappas[(*p_selection)[i]];
		double v = Database::Inst()->volumes.get(atomID);
		double efactor = Database::Inst()->exclusionfactors.get(atomID,k*v,ql);

        efactor_sum+=efactor;
        sf_sum+=sf;            
	}
    return sf_sum/efactor_sum;
}


void ScatterFactors::set_background(bool status) {
	m_background = status;
}

// end of file
