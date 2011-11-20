/** \file
This file contains the class used to encapsulate the paritioning logic. It will find a reasonable partitioning scheme based on some job dependent conditions (number of nodes, number of frames and atoms..)

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/


// direct header
#include "decomposition/decomposition_plan.hpp"

// standard header
#include <complex>
#include <fstream>

// special library headers
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/lexical_cast.hpp>

// other headers
#include "math/coor3d.hpp"
#include "control.hpp"
#include "log.hpp"

using namespace std;

DecompositionParameters::DecompositionParameters(size_t NN,size_t NQ,size_t NAF,size_t NNpP,size_t elbytesize) {
    size_t NP = NN/NNpP;
    size_t NPused = NP;
    if (NQ<NP) NPused = NQ;
    size_t NNnotused = NN - NPused*NNpP;
    
    size_t NQcycles;
    if ((NQ%NP)==0) {
        NQcycles = NQ/NP;
    } else {
        NQcycles = NQ/NP + 1;
    }
    
    size_t NAFcycles;
    if ((NAF%NNpP)==0) {
        NAFcycles = NAF/NNpP;
    } else {
        NAFcycles =  NAF/NNpP + 1;
    }
    
    size_t penalty = NNnotused*NQcycles*NAFcycles;
    penalty +=  (NPused*NQcycles - NQ)*(NNpP*NAFcycles);
    penalty += (NNpP*NAFcycles - NAF)*NQ;
    
    // store for future reference
    m_penalty = penalty;
    m_NAFcycles = NAFcycles;
    m_NQcycles = NQcycles;
    m_NNpP = NNpP;

    m_NN = NN;
    m_NQ = NQ;
    m_NAF = NAF;
    m_NP = NP;

    m_elbytesize = elbytesize;
    m_nbytesize = NAFcycles*m_elbytesize;
}


DecompositionPlan::DecompositionPlan(size_t nn,size_t nq,size_t naf,size_t elbytesize,size_t nmaxbytesize) {

    // initialize pointer
    p_dp_best = NULL;
    
    if (naf<1) {
        Err::Inst()->write("No data to decompose.");
        throw;
    }
    
    size_t npmax = nn;
    if (npmax>nq) npmax = nq;
        
	if (Params::Inst()->limits.decomposition.partitions.automatic) {
	    Info::Inst()->write("Automatic decomposition. Searching for best utilization.");	
        size_t npmax = naf;
        if (naf>nn) npmax = nn;
        for (size_t nnpp=npmax;nnpp>=1;nnpp--) {
            DecompositionParameters* p_dp = new DecompositionParameters(nn,nq,naf,nnpp,elbytesize);
            if (p_dp->nbytesize()>nmaxbytesize) continue;
            if (p_dp_best == NULL) {
                p_dp_best = p_dp;
            } else {
                if (p_dp->penalty()<p_dp_best->penalty()) {
                    delete p_dp_best;
                    p_dp_best = p_dp;
                } else {
                    delete p_dp;
                }
            }
        } 
        
        if (p_dp_best == NULL) {
    		Err::Inst()->write("Automatic decomposition failed to match the necessary requirements.");	
    		Err::Inst()->write("Either change the partition size manually or change the number of nodes.");
    		Err::Inst()->write("Beware that size of a partition <= frames / atoms (depends)");
    		Err::Inst()->write("Limits:");
    		Err::Inst()->write(string("limits.stage.memory.data=")+boost::lexical_cast<string>(nmaxbytesize));
            Err::Inst()->write("Minimal Requirements:");         
            DecompositionParameters dp(nn,nq,naf,npmax,elbytesize);   
    		Err::Inst()->write(string("limits.stage.memory.data=")+boost::lexical_cast<string>(dp.nbytesize()));            	
    		throw;
        }
        
    } else {
        Warn::Inst()->write("Manual decomposition. This might not yield the best utilization!");	
        size_t nnpp = Params::Inst()->limits.decomposition.partitions.size;
        if (nnpp>nn) {
            nnpp=nn;
            Warn::Inst()->write("Partition size larger than NN. Setting NNPP=NN.");	
            Warn::Inst()->write(string("New partition size: ")+boost::lexical_cast<string>(nnpp));	
        }
        if (Params::Inst()->limits.decomposition.partitions.size>naf) {
            nnpp=naf;
            Warn::Inst()->write("Partition size larger than NAF. Setting NNPP=NAF.");	
            Warn::Inst()->write(string("New partition size: ")+boost::lexical_cast<string>(nnpp));	
        }
        
        p_dp_best = new DecompositionParameters(nn,nq,naf,nnpp,elbytesize); 
        
        if (p_dp_best == NULL) {
    		Err::Inst()->write("Manual decomposition failed.");	
    		throw;
        }
    }
        
    Info::Inst()->write("Final decomposition parameters:");
    Info::Inst()->write(string("NN                : ")+boost::lexical_cast<string>(p_dp_best->get_NN()));
    Info::Inst()->write(string("NQ                : ")+boost::lexical_cast<string>(p_dp_best->get_NQ()));
    Info::Inst()->write(string("NAF               : ")+boost::lexical_cast<string>(p_dp_best->get_NAF()));
    Info::Inst()->write(string("NP                : ")+boost::lexical_cast<string>(p_dp_best->get_NP()));
    Info::Inst()->write(string("NNpP              : ")+boost::lexical_cast<string>(p_dp_best->get_NNpP()));
    Info::Inst()->write(string("NAFcycles         : ")+boost::lexical_cast<string>(p_dp_best->get_NAFcycles()));
    Info::Inst()->write(string("NQcycles          : ")+boost::lexical_cast<string>(p_dp_best->get_NQcycles()));
    size_t used =  (p_dp_best->get_NQ()*p_dp_best->get_NAF());
    size_t wasted = p_dp_best->penalty();
    double utilization = (used)*1.0/(used+wasted);
    Info::Inst()->write(string("CompEl (WASTE/USE): ")+boost::lexical_cast<string>(wasted)+string("/")+boost::lexical_cast<string>(used));
    Info::Inst()->write(string("utilization(1=best)    : ")+boost::lexical_cast<string>(utilization));

    if (utilization < Params::Inst()->limits.decomposition.utilization) {
		Err::Inst()->write(string("Utilization too low. Aborting. Change the number of nodes or the threshold value (")+boost::lexical_cast<string>(Params::Inst()->limits.decomposition.utilization)+string(")"));
        delete p_dp_best;
        p_dp_best = NULL;
		throw;
    }
    
}

DecompositionPlan::~DecompositionPlan() {
    if (p_dp_best!=NULL) delete p_dp_best;
}

size_t DecompositionParameters::PID(size_t rank) {
    return (rank / m_NNpP);
}

size_t DecompositionParameters::penalty() {
    return (m_penalty);
}

double DecompositionParameters::penalty_percent() {
    return (m_penalty / (1.0*m_NAF*m_NQ));
}

std::vector<size_t> DecompositionPlan::colors() {
		
	std::vector<size_t> colors;
    size_t NN = p_dp_best->get_NN();
    for(size_t i=0;i<NN;i++) {
		colors.push_back( p_dp_best->PID(i) );
	}
	return colors;
}

double DecompositionPlan::static_imbalance() {
    return p_dp_best->penalty_percent();
}

size_t DecompositionPlan::partitions() {
	return p_dp_best->get_NP();
}

size_t DecompositionPlan::partitionsize() {
	return p_dp_best->get_NNpP();
}

size_t DecompositionPlan::nbytesize() {
    return p_dp_best->nbytesize();
}

// end of file
