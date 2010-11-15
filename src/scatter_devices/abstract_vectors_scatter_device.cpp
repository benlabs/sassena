/*
 *  scatterdevices.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

// direct header
#include "scatter_devices/abstract_vectors_scatter_device.hpp"

// standard header
#include <algorithm>
#include <complex>
#include <fstream>

// special library headers
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/lexical_cast.hpp>

// other headers
#include "math/coor3d.hpp"
#include "decomposition/decompose.hpp"
#include "control.hpp"
#include "log.hpp"
#include "sample.hpp"
#include "scatter_devices/data_stager.hpp"

using namespace std;

AbstractVectorsScatterDevice::AbstractVectorsScatterDevice(
    boost::mpi::communicator allcomm,
    boost::mpi::communicator partitioncomm,
    Sample& sample,
    std::vector<CartesianCoor3D> vectors,
    size_t NAF,
    boost::asio::ip::tcp::endpoint fileservice_endpoint,
	boost::asio::ip::tcp::endpoint monitorservice_endpoint
) :
    AbstractScatterDevice(
        allcomm,
        partitioncomm,
        sample,
        vectors,
        NAF,
        fileservice_endpoint,
        monitorservice_endpoint
    ),
    current_subvector_(0)
{
	NM = Params::Inst()->scattering.average.orientation.vectors.size();
    if (NM==0) NM=1;
	sample_.coordinate_sets.set_representation(CARTESIAN);		
}


void AbstractVectorsScatterDevice::print_pre_stage_info() {
    if (allcomm_.rank()==0) {
        Info::Inst()->write("Staging data...");
        Info::Inst()->write(("limits.stage.nodes=")+boost::lexical_cast<string>(Params::Inst()->limits.stage.nodes));
    }
}

void AbstractVectorsScatterDevice::print_post_stage_info() { }

void AbstractVectorsScatterDevice::print_pre_runner_info() {
    
    if (Params::Inst()->debug.print.orientations) {     
    if (Params::Inst()->scattering.average.orientation.vectors.size()>0) {
        if (allcomm_.rank()==0) {		
    		Info::Inst()->write("Qvectors orientations used for averaging: ");
            for (size_t i=0;i<Params::Inst()->scattering.average.orientation.vectors.size();i++) {
                string qvector = "";
                qvector += boost::lexical_cast<string>(Params::Inst()->scattering.average.orientation.vectors[i].x);
                qvector += " ";
                qvector += boost::lexical_cast<string>(Params::Inst()->scattering.average.orientation.vectors[i].y);
                qvector += " ";
                qvector += boost::lexical_cast<string>(Params::Inst()->scattering.average.orientation.vectors[i].z);
    		    Info::Inst()->write(qvector);                                
            }
    	}	    
    }    
    }
    
    if (allcomm_.rank()==0) {
        Info::Inst()->write("Starting computation...");
    }
}

void AbstractVectorsScatterDevice::print_post_runner_info() { }

double AbstractVectorsScatterDevice::progress() {
    double scale1 = 1.0/vectors_.size();
    double scale2 = 1.0/subvector_index_.size();
    
    double base1 =  current_vector_*scale1;
    double base2 =  current_subvector_*scale1*scale2;
    return base1 + base2;
}



void AbstractVectorsScatterDevice::init_subvectors(CartesianCoor3D& q) {
    subvector_index_.clear();	
	
	if (Params::Inst()->scattering.average.orientation.vectors.size()>0) {
		if (Params::Inst()->scattering.average.orientation.vectors.type=="sphere") {
			double ql = q.length();
			
			for(size_t i = 0; i < Params::Inst()->scattering.average.orientation.vectors.size(); ++i)
			{
				subvector_index_.push_back(ql*Params::Inst()->scattering.average.orientation.vectors[i]);
			}					
		} else if (Params::Inst()->scattering.average.orientation.vectors.type=="cylinder") {
			CartesianCoor3D o = Params::Inst()->scattering.average.orientation.vectors.axis;
			// make sure o is normalized;
			o = o / o.length();
			
			// get the part of the scattering vector perpenticular to the o- orientation
			CartesianCoor3D qparallel = (o*q)*o; 
			CartesianCoor3D qperpenticular = q - qparallel; 			
			double qperpenticular_l = qperpenticular.length();
			//double qparallel_l = qparallel.length();

			CartesianCoor3D e1 = o.cross_product(qperpenticular) ;
			CartesianCoor3D e2 = qperpenticular;
			
			if (qperpenticular_l==0.0)  { 
				subvector_index_.push_back( qparallel );
			}
			else {
				for(size_t i = 0; i < Params::Inst()->scattering.average.orientation.vectors.size(); ++i)
				{
					CartesianCoor3D& vec = Params::Inst()->scattering.average.orientation.vectors[i];
					CartesianCoor3D qnew = qparallel + vec.x * e1 + vec.y * e2;					
					subvector_index_.push_back(qnew);
				}
			}				
		}
	} else {
		subvector_index_.push_back(q);
	}
}
