/** \file
This file contains an refined version of the abstract scatter device, used for performing vector based orientationally averaged scattering calculations.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
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
#include "control.hpp"
#include "log.hpp"
#include "sample.hpp"
#include "stager/data_stager.hpp"

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

bool AbstractVectorsScatterDevice::ram_check() {
	// inherit ram requirements for parent class
	bool state = AbstractScatterDevice::ram_check();
	return state;
}

void AbstractVectorsScatterDevice::print_pre_stage_info() {
    if (allcomm_.rank()==0) {
        Info::Inst()->write("Staging data...");
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
    		    cout << qvector << endl;                                
            }
    	}	    
    }    
    }
    
    if (allcomm_.rank()==0) {
        Info::Inst()->write("Starting computation...");
        if (Params::Inst()->scattering.average.orientation.vectors.size()>0) {
            Info::Inst()->write(string("Orientational averaging includes ")+boost::lexical_cast<string>(Params::Inst()->scattering.average.orientation.vectors.size())+string(" vectors"));            
        }
    }
}

void AbstractVectorsScatterDevice::print_post_runner_info() { }

double AbstractVectorsScatterDevice::progress() {
    double scale1 = 1.0/vectors_.size();
    // scale2 is not preserved for cases where the size is dependent on the vector (e.g. 0 does not require any "orientations")
    //double scale2 = 1.0/subvector_index_.size();
    double scale2 = 1.0/NM;
    
    double base1 =  current_vector_*scale1;
    double base2 =  current_subvector_*scale1*scale2;
    return base1 + base2;
}



void AbstractVectorsScatterDevice::init_subvectors(CartesianCoor3D& q) {
    subvector_index_.clear();	
	
	if (Params::Inst()->scattering.average.orientation.vectors.size()>0) {
		if (Params::Inst()->scattering.average.orientation.vectors.type=="file") {
			double ql = q.length();
			
			for(size_t i = 0; i < Params::Inst()->scattering.average.orientation.vectors.size(); ++i)
			{
				subvector_index_.push_back(ql*Params::Inst()->scattering.average.orientation.vectors[i]);
			}					
		} else if (Params::Inst()->scattering.average.orientation.vectors.type=="sphere") {
		   double ql = q.length();
           
		   for(size_t i = 0; i < Params::Inst()->scattering.average.orientation.vectors.size(); ++i)
		   {
		   	subvector_index_.push_back(ql*Params::Inst()->scattering.average.orientation.vectors[i]);
		   }					
		} else if (Params::Inst()->scattering.average.orientation.vectors.type=="cylinder") {			
			CartesianCoor3D o = Params::Inst()->scattering.average.orientation.axis;
			
			// constructs a base out of thin air
            CartesianVectorBase base(o);
            CartesianCoor3D qprojected = base.project(q);
            CylinderCoor3D qcylinder(qprojected);
            
            // components are r,phi,z (order of base vectors for cylinder)
            if (qcylinder.r==0) {
                subvector_index_.push_back( qprojected );
            } else {
				for(size_t i = 0; i < Params::Inst()->scattering.average.orientation.vectors.size(); ++i)
				{
					CartesianCoor3D& vec = Params::Inst()->scattering.average.orientation.vectors[i];
					// it's "easier" to just construct the vectors from the base, then to do coordinate transformations
					CartesianCoor3D qnew = qcylinder.z*base[2] + qcylinder.r*(vec.x * base[0] + vec.y * base[1]);
					subvector_index_.push_back(qnew);
				}
            }    
			
			if (Params::Inst()->debug.print.orientations) {     
            if (subvector_index_.size()>0) {
                if (allcomm_.rank()==0) {		
            		Info::Inst()->write("Final Qvectors orientations used for averaging: ");
                    for (size_t i=0;i<subvector_index_.size();i++) {
                        string qvector = "";
                        qvector += boost::lexical_cast<string>(subvector_index_[i].x);
                        qvector += " ";
                        qvector += boost::lexical_cast<string>(subvector_index_[i].y);
                        qvector += " ";
                        qvector += boost::lexical_cast<string>(subvector_index_[i].z);
            		    cout << qvector << endl;                                
                    }
            	}	    
            }    
            }
            			
		}
	} else {
		subvector_index_.push_back(q);
	}
	
	// set shortcut variable
	NM = subvector_index_.size();
}

// end of file
