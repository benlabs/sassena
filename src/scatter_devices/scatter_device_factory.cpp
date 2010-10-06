/*
 *  scatter_device_factory.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

// direct header
#include "scatter_devices/scatter_device_factory.hpp"

// other headers
#include "control.hpp"
#include "log.hpp"
#include "decomposition/decompose.hpp"
#include "decomposition/decomposition_plan.hpp"
//#include "scatter_devices/all_multipole_sphere.hpp"
//#include "scatter_devices/all_multipole_cylinder.hpp"
//#include "scatter_devices/all_vectors.hpp"
#include "scatter_devices/all_vectorsthread.hpp"
//#include "scatter_devices/self_vectors.hpp"
#include "scatter_devices/self_vectorsthread.hpp"
using namespace std;

ScatterDevice* ScatterDeviceFactory::create(
		boost::mpi::communicator& all_comm,
		Sample& sample,
        boost::asio::ip::tcp::endpoint fileservice_endpoint,
        boost::asio::ip::tcp::endpoint monitorservice_endpoint,        
		std::vector<size_t>& qindexes,
		std::vector<CartesianCoor3D>& qvectors)
{
    
    ScatterDevice* p_ScatterDevice = NULL;

    size_t NN = all_comm.size();
    size_t NF = sample.coordinate_sets.size();
    string target = Params::Inst()->scattering.target;
    size_t NA = sample.atoms.selections[target].indexes.size();
    size_t NQ = qindexes.size();

    // check sample for data
    if (NF<1) {
        Err::Inst()->write("No frames available. Aborting");
        throw;
    }

    if (NA<1) {
        Err::Inst()->write("No atoms available. Aborting");
        throw;
    }

    if (NQ<1) {
        Err::Inst()->write("No qvectors left to compute. Aborting");
        throw;
    }

    std::vector<size_t> colors(all_comm.size());    
    size_t partitions;
        
    ////////////////////////////////////////////////////////////
    // Decomposition
    ////////////////////////////////////////////////////////////
    if (all_comm.rank()==0) {

        // decompose parallel space into independent partitions
        // each operating on a distinct set of qvectors.
        
    	Info::Inst()->write(string("Searching for decomposition plan: "));
	    Info::Inst()->write(string("nodes    = ")+ boost::lexical_cast<string>(all_comm.size()));
	    Info::Inst()->write(string("qvectors = ")+ boost::lexical_cast<string>(NQ));
	    Info::Inst()->write(string("frames   = ")+ boost::lexical_cast<string>(NF));
	    Info::Inst()->write(string("atoms   = ")+ boost::lexical_cast<string>(NA));

        if (Params::Inst()->scattering.type == "self") {
    	    Info::Inst()->write(string("Self interference scattering detected. Applying atom decomposition."));            
        } else if (Params::Inst()->scattering.type == "all") {
    	    Info::Inst()->write(string("All interference scattering detected. Applying frame decomposition."));            
        } else {
    	    Err::Inst()->write(string("Scattering Interference type not understood. Must be 'self' or 'all'."));            
            throw;
        }

        // for coherent scattering:
        size_t NAF = NA;
        if (Params::Inst()->scattering.type == "self") {
            NAF = NA;
        } else if (Params::Inst()->scattering.type == "all") {
            NAF = NF;
        }

        size_t ELBYTESIZE = 24;
        size_t NMAXBYTESIZE = Params::Inst()->limits.memory.data;
        if (Params::Inst()->scattering.type == "self") {
            ELBYTESIZE=24; // 3 times double: xyz
        } else if (Params::Inst()->scattering.type == "all") {
            ELBYTESIZE=NA*24; // 3 times double times number of atoms = frame
        }
        
		DecompositionPlan dplan(NN,NQ,NAF,ELBYTESIZE,NMAXBYTESIZE);
		colors = dplan.colors();        

        partitions = dplan.partitions();
    }

    broadcast(all_comm,partitions,0);

    //broadcast colors
	broadcast(all_comm,reinterpret_cast<size_t*>(&colors[0]),colors.size(),0);
    
    boost::mpi::communicator partition_comm = all_comm.split(colors[all_comm.rank()]);

    // exchange absolute values for q vectors
    vector<CartesianCoor3D> finalqvectors = qvectors;

    EvenDecompose qindex_decomposition(qindexes.size(),partitions);
	vector<pair<size_t,CartesianCoor3D> > thispartition_QIV;

	// determine the partition this node lives in:
	size_t mypartition = colors[all_comm.rank()];

	// don't include any "leftover" worlds
	if (mypartition<partitions) {
		vector<size_t> qii = qindex_decomposition.indexes_for(mypartition);
		for(size_t j=0;j<qii.size();j++) {
			thispartition_QIV.push_back(make_pair(qindexes[qii[j]],finalqvectors[qii[j]]));
		}
	}
	
	// compute assignments
    std::vector<size_t> assignment;
	if (Params::Inst()->scattering.type == "self") {
        assignment = EvenDecompose(NA,partition_comm.size()).indexes_for(partition_comm.rank());
    } else if (Params::Inst()->scattering.type == "all") {
        assignment = EvenDecompose(NF,partition_comm.size()).indexes_for(partition_comm.rank());
    }
	
	/////////////////////////////////////////////////////////////
    // Creating of scattering devices
    ////////////////////////////////////////////////////////////

    // all_comm for inter-partition communication, parition_comm for intra-partition communication

    if (Params::Inst()->scattering.type == "self") {
    	p_ScatterDevice = new SelfVectorsThreadScatterDevice(
    			all_comm,
    			partition_comm,
    			sample,
    			thispartition_QIV,
                assignment,
		        fileservice_endpoint,
		        monitorservice_endpoint);
    }
    else if (Params::Inst()->scattering.type == "all"){
    	if (Params::Inst()->scattering.average.orientation.type == "vectors") {
        	p_ScatterDevice = new AllVectorsThreadScatterDevice(
        			all_comm,
        			partition_comm,
        			sample,
        			thispartition_QIV,
	                assignment,
        			fileservice_endpoint,
        			monitorservice_endpoint);
//    	} else if (Params::Inst()->scattering.average.orientation.type == "multipole") {
//    		if (Params::Inst()->scattering.average.orientation.multipole.type == "sphere") {
//            	p_ScatterDevice = new AllMSScatterDevice(
//            			all_comm,
//            			partition_comm,
//            			sample,
//            			thispartition_QIV,
//            			fileservice_endpoint,
//            			monitorservice_endpoint);
//    		} else if (Params::Inst()->scattering.average.orientation.multipole.type == "cylinder") {
//            	p_ScatterDevice = new AllMCScatterDevice(
//            			all_comm,
//            			partition_comm,
//            			sample,
//            			thispartition_QIV,
//            			fileservice_endpoint,
//            			monitorservice_endpoint);
//    		} else {
//    			Err::Inst()->write(string("scattering.average.orientation.multipole.type not understood: ")+Params::Inst()->scattering.average.orientation.multipole.type);
//    			throw;
//    		}
    	} else if (Params::Inst()->scattering.average.orientation.type == "none") {
        	p_ScatterDevice = new AllVectorsThreadScatterDevice(
        			all_comm,
        			partition_comm,
        			sample,
        			thispartition_QIV,
        			assignment,
        			fileservice_endpoint,
        			monitorservice_endpoint);
    	}
    }    

    if (p_ScatterDevice==NULL) {
    	Err::Inst()->write("Error initializing ScatterDevice");
    	throw;
    }

    return p_ScatterDevice;
}
// end of file
