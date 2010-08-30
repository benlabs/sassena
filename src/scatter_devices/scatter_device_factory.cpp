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
#include "scatter_devices/all_multipole_sphere.hpp"
#include "scatter_devices/all_multipole_cylinder.hpp"
#include "scatter_devices/all_vectors.hpp"
#include "scatter_devices/all_vectorsthread.hpp"
#include "scatter_devices/self_vectors.hpp"
#include "scatter_devices/io/h5_fqt_interface.hpp"

using namespace std;

ScatterDevice* ScatterDeviceFactory::create(
		boost::mpi::communicator& scatter_comm,
		Sample& sample,
		std::vector<CartesianCoor3D>& qvectors,
		std::string fqt_filename)
{
    
    ScatterDevice* p_ScatterDevice = NULL;

    // initialize associated data file and use checkpoints
    vector<size_t> qindexes;
    if (scatter_comm.rank()==0) {
        qindexes = H5FQTInterface::init(fqt_filename,qvectors,sample.coordinate_sets.size());

        size_t nq = qindexes.size();
        broadcast(scatter_comm,nq,0);

        broadcast(scatter_comm,reinterpret_cast<size_t*>(&qindexes[0]),qindexes.size(),0);

        scatter_comm.barrier();
    } else {
        size_t nq=0;
        broadcast(scatter_comm,nq,0);
        qindexes.resize(nq);
        broadcast(scatter_comm,reinterpret_cast<size_t*>(&qindexes[0]),qindexes.size(),0);
        scatter_comm.barrier();
    }
    // qindexes contain absolute index positions which have to be used when writing to the data file
    

    // let the first node read the qvector values from file
    std::vector<size_t> colors(scatter_comm.size());
    vector<CartesianCoor3D> finalqvectors;
    size_t partitions;
    
    if (scatter_comm.rank()==0) {

        // decompose parallel space into independent partitions
        // each operating on a distinct set of qvectors.
        
    	Info::Inst()->write(string("Searching for decomposition plan: "));
	    Info::Inst()->write(string("nodes    = ")+ to_s(scatter_comm.size()));
	    Info::Inst()->write(string("qvectors = ")+ to_s(qindexes.size()));
	    Info::Inst()->write(string("frames   = ")+ to_s(sample.coordinate_sets.size()));

		DecompositionPlan dplan(scatter_comm.size(),(qindexes.size()>0 ? qindexes.size():1),sample.coordinate_sets.size());
		colors = dplan.colors();
        
		Info::Inst()->write(string("Decomposition has ")+to_s(dplan.partitions())+string(" partitions"));
		Info::Inst()->write(string("Static imbalance factor (0 is best): ")+ to_s(dplan.static_imbalance()));

        finalqvectors = H5FQTInterface::get_qvectors(fqt_filename,qindexes);
        partitions = dplan.partitions();
    }

    broadcast(scatter_comm,partitions,0);

    //broadcast colors
	broadcast(scatter_comm,reinterpret_cast<size_t*>(&colors[0]),colors.size(),0);
    
    boost::mpi::communicator fqt_comm = scatter_comm.split(colors[scatter_comm.rank()]);

    size_t nfinalqvectors= finalqvectors.size();
	//broadcast qindexes,qvectors
    broadcast(scatter_comm,nfinalqvectors,0);

    if (scatter_comm.rank()!=0) {
    	finalqvectors.resize(nfinalqvectors);
    	qindexes.resize(nfinalqvectors);
    }
	broadcast(scatter_comm,&(finalqvectors[0].x),nfinalqvectors*3,0);
	broadcast(scatter_comm,&(qindexes[0]),nfinalqvectors,0);
    EvenDecompose qindex_decomposition(qindexes.size(),partitions);
	vector<pair<size_t,CartesianCoor3D> > thispartition_QIV;

	// determine the partition this node lives in:
	size_t mypartition = colors[scatter_comm.rank()];

	// don't include any "leftover" worlds
	if (mypartition<partitions) {
		vector<size_t> qii = qindex_decomposition.indexes_for(mypartition);
		for(size_t j=0;j<qii.size();j++) {
			thispartition_QIV.push_back(make_pair(qindexes[qii[j]],finalqvectors[qii[j]]));
		}
	}

    // create scatter device
    // scatter_comm for inter-partition communication, fqt_comm for intra-partition communication

    if (Params::Inst()->scattering.interference.type == "self") {
    	p_ScatterDevice = new SelfVectorsScatterDevice(
    			scatter_comm,
    			fqt_comm,
    			sample,
    			thispartition_QIV,
    			fqt_filename);
    }
    else if (Params::Inst()->scattering.interference.type == "all"){
    	if (Params::Inst()->scattering.average.orientation.type == "vectors") {
        	p_ScatterDevice = new AllVectorsScatterDevice(
        			scatter_comm,
        			fqt_comm,
        			sample,
        			thispartition_QIV,
        			fqt_filename);
    	} else if (Params::Inst()->scattering.average.orientation.type == "vectorsthread") {
        	p_ScatterDevice = new AllVectorsThreadScatterDevice(
        			scatter_comm,
        			fqt_comm,
        			sample,
        			thispartition_QIV,
        			fqt_filename);
    	} else if (Params::Inst()->scattering.average.orientation.type == "multipole") {
    		if (Params::Inst()->scattering.average.orientation.multipole.type == "sphere") {
            	p_ScatterDevice = new AllMSScatterDevice(
            			scatter_comm,
            			fqt_comm,
            			sample,
            			thispartition_QIV,
            			fqt_filename);
    		} else if (Params::Inst()->scattering.average.orientation.multipole.type == "cylinder") {
            	p_ScatterDevice = new AllMCScatterDevice(
            			scatter_comm,
            			fqt_comm,
            			sample,
            			thispartition_QIV,
            			fqt_filename);
    		} else {
    			Err::Inst()->write(string("scattering.average.orientation.multipole.type not understood: ")+Params::Inst()->scattering.average.orientation.multipole.type);
    			throw;
    		}
    	} else if (Params::Inst()->scattering.average.orientation.type == "none") {
        	p_ScatterDevice = new AllVectorsScatterDevice(
        			scatter_comm,
        			fqt_comm,
        			sample,
        			thispartition_QIV,
        			fqt_filename);
    	}
    }

    if (p_ScatterDevice==NULL) {
    	Err::Inst()->write("Error initializing ScatterDevice");
    	throw;
    }

    return p_ScatterDevice;
}
// end of file
