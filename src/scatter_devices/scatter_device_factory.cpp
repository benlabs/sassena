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
		boost::mpi::communicator& all_comm,
		Sample& sample,
		std::vector<CartesianCoor3D>& qvectors)
{
    
    ScatterDevice* p_ScatterDevice = NULL;

    size_t NN = all_comm.size();
    size_t NF = sample.coordinate_sets.size();
    string target = Params::Inst()->scattering.target;
    size_t NA = sample.atoms.selections[target].indexes.size();

    // check sample for data
    if (NF<1) {
        Err::Inst()->write("No frames available. Aborting");
        throw;
    }

    if (NA<1) {
        Err::Inst()->write("No atoms available. Aborting");
        throw;
    }

    // initialize associated data file and use checkpoints
    vector<size_t> qindexes;
    if (all_comm.rank()==0) {
        qindexes = H5FQTInterface::init(Params::Inst()->scattering.data.file,qvectors,sample.coordinate_sets.size());

        size_t nq = qindexes.size();
        broadcast(all_comm,nq,0);

        broadcast(all_comm,reinterpret_cast<size_t*>(&qindexes[0]),qindexes.size(),0);

        all_comm.barrier();
    } else {
        size_t nq=0;
        broadcast(all_comm,nq,0);
        qindexes.resize(nq);
        broadcast(all_comm,reinterpret_cast<size_t*>(&qindexes[0]),qindexes.size(),0);
        all_comm.barrier();
    }
    // qindexes contain absolute index positions which have to be used when writing to the data file
    size_t NQ = qindexes.size();

    if (NQ<1) {
        Err::Inst()->write("No qvectors left to compute. Aborting");
        throw;
    }

    // let the first node read the qvector values from file
    std::vector<size_t> colors(all_comm.size());
    vector<CartesianCoor3D> finalqvectors;
    size_t partitions;
        
    ////////////////////////////////////////////////////////////
    // Decomposition
    ////////////////////////////////////////////////////////////
    if (all_comm.rank()==0) {

        // decompose parallel space into independent partitions
        // each operating on a distinct set of qvectors.
        
    	Info::Inst()->write(string("Searching for decomposition plan: "));
	    Info::Inst()->write(string("nodes    = ")+ to_s(all_comm.size()));
	    Info::Inst()->write(string("qvectors = ")+ to_s(NQ));
	    Info::Inst()->write(string("frames   = ")+ to_s(NF));
	    Info::Inst()->write(string("atoms   = ")+ to_s(NA));

        if (Params::Inst()->scattering.interference.type == "self") {
    	    Info::Inst()->write(string("Self interference scattering detected. Applying atom decomposition."));            
        } else if (Params::Inst()->scattering.interference.type == "all") {
    	    Info::Inst()->write(string("All interference scattering detected. Applying frame decomposition."));            
        } else {
    	    Err::Inst()->write(string("Scattering Interference type not understood. Must be 'self' or 'all'."));            
            throw;
        }

        // for coherent scattering:
        size_t NAF = NA;
        if (Params::Inst()->scattering.interference.type == "self") {
            NAF = NA;
        } else if (Params::Inst()->scattering.interference.type == "all") {
            NAF = NF;
        }

        size_t ELBYTESIZE = 24;
        size_t NMAXBYTESIZE = Params::Inst()->limits.memory.data;
        if (Params::Inst()->scattering.interference.type == "self") {
            ELBYTESIZE=24; // 3 times double: xyz
        } else if (Params::Inst()->scattering.interference.type == "all") {
            ELBYTESIZE=NA*24; // 3 times double times number of atoms = frame
        }
        
		DecompositionPlan dplan(NN,NQ,NAF,ELBYTESIZE,NMAXBYTESIZE);
		colors = dplan.colors();        

        finalqvectors = H5FQTInterface::get_qvectors(Params::Inst()->scattering.data.file,qindexes);
        partitions = dplan.partitions();
    }

    broadcast(all_comm,partitions,0);

    //broadcast colors
	broadcast(all_comm,reinterpret_cast<size_t*>(&colors[0]),colors.size(),0);
    
    boost::mpi::communicator partition_comm = all_comm.split(colors[all_comm.rank()]);

    size_t nfinalqvectors= finalqvectors.size();
	//broadcast qindexes,qvectors
    broadcast(all_comm,nfinalqvectors,0);

    if (all_comm.rank()!=0) {
    	finalqvectors.resize(nfinalqvectors);
    	qindexes.resize(nfinalqvectors);
    }
	broadcast(all_comm,&(finalqvectors[0].x),nfinalqvectors*3,0);
	broadcast(all_comm,&(qindexes[0]),nfinalqvectors,0);
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

    ////////////////////////////////////////////////////////////
    // Creating of scattering devices
    ////////////////////////////////////////////////////////////

    // all_comm for inter-partition communication, parition_comm for intra-partition communication

    if (Params::Inst()->scattering.interference.type == "self") {
    	p_ScatterDevice = new SelfVectorsScatterDevice(
    			all_comm,
    			partition_comm,
    			sample,
    			thispartition_QIV,
    			Params::Inst()->scattering.data.file);
    }
    else if (Params::Inst()->scattering.interference.type == "all"){
    	if (Params::Inst()->scattering.average.orientation.type == "vectors") {
        	p_ScatterDevice = new AllVectorsScatterDevice(
        			all_comm,
        			partition_comm,
        			sample,
        			thispartition_QIV,
        			Params::Inst()->scattering.data.file);
    	} else if (Params::Inst()->scattering.average.orientation.type == "vectorsthread") {
        	p_ScatterDevice = new AllVectorsThreadScatterDevice(
        			all_comm,
        			partition_comm,
        			sample,
        			thispartition_QIV,
        			Params::Inst()->scattering.data.file);
    	} else if (Params::Inst()->scattering.average.orientation.type == "multipole") {
    		if (Params::Inst()->scattering.average.orientation.multipole.type == "sphere") {
            	p_ScatterDevice = new AllMSScatterDevice(
            			all_comm,
            			partition_comm,
            			sample,
            			thispartition_QIV,
            			Params::Inst()->scattering.data.file);
    		} else if (Params::Inst()->scattering.average.orientation.multipole.type == "cylinder") {
            	p_ScatterDevice = new AllMCScatterDevice(
            			all_comm,
            			partition_comm,
            			sample,
            			thispartition_QIV,
            			Params::Inst()->scattering.data.file);
    		} else {
    			Err::Inst()->write(string("scattering.average.orientation.multipole.type not understood: ")+Params::Inst()->scattering.average.orientation.multipole.type);
    			throw;
    		}
    	} else if (Params::Inst()->scattering.average.orientation.type == "none") {
        	p_ScatterDevice = new AllVectorsScatterDevice(
        			all_comm,
        			partition_comm,
        			sample,
        			thispartition_QIV,
        			Params::Inst()->scattering.data.file);
    	}
    }    

    if (p_ScatterDevice==NULL) {
    	Err::Inst()->write("Error initializing ScatterDevice");
    	throw;
    }

    return p_ScatterDevice;
}
// end of file
