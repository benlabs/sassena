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
#include "scatter_devices/all_multipole_sphere.hpp"
#include "scatter_devices/all_multipole_cylinder.hpp"
#include "scatter_devices/all_vectors.hpp"
#include "scatter_devices/all_vectorsthread.hpp"
#include "scatter_devices/self_vectors.hpp"

using namespace std;

ScatterDevice* ScatterDeviceFactory::create(boost::mpi::communicator& local,Sample& sample) {
    
    ScatterDevice* p_ScatterDevice = NULL;
    
    if (Params::Inst()->scattering.interference.type == "self") {
    	p_ScatterDevice = new SelfVectorsScatterDevice(local,sample);
    }
    else if (Params::Inst()->scattering.interference.type == "all"){
    	if (Params::Inst()->scattering.average.orientation.type == "vectors") {
    		p_ScatterDevice = new AllVectorsScatterDevice(local,sample);			
    	} else if (Params::Inst()->scattering.average.orientation.type == "vectorsthread") {
        		p_ScatterDevice = new AllVectorsThreadScatterDevice(local,sample);			
    	} else if (Params::Inst()->scattering.average.orientation.type == "multipole") {
    		if (Params::Inst()->scattering.average.orientation.multipole.type == "sphere") {
    			p_ScatterDevice = new AllMSScatterDevice(local,sample);			
    		} else if (Params::Inst()->scattering.average.orientation.multipole.type == "cylinder") {
    			p_ScatterDevice = new AllMCScatterDevice(local,sample);			
    		} else {
    			Err::Inst()->write(string("scattering.average.orientation.multipole.type not understood: ")+Params::Inst()->scattering.average.orientation.multipole.type);
    			throw;
    		}
    	} else if (Params::Inst()->scattering.average.orientation.type == "none") {
    		p_ScatterDevice = new AllVectorsScatterDevice(local,sample);					    
    	}
    }
    
    if (p_ScatterDevice==NULL) {
    	Err::Inst()->write("Error initializing ScatterDevice");
    }
    
    return p_ScatterDevice;
}
// end of file