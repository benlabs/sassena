/** \file
This file contains a class which generates the scattering device based on user input provides through the CONTROL module. The constructing routine only return a pointer to the interface.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/

// direct header
#include "scatter_devices/scatter_device_factory.hpp"

// other headers
#include "control.hpp"
#include "log.hpp"
#include "decomposition/assignment.hpp"
#include "decomposition/decomposition_plan.hpp"
#include "scatter_devices/all_vectors_scatter_device.hpp"
#include "scatter_devices/self_vectors_scatter_device.hpp"
#include "scatter_devices/multipole_scatter_device.hpp"

using namespace std;

IScatterDevice* ScatterDeviceFactory::create(
		boost::mpi::communicator& scatter_comm,
		Sample& sample,
        boost::asio::ip::tcp::endpoint fileservice_endpoint,
        boost::asio::ip::tcp::endpoint monitorservice_endpoint,        
		std::vector<CartesianCoor3D>& qvectors)
{
    
    IScatterDevice* p_ScatterDevice = NULL;

    size_t NN = scatter_comm.size();
    size_t NF = sample.coordinate_sets.size();
    string target = Params::Inst()->stager.target;
    size_t NA = sample.atoms.selections[target]->size();
    size_t NQ = qvectors.size();
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

    size_t partitions;
    size_t partitionsize;
    ////////////////////////////////////////////////////////////
    // Decomposition
    ////////////////////////////////////////////////////////////
    // for coherent scattering:
    size_t NAF = NA;
    if (Params::Inst()->scattering.type == "self") {
        NAF = NA;
    } else if (Params::Inst()->scattering.type == "all") {
        NAF = NF;
    }
        
    if (scatter_comm.rank()==0) {

        // decompose parallel space into independent partitions
        // each operating on a distinct set of qvectors.
        
    	Info::Inst()->write(string("Searching for decomposition plan: "));
	    Info::Inst()->write(string("nodes    = ")+ boost::lexical_cast<string>(scatter_comm.size()));
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

        size_t ELBYTESIZE = NF*3*sizeof(coor_t);
        size_t NMAXBYTESIZE = Params::Inst()->limits.stage.memory.data;
        if (Params::Inst()->scattering.type == "self") {
            ELBYTESIZE=NF*3*sizeof(coor_t); // 3 times coor_t times number of frames 
        } else if (Params::Inst()->scattering.type == "all") {
            ELBYTESIZE=NA*3*sizeof(coor_t); // 3 times coor_t times number of atoms = frame
        }
        
		DecompositionPlan dplan(NN,NQ,NAF,ELBYTESIZE,NMAXBYTESIZE);

        partitions = dplan.partitions();
        partitionsize = dplan.partitionsize();
    }

    broadcast(scatter_comm,&partitions,1,0);
    broadcast(scatter_comm,&partitionsize,1,0);
    
    size_t allcommsize = partitions*partitionsize;
    size_t scattercommsize = scatter_comm.size();
    size_t sparenodes = scattercommsize-allcommsize;

    if (sparenodes>0) {
        if (scatter_comm.rank()==0) {
            Warn::Inst()->write(string("Partitioning only partly successful. Number of nodes NOT used: ")+boost::lexical_cast<string>(sparenodes));
        }
    }
    
    size_t allcommflag = 0;
    if (scatter_comm.rank()<allcommsize) allcommflag = 1;
    boost::mpi::communicator all_comm = scatter_comm.split( allcommflag );
    
    if (allcommflag==0) {
        return NULL;
    }
    
    vector<size_t> partitionIDs(allcommsize);
    for(size_t i = 0; i < allcommsize; ++i)
    {
        partitionIDs[i]=(i*partitions)/allcommsize;
    }
    
	// determine the partition this node lives in:
	size_t partitionID = partitionIDs[all_comm.rank()];
    
    boost::mpi::communicator partition_comm = all_comm.split(partitionID);
    
    DivAssignment qindex_assignment(partitions,partitionID,qvectors.size());
	vector<CartesianCoor3D> thispartition_QIV;

	for(size_t i=0;i<qindex_assignment.size();i++) {
		thispartition_QIV.push_back(qvectors[qindex_assignment[i]]);
	}
		
	/////////////////////////////////////////////////////////////
    // Creating of scattering devices
    ////////////////////////////////////////////////////////////

    // all_comm for inter-partition communication, parition_comm for intra-partition communication
    if (Params::Inst()->scattering.type == "self") {
        if (scatter_comm.rank()==0) Info::Inst()->write("Initializing Scatter Device, Vectors (self)");
    	p_ScatterDevice = new SelfVectorsScatterDevice(
    			all_comm,
    			partition_comm,
    			sample,
    			thispartition_QIV,
    			NAF,
		        fileservice_endpoint,
		        monitorservice_endpoint);
    }
    else if (Params::Inst()->scattering.type == "all"){
    	if (Params::Inst()->scattering.average.orientation.type == "vectors") {
            if (scatter_comm.rank()==0) Info::Inst()->write("Initializing Scatter Device, Vectors (all)");
        	p_ScatterDevice = new AllVectorsScatterDevice(
        			all_comm,
        			partition_comm,
        			sample,
        			thispartition_QIV,
        			NAF,        			
        			fileservice_endpoint,
        			monitorservice_endpoint);
    	} else if (Params::Inst()->scattering.average.orientation.type == "multipole") {
    		if (Params::Inst()->scattering.average.orientation.multipole.type == "sphere") {
                if (scatter_comm.rank()==0) Info::Inst()->write("Initializing Scatter Device, Multipole Sphere");
            	p_ScatterDevice = new MPSphereScatterDevice(
        			all_comm,
        			partition_comm,
        			sample,
        			thispartition_QIV,
        			NAF,        			
        			fileservice_endpoint,
        			monitorservice_endpoint);
    		} else if (Params::Inst()->scattering.average.orientation.multipole.type == "cylinder") {
                if (scatter_comm.rank()==0) Info::Inst()->write("Initializing Scatter Device, Multipole Cylinder");
            	p_ScatterDevice = new MPCylinderScatterDevice(
        			all_comm,
        			partition_comm,
        			sample,
        			thispartition_QIV,
        			NAF,        			
        			fileservice_endpoint,
        			monitorservice_endpoint);
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
        			NAF,
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
