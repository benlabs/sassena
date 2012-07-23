/** \file
This file is the main executable for the software sassena.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/

/** @mainpage Project: Sassena
\par Description:
Sassena is a software to compute X-ray and Neutron Scattering Intensities from Molecular Dynamics Simulation Trajectories. It's purpose is to allow high performance scattering calculations using massively parallel computers which feature superior network infrastructures, e.g. infiniband.

\par License:
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as
published by the Free Software Foundation; either version 2 of
the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details at
http://www.gnu.org/copyleft/gpl.html
If you decide to redistribute modifications, please follow the
scientific honor code and refer to any derived work appropriately.
Also consider to consult with the original developers first, before
you start to fork the project.

\par Further Information:
- If you find bugs file them at http://bugzilla.sassena.org
- If you want to consult the user community goto http://forum.sassena.org
- If you want to download support material goto http://www.sassena.org

\par Contact:
- Benjamin Lindner <ben@benlabs.net> (Original Developer)
*/


// direct header
#include "common.hpp"

// standard header
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <sys/time.h>
#include <vector>

// special library headers
#include <boost/asio.hpp>
#include <boost/date_time.hpp>
#include <boost/regex.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/exception/all.hpp>
#include <boost/filesystem.hpp>
#include <boost/mpi.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>

// other headers
#include "exceptions/exceptions.hpp"
#include "math/coor3d.hpp"
#include "decomposition/decomposition_plan.hpp"
#include "control.hpp"
#include "log.hpp"
#include "report/performance_analyzer.hpp"
#include "report/timer.hpp"
#include "mpi/wrapper.hpp"
#include "sample/sample.hpp"
#include "scatter_devices/scatter_device_factory.hpp"
#include "services.hpp"

#include "SassenaConfig.hpp"

using namespace std;

void print_title() {
    
	Info::Inst()->write("This software is being developed by Benjamin Lindner.                    ");
	Info::Inst()->write("For help, suggestions or correspondense use:                             ");
	Info::Inst()->write("ben@benlabs.net, Benjamin Lindner (Main Developer, Impl. & Maintenance)  ");
	Info::Inst()->write("franc@cmm.ki.si, Franci Merzel (Methodology)                             ");
	Info::Inst()->write("For publications include the following references:                       ");
	Info::Inst()->write(".........................................................................");
	Info::Inst()->write("1. Sassena - X-ray and Neutron Scattering Calculated                     ");
	Info::Inst()->write("  from Molecular Dynamics Trajectories using Massively Parallel Computers");		
	Info::Inst()->write("  published as a Computer Program in Physics Paper                       ");		
	Info::Inst()->write(".........................................................................");
    Info::Inst()->write(string("Version Information: ") + string(Sassena_VERSIONSTRING));
	Info::Inst()->write("");

}
void print_description() {
    Info::Inst()->write(".................................................................");
	Info::Inst()->write("......................D.E.S.C.R.I.P.T.I.O.N......................");
	Info::Inst()->write(".................................................................");	

	Info::Inst()->write("This binary computes the scattering intensities directly from");	
	Info::Inst()->write("a molecular dynamics trajectory. ");	
	Info::Inst()->write(".................................................................");	
}



void print_initialization() {
    Info::Inst()->write(".................................................................");
	Info::Inst()->write("...................I.N.I.T.I.A.L.I.Z.A.T.I.O.N...................");
	Info::Inst()->write(".................................................................");	
}

void print_computation() {
    Info::Inst()->write(".................................................................");
	Info::Inst()->write("......................C.O.M.P.U.T.A.T.I.O.N......................");
	Info::Inst()->write(".................................................................");	
}

void print_analysis() {
    Info::Inst()->write(".................................................................");
	Info::Inst()->write(".....................R.U.N...A.N.A.L.Y.S.I.S.....................");
	Info::Inst()->write(".................................................................");	
}

/** 
Entry point for the software sassena
*/
int main(int argc,char* argv[]) {

	//------------------------------------------//
	//
	// MPI Initialization
	//
	//------------------------------------------//
  	boost::mpi::environment env(argc, argv);
  	boost::mpi::communicator world;

    // The rank 0 node is responsible for the progress output and to inform the user
    // compute nodes should be silent all the time, except when errors occur.
    // In that case the size and the rank should be included into the error message in the following way:

	// make sure any singleton class exists:
	Info::Inst();
	Err::Inst();
	Warn::Inst();
	
	Info::Inst()->set_prefix(boost::lexical_cast<string>(world.rank())+string(".Info>>"));
	Warn::Inst()->set_prefix(boost::lexical_cast<string>(world.rank())+string(".Warn>>"));
	Err::Inst()->set_prefix(boost::lexical_cast<string>(world.rank())+string(".Err>>"));
	
	Params* params = Params::Inst();
	Database* database = Database::Inst();

	Sample sample;
	
	Timer timer;
	timer.start("total");

    if (world.rank()==0) {
        print_title();
        print_description();
    }

    bool initstatus = true;

    if (world.rank()==0) print_initialization();
	if (world.rank()==0) {
	    
		try {

    		timer.start("sample::setup");
	    
            params->init(argc,argv);
            database->init();
	    
            sample.init();
		
		    timer.stop("sample::setup");
		} catch (sassena::terminate_request const& e) {
            Info::Inst()->write("Hangup requested");
            initstatus = false; 		    
        } catch (boost::exception const& e ) {
            initstatus = false; 
            Err::Inst()->write("Caught BOOST error, sending hangup to all nodes");
            stringstream ss; ss << diagnostic_information(e);
            Err::Inst()->write(string("Diagnotic information: ") + ss.str());
        } catch (std::exception const& e) {
            initstatus = false; 
            Err::Inst()->write("Caught STD error, sending hangup to all nodes");
            Err::Inst()->write(string("what() : ") + e.what());
        } catch (...) {
            initstatus = false; 
            Err::Inst()->write("Caught error: UNKNOWN sending hangup to all nodes");
        }
    }
    
    broadcast(world,&initstatus,1,0);            	
	// if something went wrong during initialization, exit now.
    if (!initstatus) {
        world.barrier();
        return 0;
    }
    
	if (world.rank()==0) Info::Inst()->write(string("Set background scattering length density set to ")+boost::lexical_cast<string>(Params::Inst()->scattering.background.factor));
		
	//------------------------------------------//
	//
	// Communication of the sample
	// At this point it is ILLEGAL to change anything within the sample.
	//
	//------------------------------------------//

	if (world.rank()==0) Info::Inst()->write("Exchanging sample, database & params information with compute nodes... ");
    
	world.barrier();

	timer.start("sample::communication");
    if (world.rank()==0) Info::Inst()->write("params... ");
    mpi::wrapper::broadcast_class<Params>(world,*params,0);
	world.barrier();
    if (world.rank()==0) Info::Inst()->write("database... ");
    mpi::wrapper::broadcast_class<Database>(world,*database,0);
	world.barrier();
    if (world.rank()==0) Info::Inst()->write("sample... ");
    mpi::wrapper::broadcast_class<Sample>(world,sample,0);
	world.barrier();
	timer.stop("sample::communication");

	//------------------------------------------//
	//
	// Scattering calculation
	//
	//------------------------------------------//

    if (world.rank()==0) print_computation();

    boost::mpi::communicator scatter_comm = boost::mpi::communicator(world,world.group());

	world.barrier();

    if (world.rank()==0) Info::Inst()->write("Setting up services...");	
    
    boost::asio::io_service io_service;
    
    // setup hdf5writer service and broadcast service to clients
    HDF5WriterService* p_hdf5writer = NULL;    
    MonitorService* p_monitorservice = NULL;    

    // prepare services
    if (world.rank()==0) {
        Info::Inst()->write("Initializing data file service...");	
        p_hdf5writer = new HDF5WriterService(io_service, Params::Inst()->scattering.signal.filepath,sample.coordinate_sets.size());
        Info::Inst()->write("Initializing monitor service...");	
        p_monitorservice = new MonitorService(io_service, 0, scatter_comm.size());        
    }
    
    // start services
    if (world.rank()==0) {
        Info::Inst()->write("Starting data file service...");	
        p_hdf5writer->run();
        Info::Inst()->write("Starting monitor service...");	
        p_monitorservice->run();
        Info::Inst()->write("Services setup and running...");	
    }

    // only compute those q vectors which have not been written so far
    // shutdown gracefully if no work is to be done
    initstatus = true;
    vector<CartesianCoor3D> qvectors;
    if (world.rank()==0) {
        
        Info::Inst()->write(string("Checking ")+Params::Inst()->scattering.signal.filepath + string(" for old results."));
        vector<CartesianCoor3D> oldqvectors = p_hdf5writer->get_qvectors();
        set<CartesianCoor3D> oldqvectors_set;
        for(size_t i = 0; i < oldqvectors.size(); ++i)
        {
            oldqvectors_set.insert(oldqvectors[i]);
        }

        for(size_t i = 0; i < Params::Inst()->scattering.qvectors.size(); ++i)
        {
            if (oldqvectors_set.find(Params::Inst()->scattering.qvectors[i])==oldqvectors_set.end()) {
                qvectors.push_back(Params::Inst()->scattering.qvectors[i]);
            }
        }
        
        if (qvectors.size()!=Params::Inst()->scattering.qvectors.size()) {
            Info::Inst()->write(string("Found ")+ boost::lexical_cast<string>(Params::Inst()->scattering.qvectors.size()-qvectors.size()) + string(" old qvectors."));            
            Info::Inst()->write(string("Total number of qvectors to compute reduced from ")+boost::lexical_cast<string>(Params::Inst()->scattering.qvectors.size()) + string(" to ") + boost::lexical_cast<string>(qvectors.size()) + string("."));            
        }

        if (qvectors.size()==0) {
            Info::Inst()->write("Very good. No qvectors left to compute.");

            Info::Inst()->write("Stopping services.");
            p_monitorservice->hangup();
            p_hdf5writer->hangup();

            initstatus = false;
            Info::Inst()->write("Broadcasting hangup to all nodes. Have a nice day!");
        }
    }
    broadcast(world,&initstatus,1,0);            	
    if (!initstatus) return 1;
    
    // communicate some global parameters
    size_t nq=qvectors.size();
    broadcast(scatter_comm,&nq,1,0);
    if (scatter_comm.rank()!=0) {
        qvectors.resize(nq);
    }
    coor2_t* p_qvectors = &(qvectors[0].x);
    broadcast(scatter_comm,p_qvectors,3*nq,0);
       
        
    // communicate services to all nodes
    std::string host_str,fileservice_port_str,monitorservice_port_str;
    if (scatter_comm.rank()==0) {
        host_str = boost::asio::ip::host_name();
        fileservice_port_str = boost::lexical_cast<string>(p_hdf5writer->get_endpoint().port());
        monitorservice_port_str = boost::lexical_cast<string>(p_monitorservice->get_endpoint().port());
    }    
    if (world.rank()==0) Info::Inst()->write("Broadcasting service information");	

    mpi::wrapper::broadcast_class<std::string>(scatter_comm,host_str,0);
    mpi::wrapper::broadcast_class<std::string>(scatter_comm,fileservice_port_str,0);
    mpi::wrapper::broadcast_class<std::string>(scatter_comm,monitorservice_port_str,0);

    if (world.rank()==0) {
        Info::Inst()->write(string("Server host name    : ")+host_str);	
        Info::Inst()->write(string("FileService Port    : ")+fileservice_port_str);	
        Info::Inst()->write(string("MonitorService Port : ")+monitorservice_port_str);	
    }


    boost::asio::ip::tcp::resolver resolver(io_service);
    //fileserver
    boost::asio::ip::tcp::resolver::query query(boost::asio::ip::tcp::v4(),host_str,fileservice_port_str,boost::asio::ip::resolver_query_base::numeric_service);
    boost::asio::ip::tcp::resolver::iterator it = resolver.resolve(query);
    boost::asio::ip::tcp::endpoint fileservice_endpoint(*it);
    //monitor
    boost::asio::ip::tcp::resolver::query query2(boost::asio::ip::tcp::v4(),host_str,monitorservice_port_str,boost::asio::ip::resolver_query_base::numeric_service);
    it = resolver.resolve(query2);
    boost::asio::ip::tcp::endpoint monitorservice_endpoint(*it);
    
    if (world.rank()==0) Info::Inst()->write("Setting up parallel environment...");	
    
    // now create scattering device w/ reference to file server
    IScatterDevice* p_ScatterDevice =NULL;
    p_ScatterDevice = ScatterDeviceFactory::create(
    		scatter_comm,
    		sample,
    	    fileservice_endpoint,
            monitorservice_endpoint,
            qvectors);

	world.barrier();	

    // start computation tasks
    if (world.rank()==0) Info::Inst()->write("Starting scattering...");

	try {
	    if (p_ScatterDevice!=NULL) p_ScatterDevice->run();		
	} catch(sassena::terminate_request) {
		Err::Inst()->write("Terminatation requested.");
		// do a clean hangup here...
		throw;
	}

    world.barrier();
    
    if (world.rank()==0) {
        // grace period for the monitor thread to update last progress
        boost::this_thread::sleep(boost::posix_time::milliseconds(25));
        Info::Inst()->write("Scattering finished...");
    }
    
    // shutdown services
    if (scatter_comm.rank()==0) {
        p_monitorservice->hangup();
        p_hdf5writer->hangup();
    }
    
    world.barrier();

    if (world.rank()==0) print_analysis();
    
	if (world.rank()==0) {
		Info::Inst()->write(string("Aggregating timing information for performance analysis..."));		
	}
	
	timer.stop("total");
	
    std::map<boost::thread::id,Timer> performance_timer;
    if (p_ScatterDevice!=NULL) performance_timer = p_ScatterDevice->getTimer();
    PerformanceAnalyzer perfanal(world,performance_timer); // collect timing information from everybody.
    if (p_ScatterDevice!=NULL) delete p_ScatterDevice;
        
    if (world.rank()==0) {
		perfanal.report();
		perfanal.report_relative(timer.sum("total")*world.size());
    	Info::Inst()->write(string("Total runtime (s): ")+boost::lexical_cast<string>(timer.sum("total")));
    	Info::Inst()->write("Successfully finished... Have a nice day!");
	}

	//------------------------------------------//
	//
	// Finished
	//
	//------------------------------------------//	

	return 0;
}


// end of file
