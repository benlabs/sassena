/*
 *  sassena.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
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
#include <boost/regex.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/filesystem.hpp>
#include <boost/mpi.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>

// other headers
#include "analysis.hpp"
#include "coor3d.hpp"
#include "decompose.hpp"
#include "decomposition_plan.hpp"
#include "control.hpp"
#include "performance_analyzer.hpp"
#include "progress_reporter.hpp"
#include "sample/sample.hpp"
#include "scatter_devices.hpp"
#include "scatter_spectrum.hpp"
#include "timer.hpp"

using namespace std;

int main(int argc,char** argv) {

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
	Params::Inst();
	
	Info::Inst()->set_prefix(to_s(world.rank())+string(".Info>>"));
	Warn::Inst()->set_prefix(to_s(world.rank())+string(".Warn>>"));
	 Err::Inst()->set_prefix(to_s(world.rank())+string(".Err>>"));
	
	Params* params = Params::Inst();
	Database* database = Database::Inst();

	Sample sample;
	
	Timer timer;
	timer.start("total");
	
	if (world.rank()==0) {
		//------------------------------------------//
		//
		// Some welcome message....
		//
		//------------------------------------------//	
	
		Info::Inst()->write("This software is being developed by Benjamin Lindner.                    ");
		Info::Inst()->write("For help, suggestions or correspondense use:                             ");
		Info::Inst()->write("ben@benlabs.net, Benjamin Lindner (Main Developer, Impl. & Maintenance)  ");
		Info::Inst()->write("franc@cmm.ki.si, Franci Merzel (Methodology)                             ");
		Info::Inst()->write("For publications include the following references:                       ");
		Info::Inst()->write(".........................................................................");
		Info::Inst()->write("1. Sassena - Elastic Scattering Calculations on Parallel Computers       ");
		Info::Inst()->write("   to be published                                                       ");		
		Info::Inst()->write(".........................................................................");
		Info::Inst()->write("");

		//------------------------------------------//
		//
		// Test the mpi environment
		//
		//------------------------------------------//	
	
		//------------------------------------------//
		//
		// Setting up the sample
		//
		//------------------------------------------//	
	
		timer.start("sample::setup");
	    
	    Info::Inst()->write("reading params");
        params->init(string(argv[1]));

        Info::Inst()->write("reading database");
        database->init(string(argv[2]));
	    
        sample.init();
	    
		timer.stop("sample::setup");

		//------------------------------------------//
		//
		// Preparation, Analysis of the system
		//
		//------------------------------------------//
	
		Info::Inst()->write(string("Set background scattering length density set to ")+to_s(Params::Inst()->scattering.background.factor));
		
		//------------------------------------------//
		//
		// Communication of the sample
		// At this point it is ILLEGAL to change anything within the sample.
		//
		//------------------------------------------//

		// before calculating anything we need to communicate the sample to any node
		// this is a one-to-all communication

		sample.coordinate_sets.clear_cache(); // reduce overhead

		Info::Inst()->write("Exchanging sample, database & params information with compute nodes... ");

		world.barrier();

		timer.start("sample::communication");

		broadcast(world,*params,0);

		world.barrier();

		broadcast(world,*database,0);

		world.barrier();

		broadcast(world,sample,0);

		world.barrier();
	
		timer.stop("sample::communication");
		
	}
	else {
		
		world.barrier();
	
		world.recv(0,boost::mpi::any_tag, *params);			

		world.barrier();
	
		world.recv(0,boost::mpi::any_tag, *database);			

		world.barrier();
	
		world.recv(0,boost::mpi::any_tag, sample);	

		world.barrier();
		
	}


	//------------------------------------------//
	//
	// Scattering calculation
	//
	//------------------------------------------//

	// prepare an array for all qqqvectors here:
	std::vector<CartesianCoor3D>& qvectors = params->scattering.qvectors;

	if (world.rank()==0) {
		Info::Inst()->write(string("Searching for decomposition plan: "));
		Info::Inst()->write(string("nodes    = ")+ to_s(world.size()));
		Info::Inst()->write(string("qvectors = ")+ to_s(qvectors.size()));
		Info::Inst()->write(string("frames   = ")+ to_s(sample.coordinate_sets.size()));	
	}
	
	vector<size_t> frames;
		
	for(size_t i = 0; i < sample.coordinate_sets.size(); ++i)
	{
		frames.push_back(i);
	}
	DecompositionPlan dplan(world,qvectors,frames);
	if (world.rank()==0) {	
		Info::Inst()->write(string("Decomposition has ")+to_s(dplan.worlds())+string(" partitions"));
		Info::Inst()->write(string("Static imbalance factor (0 is best): ")+ to_s(dplan.static_imbalance()));
	}
	
	boost::mpi::communicator local = dplan.split();
	
	vector<CartesianCoor3D> myqvectors = dplan.qvectors();
	if (world.rank()==0) {
		Info::Inst()->write(string("Scattering target selection: ")+params->scattering.target);
	}

	ScatterDevice* p_ScatterDevice = NULL;

	if (params->scattering.interference.type == "self") {
		p_ScatterDevice = new SelfScatterDevice(local,sample);
	}
	else if (params->scattering.interference.type == "all"){
		if (params->scattering.average.orientation.type == "vectors") {
			p_ScatterDevice = new AllScatterDevice(local,sample);			
		} else if (params->scattering.average.orientation.type == "multipole") {
			//Err::Inst()->write("Multipole method currently not supported");
			//throw;
			if (params->scattering.average.orientation.multipole.type == "sphere") {
				p_ScatterDevice = new AllMSScatterDevice(local,sample);			
			} else if (params->scattering.average.orientation.multipole.type == "cylinder") {
				p_ScatterDevice = new AllMCScatterDevice(local,sample);			
			} else {
				Err::Inst()->write(string("scattering.average.orientation.multipole.type not understood: ")+params->scattering.average.orientation.multipole.type);
				throw;
			}
		} else if (params->scattering.average.orientation.type == "exact") {
			p_ScatterDevice = new AllExactScatterDevice(local,sample);					    
		} else if (params->scattering.average.orientation.type == "none") {
			p_ScatterDevice = new AllScatterDevice(local,sample);					    
		}
	}
	
	if (p_ScatterDevice==NULL) {
		Err::Inst()->write("Error initializing ScatterDevice");
	}

	ScatterSpectrum scatter_spectrum;

	if (world.rank()==0) { // only let the world head report...
		Info::Inst()->write("Progress monitors head node only.");
	}
	ProgressReporter progress_reporter(myqvectors.size(),0); // this acts like a progress bar
	for(size_t qi = 0; qi < myqvectors.size(); ++qi)
	{			

			// progress indicator
			if (world.rank()==0) { // only let the world head report...
				progress_reporter.set(qi);
				progress_reporter.report(); // only reports if significant progress has been made!
			}
			
			timer.start("scatter");
			
			p_ScatterDevice->execute(myqvectors[qi]);
			
			if (local.rank()==0) {
				// results are aggregated in 0-node
				scatter_spectrum.push_back( make_pair(myqvectors[qi],p_ScatterDevice->get_spectrum() ) ) ;
			}
					
			timer.stop("scatter");
	}

	if (world.rank()==0) {
		Info::Inst()->write(string("Waiting for all nodes to finish calculations..."));		
	}
	world.barrier();	
	
	//------------------------------------------//
	//
	// Output
	//
	//------------------------------------------//	
	
	// now all local.rank()==0 have aggregated the spectra
	// we have to aggregrate it on world scope as well

	size_t headindicator = 0;
	if (local.rank()==0) headindicator = 1;
	// do another split, based on the headindicator
	boost::mpi::communicator headcomm = world.split(headindicator);
	ScatterSpectrum final_scatter_spectrum;
	if (headindicator==1) {
		// gather results on head of heads
		vector<ScatterSpectrum> scatspecs;
		Info::Inst()->write(string("Aggregrating spectra, local count:")+to_s(scatter_spectrum.size()));

		boost::mpi::gather(headcomm,scatter_spectrum,scatspecs,0);
		
		if (headcomm.rank()==0) {
			ScatterSpectrum fss(scatspecs); // fss = final scatter spectrum

			timer.start("result::output");
			std::string prefix = Params::Inst()->output.prefix;
			for(std::vector<OutputFileParameters>::iterator ofi = Params::Inst()->output.files.begin(); ofi != Params::Inst()->output.files.end(); ++ofi)
			{
				Info::Inst()->write(string("Writing to: ") + ofi->filename);
				if (ofi->method=="plain") {
					fss.write_plain(ofi->filename,ofi->format);				
				} else if (ofi->method=="average") {
					fss.write_average(ofi->filename,ofi->format);
				} else {
					Warn::Inst()->write(string("Output method not understood: ")+ofi->method);
				}
			}

			timer.stop("result::output");

		}
	}

	// wait before deleting anything...
	world.barrier();

	if (world.rank()==0) {
		Info::Inst()->write(string("Aggregating timing information for performance analysis..."));		
	}

	PerformanceAnalyzer perfanal(world,p_ScatterDevice->timer); // collect timing information from everybody.
	
	delete p_ScatterDevice;

	// make nodes empty, computation finished
	sample.coordinate_sets.clear_cache();


	//------------------------------------------//
	//
	// Finished
	//
	//------------------------------------------//	
	
	timer.stop("total");
		
	if (world.rank()==0) {
		perfanal.report();
		perfanal.report_relative(timer.sum("total"));
		
		Info::Inst()->write(string("Total runtime (s): ")+to_s(timer.sum("total")));
		Info::Inst()->write("Successfully finished... Have a nice day!");
	}

	return 0;
}

// end of file
