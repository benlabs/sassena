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
#include "log.hpp"
#include "parameters.hpp"
#include "database.hpp"
#include "performance_analyzer.hpp"
#include "progress_reporter.hpp"
#include "sample.hpp"
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

	Sample sample;

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

	    // create the sample via structure file	
		Info::Inst()->write(string("Reading structure information from file: ")+params->sample.structure.file);
		sample.add_atoms(params->sample.structure.file,params->sample.structure.format);
		Info::Inst()->write(string("Done. Atoms read: ")+to_s(sample.atoms.size()));

	   	// add selections / groups
		for(map<string,SampleGroupParameters>::iterator sgpi = params->sample.groups.begin();sgpi!=params->sample.groups.end();sgpi++)
		{
			SampleGroupParameters& sp = sgpi->second;
			Info::Inst()->write(string("Reading selection file: ")+sp.file);
			sample.atoms.add_selection(sp.name,sp.file,sp.format,sp.select,sp.select_value);
		}

		// add custom selections here ( if not already set! )
		if (sample.atoms.selections.find("system")==sample.atoms.selections.end()) {
			// this shortcut creates a full selection
			sample.atoms.add_selection("system",true);		
		}
	
		// apply deuteration
		for(size_t i = 0; i < params->sample.deuter.size(); ++i)
		{
			sample.deuter(params->sample.deuter[i]);
		}
			
		// read in frame information
		for(size_t i = 0; i < params->sample.frames.size(); ++i)
		{
			SampleFramesetParameters& f = params->sample.frames[i];
			Info::Inst()->write(string("Reading frames from: ")+f.filename);
			size_t nof = sample.frames.add_frameset(f.filename,f.type,f.first,f.last,f.last_set,f.stride,sample.atoms);			
			Info::Inst()->write(string("Found ")+to_s(nof)+string(" frames"));			
		}
		Info::Inst()->write(string("Total number of frames found: ")+to_s(sample.frames.size()));
		
		// select wrapping behaviour
		if (params->sample.pbc.wrapping) {
			sample.frames.wrapping=true;
			string centergroup = params->sample.pbc.center;
			sample.atoms.assert_selection(centergroup);
			sample.frames.centergroup_selection = sample.atoms.selections[centergroup];			
		}
		else {
			sample.frames.wrapping = false;
		}

		timer.stop("sample::setup");

		//------------------------------------------//
		//
		// Preparation, Analysis of the system
		//
		//------------------------------------------//
	
		timer.start("sample::preparation");

		if (params->scattering.background.type=="auto") {
			Err::Inst()->write("Automatic background support disabled for now.");
			throw;

//			// determine background scattering density from grid based analysis
//			// necessary:
//			// group = name of the group which accounts for the background
//			// resolution = grid resolution
//			// hydration = hydration layer around each particle, 
//			//    ie. solvent molecules within this grid distance are not counted towards bulk
//			int r     = params->scattering.background.resolution; // resolution of grid
//			double h  = params->scattering.background.hydration; // hydration layer of each particle (not background)
//		
//			Info::Inst()->write("Calculate background scattering length density");
//			Info::Inst()->write(string("using grid resolution ")+to_s(r)+string(" and hydration of ")+to_s(h));
//			pair<double,double>	b0 = Analysis::background_avg(sample,r,h,CartesianCoor3D(0,0,0));
//			Info::Inst()->write(string("Average background scattering length: ")+to_s(b0.first)+string(" +- ")+to_s(b0.second));
//			params->scattering.background.factor = b0.first;
			
		}
		
		sample.background = params->scattering.background.factor;
		Info::Inst()->write(string("Set background scattering length density to value of ")+to_s(sample.background));
		
		
		timer.stop("sample::preparation");
		
		//------------------------------------------//
		//
		// Communication of the sample
		// At this point it is ILLEGAL to change anything within the sample.
		//
		//------------------------------------------//

		// before calculating anything we need to communicate the sample to any node
		// this is a one-to-all communication

		sample.frames.clear_cache(); // reduce overhead

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
		Info::Inst()->write(string("frames   = ")+ to_s(sample.frames.size()));	
	}
	
	vector<size_t> frames;
	for(size_t i = 0; i < sample.frames.size(); ++i)
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
		if (params->scattering.average.orientation.method == "bruteforce") {
			p_ScatterDevice = new AllScatterDevice(local,sample);			
		} else if (params->scattering.average.orientation.method == "multipole") {
			//Err::Inst()->write("Multipole method currently not supported");
			//throw;
			p_ScatterDevice = new AllMScatterDevice(local,sample);			
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
	sample.frames.clear_cache();


	//------------------------------------------//
	//
	// Finished
	//
	//------------------------------------------//	
	
	timer.stop("total");
		
	if (world.rank()==0) {
		perfanal.report();
		Info::Inst()->write(string("Total runtime (s): ")+to_s(timer.sum("total")));
		Info::Inst()->write("Successfully finished... Have a nice day!");
	}

	return 0;
}

// end of file
