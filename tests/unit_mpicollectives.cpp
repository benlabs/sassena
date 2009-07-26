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
#include "log.hpp"
#include "parameters.hpp"
#include "sample.hpp"
#include "scatterdata.hpp"
#include "settings.hpp"
#include "tasks.hpp"
#include "timer.hpp"

using namespace std;

int main (int argc, char **argv)
{
	boost::mpi::environment env(argc, argv);
  	boost::mpi::communicator world;
    // Get the number of processors this job is using:
	int rank = world.rank();

    // Get the rank of the processor this thread is running on.  (Each
    // processor has a unique rank.)
	int size = world.size();
  	

	world.barrier();

	Timer timer;
	timer.start("total");
	
	if (rank==0) clog << "Testing MPI collectives" << endl;
	
	vector<double> myvectors;
	vector<vector<double> > allvectors;

//	if (rank==0) clog << "allgather 10M ..." << endl;
//
//	myvectors.resize(10*1000*1000,rank);
//	for(size_t i = 0; i < 10; ++i)
//	{
//		timer.start("allgather 10M");
//		boost::mpi::all_gather(world,myvectors,allvectors);
//		timer.stop("allgather 10M");
//	}
//
//	if (rank==0) clog << "allgather 1M ..." << endl;
//	myvectors.resize(1000*1000,rank);
//	for(size_t i = 0; i < 10; ++i)
//	{
//		timer.start("allgather 1M");
//		boost::mpi::all_gather(world,myvectors,allvectors);
//		timer.stop("allgather 1M");
//	}
//
//	if (rank==0) clog << "allgather 100k ..." << endl;
//	myvectors.resize(100*1000,rank);
//	for(size_t i = 0; i < 10; ++i)
//	{
//		timer.start("allgather 100k");
//		boost::mpi::all_gather(world,myvectors,allvectors);
//		timer.stop("allgather 100k");
//	}
//
//	if (rank==0) clog << "allgather 10k ..." << endl;
//	myvectors.resize(10*1000,rank);
//	for(size_t i = 0; i < 10; ++i)
//	{
//		timer.start("allgather 10k");
//		boost::mpi::all_gather(world,myvectors,allvectors);
//		timer.stop("allgather 10k");
//	}
//
//	if (rank==0) clog << "allgather 1k ..." << endl;
//	myvectors.resize(1000,rank);
//	for(size_t i = 0; i < 10; ++i)
//	{
//		timer.start("allgather 1k");
//		boost::mpi::all_gather(world,myvectors,allvectors);
//		timer.stop("allgather 1k");
//	}
	
	myvectors.resize(0);
	if (rank==0) clog << "allgather DT 10M ..." << endl;
	
	double* parray = (double*) malloc(sizeof(double)*10*1000*1000);
	for(size_t i = 0; i < 10*1000*1000; ++i)
	{
		parray[i]=rank;
	}
	double* pallarray = (double*) calloc(size*10*1000*1000,sizeof(double));


	for(size_t i = 0; i < 10; ++i)
	{
		timer.start("allgather DT 10M");
		boost::mpi::all_gather(world,parray,10*1000*1000,pallarray);
		timer.stop("allgather DT 10M");
	}
	
	for(size_t i = 0; i < size; ++i)
	{
		for(size_t j = 0; j < 10*1000*1000; ++j)
		{
			if (i!=pallarray[(i*10*1000*1000)+j]) {
				cerr << "PAEEEE"<< endl;
				throw;
			}
		}
	}
	
	free(parray);
	free(pallarray);


	vector<double> vsarray;
	vsarray.resize(10*1000*1000,rank);
	vector<double> vsallarray;
	vsallarray.resize(size*10*1000*1000,0);

	for(size_t i = 0; i < 10; ++i)
	{
		timer.start("allgather VsDT 10M");
		boost::mpi::all_gather(world,&vsarray[0],10*1000*1000,&vsallarray[0]);
		timer.stop("allgather VsDT 10M");
	}
	
	
	for(size_t i = 0; i < size; ++i)
	{
		for(size_t j = 0; j < 10*1000*1000; ++j)
		{
			if (i!=vsallarray[(i*10*1000*1000)+j]) {
				cerr << "SVAEEEE"<< endl;
				throw;
			}
		}
		
	}


	
	
	if (rank==0) {
		vector<string> keys = timer.keys();
				
		clog << "INFO>> " << "                                                                         " << endl;
		clog << "INFO>> " << "                    Performance Analysis                                 " << endl;
		clog << "INFO>> " << "-------------------------------------------------------------------------" << endl;
		clog << "INFO>> " << " mean and total runtimes:                                                " << endl;				
		clog << "INFO>> " << "-------------------------------------------------------------------------" << endl;
		clog << "INFO>> ";
		clog << setw(31) << " measure |";
		clog << setw(12) << " total |";
		clog << setw(10) << " count |";
		clog << setw(10) << " mean |";
		clog << setw(10) << " stddev ";
		clog << endl;
		clog << "INFO>> " << "-------------------------------------------------------------------------" << endl;		
		for (vector<string>::iterator ki=keys.begin();ki!=keys.end();ki++) {
			clog << "INFO>> " << setw(29) << *ki << " |" << "\t" ;
			clog << setiosflags(ios::fixed) << setprecision(3) << setw(8) << timer.sum(*ki)            << " |";	
			clog << setiosflags(ios::fixed) << setprecision(0) << setw(8) << timer.count(*ki)          << " |";		
			clog << setiosflags(ios::fixed) << setprecision(3) << setw(8) << timer.mean(*ki)           << " |";
			clog << setiosflags(ios::fixed) << setprecision(3) << setw(8) << sqrt(timer.variance(*ki)) << endl;
		}
		clog << "INFO>> " << "-------------------------------------------------------------------------" << endl;		

		clog << "INFO>> " << "-------------------------------------------------------------------------" << endl;
		clog << "INFO>> " << " watermarks:                                                " << endl;				
		clog << "INFO>> " << "-------------------------------------------------------------------------" << endl;
		clog << "INFO>> ";
		clog << setw(31) << " measure |";
		clog << setw(12) << " min |";
		clog << setw(10) << " max ";
		clog << endl;
		clog << "INFO>> " << "-------------------------------------------------------------------------" << endl;		
		for (vector<string>::iterator ki=keys.begin();ki!=keys.end();ki++) {
			clog << "INFO>> " << setw(29) << *ki << " |" << "\t" ;
			clog << setiosflags(ios::fixed) << setprecision(3) << setw(8) << timer.min(*ki)            << " |";	
			clog << setiosflags(ios::fixed) << setprecision(3) << setw(8) << timer.max(*ki)            << endl;
		}
		clog << "INFO>> " << "-------------------------------------------------------------------------" << endl;		
	
	}
	
	
	return 0;
}