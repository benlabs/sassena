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
#include "control.hpp"
#include "decomposition/decompose.hpp"
#include "measures/center_of_mass.hpp"
#include "sample.hpp"


#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
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
	
		Info::Inst()->write("Exchanging sample, database & params information with compute nodes... ");
	}

    if (world.rank()==0) {
        Info::Inst()->write("reading params");
        params->init(string(argv[1]));

        Info::Inst()->write("reading database");
        database->init(string(argv[2]));

        sample.init();
   
		broadcast(world,*params,0);
		broadcast(world,*database,0);
	    broadcast(world,sample,0);	        
    } else {
	    world.recv(0,boost::mpi::any_tag, *params);			
		world.recv(0,boost::mpi::any_tag, *database);			
        world.recv(0,boost::mpi::any_tag, sample);	
    }


		//------------------------------------------//
		//
		// Preparation, Analysis of the system
		//
		//------------------------------------------//
	
		//------------------------------------------//
		//
		// Communication of the sample
		// At this point it is ILLEGAL to change anything within the sample.
		//
		//------------------------------------------//

        size_t NN = world.size();
        size_t NF = sample.coordinate_sets.size();
        
        sample.coordinate_sets.set_selection(sample.atoms.system_selection);
    	sample.coordinate_sets.set_representation(CARTESIAN);	
        
        EvenDecompose edecomp(NF,NN);
        vector<size_t> myframes = edecomp.indexes_for(world.rank());
        size_t NMYF = myframes.size();
        
        vector<CartesianCoor3D> mycofm;
        for(size_t fi = 0; fi < NMYF; ++fi)
        {
            CoordinateSet& cs = sample.coordinate_sets.load(myframes[fi]);
            CartesianCoor3D cofm = CenterOfMass(sample.atoms,sample.atoms.system_selection,sample.coordinate_sets.get_selection(),cs);
            mycofm.push_back(cofm);
        }
        
        std::list< boost::mpi::request > requests;

        vector< vector<CartesianCoor3D> > allcofm(NN);        
        if (world.rank()==0) {

            for(size_t ni = 0; ni < NN; ++ni)
            {
                requests.push_back( world.irecv(ni,0,allcofm[ni]) );                
            }
        }

        requests.push_back( world.isend(0,0,mycofm) );                
        
        boost::mpi::wait_all(requests.begin(),requests.end());

        if (world.rank()==0) {

            vector<CartesianCoor3D> allcofmdelta;
            size_t count=0;
            CartesianCoor3D* p_oldcofm =NULL;
            CartesianCoor3D* p_newcofm =NULL;
            for(size_t i = 0; i < allcofm.size(); ++i)
            {
                for(size_t j = 0; j < allcofm[i].size(); ++j)
                {
                    p_newcofm = &(allcofm[i][j]);

                    if (count!=0) {
                        allcofmdelta.push_back((*p_newcofm)-(*p_oldcofm));
                    }

                    p_oldcofm = p_newcofm;
                    count++;
                }
            }

            ofstream ofile("cofm.txt");
            
            for(size_t i = 0; i < allcofm.size(); ++i)
            {
                for(size_t j = 0; j < allcofm[i].size(); ++j)
                {
                    ofile << allcofm[i][j].x << "\t" << allcofm[i][j].y << "\t" << allcofm[i][j].z << "\t" << endl;
                }
            }
            
            ofstream ofiledelta("cofm-delta.txt");
            ofstream ofileabsdelta("cofm-absdelta.txt");
            
            for(size_t i = 0; i < allcofmdelta.size(); ++i)
            {
                ofiledelta << allcofmdelta[i].x << "\t"<< allcofmdelta[i].y << "\t"<< allcofmdelta[i].z << "\t" << endl;
                ofileabsdelta << allcofmdelta[i].length() << endl;                
            }
                        
        }

	return 0;
}

// end of file
