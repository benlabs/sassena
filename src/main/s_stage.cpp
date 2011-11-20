/** \file
This executable routine is used to utilize the parallel data staging capabilities of the software sassena. It sets the partition size to all available nodes and reads and write any trajectory data in parallel. Can be used to benchmark the IO capacities of the hardware or to rewrite trajectory data using the capabilities offered through sassena.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
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
#include "stager/data_stager.hpp"
#include "control.hpp"
#include "log.hpp"
#include "report/performance_analyzer.hpp"
#include "report/timer.hpp"
#include "mpi/wrapper.hpp"
#include "sample/sample.hpp"
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
	Info::Inst()->write("1. Sassena - Scattering Calculations on Parallel Computers               ");
	Info::Inst()->write("   to be published                                                       ");		
	Info::Inst()->write(".........................................................................");
    Info::Inst()->write(string("Version Information: ") + string(Sassena_VERSIONSTRING));
	Info::Inst()->write("");

}

void print_description() {
    Info::Inst()->write(".................................................................");
	Info::Inst()->write("......................D.E.S.C.R.I.P.T.I.O.N......................");
	Info::Inst()->write(".................................................................");	

	Info::Inst()->write("This binary stages the trajectory data on a parallel partition");	
	Info::Inst()->write("It may be used to test I/O and to generate post-processed");	
	Info::Inst()->write("trajectories by setting stager.dump=true");		
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
    		
	//------------------------------------------//
	//
	// Communication of the sample
	// At this point it is ILLEGAL to change anything within the sample.
	//
	//------------------------------------------//

	if (world.rank()==0) Info::Inst()->write("Exchanging sample, database & params information with compute nodes... ");
    
	world.barrier();

    if (world.rank()==0) Info::Inst()->write("params... ");
	timer.start("params::comm");
    mpi::wrapper::broadcast_class<Params>(world,*params,0);
	world.barrier();
	timer.stop("params::comm");
    if (world.rank()==0) Info::Inst()->write("database... ");
	timer.start("db::comm");
    mpi::wrapper::broadcast_class<Database>(world,*database,0);
	world.barrier();
	timer.stop("db::comm");
    if (world.rank()==0) Info::Inst()->write("sample... ");
	timer.start("sample::comm");
    mpi::wrapper::broadcast_class<Sample>(world,sample,0);
	world.barrier();
	timer.stop("sample::comm");

	//------------------------------------------//
	//
	// Scattering calculation
	//
	//------------------------------------------//
        
    if (world.rank()==0) Info::Inst()->write("Setting up parallel environment...");	

    if (world.rank()==0) Info::Inst()->write("Starting staging...");	
    
    // prepare sample
    size_t NF = sample.coordinate_sets.size();
    string target = Params::Inst()->stager.target;
    size_t NA = sample.atoms.selections[target]->size();
    {
        sample.coordinate_sets.set_selection(sample.atoms.selections[target]);
        sample.coordinate_sets.set_representation(CARTESIAN);        
    }

    Info::Inst()->write(string("stager.mode=")+Params::Inst()->stager.mode);
	
    if (Params::Inst()->stager.mode=="frames") {
        DataStagerByFrame data_stager(sample,world,world,timer);
        coor_t* p_coordinates = data_stager.stage();
        delete p_coordinates;
    } else if (Params::Inst()->stager.mode=="atoms") {
        DataStagerByAtom data_stager(sample,world,world,timer);
        coor_t* p_coordinates = data_stager.stage();
        delete p_coordinates;        
    } else {
        Err::Inst()->write(string("Staging mode not understood stager.mode=")+Params::Inst()->stager.mode);
        Err::Inst()->write(string("Use 'frames' or 'atoms'"));
        throw;
    }
    
    if (world.rank()==0) Info::Inst()->write("Staging finishe...");	
    
    
    // start computation tasks

    world.barrier();

    if (world.rank()==0) print_analysis();
    
	if (world.rank()==0) {
		Info::Inst()->write(string("Aggregating timing information for performance analysis..."));		
	}
	
	timer.stop("total");
	
    std::map<boost::thread::id,Timer> performance_timer;
    performance_timer[boost::this_thread::get_id()]=timer;
    PerformanceAnalyzer perfanal(world,performance_timer); // collect timing information from everybody.
        
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
