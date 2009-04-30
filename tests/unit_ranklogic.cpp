#include <iostream>
#include <vector>
#include <map>
#include <deque>

#include <boost/mpi.hpp>

#include "tasks.hpp"
#include "coor3d.hpp"
#include "scatterdata.hpp"

using namespace std;
namespace mpi = boost::mpi;
using namespace boost;

void print_tasks(Tasks tasks,int rank) {
	
		for(size_t i = 0; i < tasks.size(); ++i)
		{
			cout << tasks[i].id << "," << tasks[i].color << "|"<< tasks[i].timeunits() << "|" << "q=" << tasks[i].q << "\t" << "{";
			for(size_t j = 0; j < tasks[i].rftable.size(); ++j)
			{
				pair<int,int> rf = tasks[i].rftable[j];
				cout << rf.first << ",";
			}
			cout << "}" << "=>" << "{";
			for(size_t j = 0; j < tasks[i].rftable.size(); ++j)
			{
				pair<int,int> rf = tasks[i].rftable[j];
				cout << rf.second << ",";
			}	
			cout << "}" ;
			cout << endl;
		}
}

int main (int argc, char  **argv)
{
	//------------------------------------------//
	//
	// MPI Initialization
	//
	//------------------------------------------//	
	
  	boost::mpi::environment env(argc, argv);
  	boost::mpi::communicator world;
    // Get the number of processors this job is using:
	int rank = world.rank();

    // Get the rank of the processor this thread is running on.  (Each
    // processor has a unique rank.)
	int size = world.size();
	srand(rank);
	if (size<=0) { cout << "please run via mpirun" << endl; return -1; }
	
	if (argc!=4) {
		cout << "please supply #f and #q and mode (frame or none)" << endl;
		return -1;
	}
	size_t nf = atoi(argv[1]);
	size_t nq = atoi(argv[2]);	
	
	vector<int> frames;
	vector<CartesianCoor3D> qqq;
	
	for(size_t i = 0; i < nf; ++i)
	{
		frames.push_back(i);
	}
	for(size_t i = 0; i < nq; ++i)
	{
		qqq.push_back(CartesianCoor3D(1.0,1.0,1.0));
	}
	
	string mode(argv[3]);
	
	Tasks tasks(frames,qqq,size,mode);

	if (rank==0) {
		cout << "created the following tasks:" << endl;
//		print_tasks(tasks,rank);
		cout << "################################################################################" << endl;	
	}

	Tasks my_tasks = tasks.scope(rank);
	
	if (rank==1) {
		sleep(1);
		cout << "rank==1 has the following tasks:" << endl;
//		print_tasks(my_tasks,rank);
		cout << "################################################################################" << endl;
	}
	
	world.barrier();
	
	world.barrier();
	
	if (my_tasks.size()<1) return 0;
	
	// create a task specific communicator
	boost::mpi::communicator local = world.split(my_tasks[0].color);
	
	vector<pair<CartesianCoor3D,double> > scat;
	map<int,vector<pair<CartesianCoor3D,double> > > scatbyframe;
	
	// do 'computation'
	for(size_t i = 0; i < my_tasks.size(); ++i)
	{
		local.barrier(); // say hello to all others
		scatbyframe.clear(); // each round has new results
		
		// iterate through all 'my' frames
		vector<int> frames_i = my_tasks[i].frames(rank); // use global rank here
		for(size_t j = 0; j < frames_i.size(); ++j)
		{
			// compute specific frame:
			// unfold q vector (here: 2)
			for(size_t k = 0; k < 2; ++k)
			{
//				scatbyframe[frames_i[j]].push_back(make_pair(my_tasks[i].q,rand()*1.0/RAND_MAX));
				scatbyframe[frames_i[j]].push_back(make_pair(my_tasks[i].q,1.0));
			}
		}
		
		//stats before communication
		local.barrier();
		
		
		if (rank==0) {
//			cout << "STATS:" << endl << flush;
//			cout << "RANK"<< rank << ": scatbyframe=" << scatbyframe.size() << endl;
			for(map<int,vector<pair<CartesianCoor3D,double> > >::iterator jj = scatbyframe.begin(); jj != scatbyframe.end(); ++jj)
			{
//				cout << jj->first << ": ";
				for(size_t k = 0; k < jj->second.size(); ++k)
				{
//					cout << jj->second[k].second << "\t";
				}
//				cout << endl;
			}
//			cout << flush;
		}
		local.barrier();
		// 
		vector< map<int,vector<pair<CartesianCoor3D,double> > > > scatbyframe_vector_in;
		vector< map<int,vector<pair<CartesianCoor3D,double> > > > scatbyframe_vector_out;

		scatbyframe_vector_in.resize(local.size(),scatbyframe);
		scatbyframe_vector_out.resize(local.size());		
		
		boost::mpi::all_to_all(local,&scatbyframe_vector_in[0],&scatbyframe_vector_out[0]);
		
		// by now we know:
		// unfolded q = 2
		// local rank = index for q,A vector
		int my_q_number = local.rank();
		// only do something here if my_q_number < size of q unfolded.
		if (my_q_number<2) {
			// ok, we do FFT 
			vector<double> aqs; aqs.resize(frames.size()); //make this arrays as long as the trajectory is
			for(size_t i = 0; i < scatbyframe_vector_out.size(); ++i)
			{
				for(map<int,vector<pair<CartesianCoor3D,double> > >::iterator sbfi = scatbyframe_vector_out[i].begin(); sbfi != scatbyframe_vector_out[i].end(); sbfi++)
				{
					int thisframe = sbfi->first;
					pair<CartesianCoor3D,double>& this_qA = (sbfi->second)[my_q_number];
					aqs[thisframe] = this_qA.second;
				}
			}
			
			if (rank==1) {
//				cout << "AQS.size()=" << aqs.size() <<endl;
//				cout << "AQS=";
				for(size_t i = 0; i < aqs.size(); ++i)
				{
//					cout << aqs[i] << ",";
				}
//				cout << endl;
			}
		}
		
		local.barrier();
		
		
		
		if (rank==0) {
//			cout << "STATS:" << endl << flush;
//			cout << "RANK"<< rank << ": scatbyframe_vector_out=" << scatbyframe_vector_out.size() << endl;
			for(size_t i = 0; i < scatbyframe_vector_out.size(); ++i)
			{
//				cout << " map entries=" << scatbyframe_vector_out[i].size() << endl;
				for(map<int,vector<pair<CartesianCoor3D,double> > >::iterator sbfi = scatbyframe_vector_out[i].begin(); sbfi != scatbyframe_vector_out[i].end(); sbfi++) {

				}
			}
			
		}

	}
	local.barrier();
//	if (rank==0) cout << "STATS:" << endl << flush;
	local.barrier();

//	cout << "RANK"<< rank << ": frames=" << scatbyframe.size() << " vectors=" << scatbyframe[0].size() << endl;

	
//	boost::mpi::all_to_all(local,&my_values[0],&other_values[0]);

	
	world.barrier();
	
	return 0;
}
