#include <iostream>
#include <vector>
#include <map>
#include <deque>

#include <boost/mpi.hpp>

#include "tasks.hpp"
#include "coor3d.hpp"

using namespace std;
namespace mpi = boost::mpi;
using namespace boost;

void print_tasks(Tasks tasks,int rank) {
	
		for(size_t i = 0; i < tasks.size(); ++i)
		{
			cout << tasks[i].id << "||" << "q=" << tasks[i].q << "\t" << "{";
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
	
	Tasks tasks(frames,qqq,size,"none");

	if (rank==0) {
		cout << "created the following tasks:" << endl;
		print_tasks(tasks,rank);
		cout << "################################################################################" << endl;	
	}

	Tasks my_tasks = tasks.scope(rank);
	
	if (rank==1) {
		sleep(1);
		cout << "rank==1 has the following tasks:" << endl;
		print_tasks(my_tasks,rank);
		cout << "################################################################################" << endl;
	}
	
	world.barrier();
	
	// now each process has a list of my_tasks
	for(size_t i = 0; i < my_tasks.size(); ++i)
	{		
		Task& thistask = my_tasks[i];
//		boost::mpi::communicator comm = world.split(0);
		
//		cout << rank << ": I am " << comm.rank() << " of " << comm.size() << " in context " << thistask.id << endl; 
//		comm.barrier();
		
//		boost::mpi::communicator local = world.split(thistask.id);
//		local.barrier();
		
		my_tasks[i].result = 1.0;
	}

	double result = 0.0;
	for(size_t i = 0; i < my_tasks.size(); ++i)
	{		
		result += my_tasks[i].result;
	}

	
//	cout << "result(rank==" << rank << ")=" << result << endl;
	
	world.barrier();
	
	return 0;
}
