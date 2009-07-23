#include <iostream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/mpi.hpp>
// special library headers
#include <boost/serialization/access.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>

using namespace std;

class A {
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {

    }
	///////////////////
	
};

class B {
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {

    }
	///////////////////
	
};

class C {
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & a;
		ar & b;
    }
	///////////////////
	
public:
	A a;
	B b;
	
	static C* Inst() { static C instance; return &instance; }	
	void test(int num) { std::cout << "I am running:"<< num << std::endl; }
};

int main (int argc, char** argv)
{
	
  	boost::mpi::environment env(argc, argv);
  	boost::mpi::communicator world;
    // Get the number of processors this job is using:
	int rank = world.rank();

    // Get the rank of the processor this thread is running on.  (Each
    // processor has a unique rank.)
	int size = world.size();
		
	C c;
	c.Inst()->test(rank);
	
	
	if (rank==0) {
		broadcast(world,c,0);		
	}
	else {
		world.recv(0,boost::mpi::any_tag, c);					
	}

	world.barrier();
	cout << rank << ": done" << endl;
	
	
	return 0;
}