/** \file
This file contains an Interface class for communicating data structures across the network.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/

// direct header
#include <sstream>

#include <boost/mpi.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include "control/database.hpp"
#include "control/parameters.hpp"
#include "sample.hpp"

using namespace std;

namespace mpi {
    namespace wrapper {
        void broadcast_stream(boost::mpi::communicator& comm,std::stringstream& stream, size_t root) {

            char* buffer = NULL;
            size_t buffersize = 0;

            if (comm.rank()==static_cast<long>(root)) {
                stream.seekg(0,ios_base::end);
                buffersize = stream.tellg();
            }
            boost::mpi::broadcast(comm,&buffersize,1,root);

            buffer = (char*) malloc(buffersize*sizeof(char));

            if (comm.rank()==static_cast<long>(root)) {
                stream.seekg(0,ios_base::beg);                                            
                stream.read(buffer,buffersize);
            }
            boost::mpi::broadcast(comm,buffer,buffersize,root);
            if (comm.rank()!=static_cast<long>(root)) {
                stream.write(buffer,buffersize);
            }
            free(buffer);

        }
        
        template <class T> void broadcast_class(boost::mpi::communicator& comm,T& any, size_t root) {
            std::stringstream stream(stringstream::in|stringstream::out|stringstream::binary);
            if (comm.rank()==static_cast<long>(root)) {
                boost::archive::binary_oarchive ar(stream); 
                ar << any;
            }
        	mpi::wrapper::broadcast_stream(comm,stream,root);
            if (comm.rank()!=static_cast<long>(root)) {
                boost::archive::binary_iarchive ar(stream); 
                ar >> any;
            }
        }
    }
}

template void mpi::wrapper::broadcast_class<Sample>(boost::mpi::communicator& comm,Sample& any, size_t root);
template void mpi::wrapper::broadcast_class<Database>(boost::mpi::communicator& comm,Database& any, size_t root);
template void mpi::wrapper::broadcast_class<Params>(boost::mpi::communicator& comm,Params& any, size_t root);
template void mpi::wrapper::broadcast_class<std::vector<std::string> >(boost::mpi::communicator& comm,std::vector<std::string>& any, size_t root);
template void mpi::wrapper::broadcast_class<std::string>(boost::mpi::communicator& comm,std::string& any, size_t root);

// end of file
