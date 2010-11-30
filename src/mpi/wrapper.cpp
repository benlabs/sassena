/*
 *  This file is part of the software sassena
 *
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008-2010 Benjamin Lindner
 *
 */

// direct header
#include <sstream>

#include <boost/mpi.hpp>
using namespace std;

namespace mpi {
    namespace wrapper {
        void broadcast_stream(boost::mpi::communicator& comm,std::stringstream& stream, size_t root) {

            char* buffer = NULL;
            size_t buffersize = 0;

            if (comm.rank()==root) {
                buffer = const_cast<char*>(stream.str().c_str());
                buffersize = stream.str().size();
            }
            boost::mpi::broadcast(comm,&buffersize,1,root);

            if (comm.rank()!=root) {
                buffer = (char*) malloc(buffersize*sizeof(char));
            }
            boost::mpi::broadcast(comm,buffer,buffersize,root);

            if (comm.rank()!=root) {
                for(size_t i = 0; i < buffersize; ++i) stream << buffer[i];
                free(buffer);
            }
        }
        
        class<T> void broadcast_class(boost::mpi::communicator& comm,T& any, size_t root) {
            std::stringstream stream;
            if (world.rank()==0) {
                boost::archive::text_oarchive ar(stream); 
                ar << *any;
            }
        	mpi::wrapper::broadcast_stream(world,stream,root);
            if (world.rank()!=0) {
                boost::archive::text_iarchive ar(stream); 
                ar >> *any;
            }
        }
    }
}

// end of file
