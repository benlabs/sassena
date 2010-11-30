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
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "control/database.hpp"
#include "control/parameters.hpp"
#include "sample.hpp"

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
        
        template <class T> void broadcast_class(boost::mpi::communicator& comm,T& any, size_t root) {
            std::stringstream stream;
            if (comm.rank()==root) {
                boost::archive::text_oarchive ar(stream); 
                ar << any;
            }
        	mpi::wrapper::broadcast_stream(comm,stream,root);
            if (comm.rank()!=root) {
                boost::archive::text_iarchive ar(stream); 
                ar >> any;
            }
        }
    }
}

template void mpi::wrapper::broadcast_class<Sample>(boost::mpi::communicator& comm,Sample& any, size_t root);
template void mpi::wrapper::broadcast_class<Database>(boost::mpi::communicator& comm,Database& any, size_t root);
template void mpi::wrapper::broadcast_class<Params>(boost::mpi::communicator& comm,Params& any, size_t root);


// end of file
