/*
 *  This file is part of the software sassena
 *
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008-2010 Benjamin Lindner
 *
 */

#ifndef MPI__WRAPPER_HPP_
#define MPI__WRAPPER_HPP_

#include "boost/mpi.hpp"

namespace mpi {
    namespace wrapper {
        void broadcast_stream(boost::mpi::communicator& comm,std::stringstream& stream, size_t root);        
    }
}

#endif

// end of file