/*
 *  This file is part of the software sassena
 *
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008-2011 Benjamin Lindner
 *
 */

#ifndef EXCEPTIONS__EXCEPTIONS_HPP_
#define EXCEPTIONS__EXCEPTIONS_HPP_

// common header
#include "common.hpp"

// standard header
#include <exception>

// special library headers

// other headers

namespace sassena {
    class terminate_request : public std::runtime_error {
    public:
        terminate_request() : std::runtime_error("Terminate requested") {}
    };
}


#endif

// end of file