/*
 *  common.hpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef COMMON_HPP_
#define COMMON_HPP_

// standard header
#include <map>

// special library headers
#include <libconfig.h++>

// other headers


// this is a declaration
extern std::map<std::string,libconfig::Setting*> settings;

#endif
