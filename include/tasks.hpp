/*
 *  task.hpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef TASKS_HPP_
#define TASKS_HPP_

// common header
#include "common.hpp"

// standard header
#include <iostream>
#include <map>
#include <string>
#include <vector>

// special library headers
#include <boost/serialization/access.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>
#include <boost/mpi.hpp>

// other headers
#include "coor3d.hpp"

class Task {
public:
	size_t id;
	CartesianCoor3D q;
	
	double result;
	
	boost::mpi::communicator comm;
	
	std::vector<int> ranks; // list of processes which are associated with this task
	std::vector<std::pair<int,int> > rftable; // maps a rank to a frame, a rank can have multiple frames
	
	std::string mode;
};

class Tasks : public std::vector<Task> {
	void generate_independent_tasks(std::vector<int>& frames,std::vector<CartesianCoor3D>& qqq,size_t nn,std::string mode);
	void generate_tasks_by_framecoupling(std::vector<int>& frames,std::vector<CartesianCoor3D>& qqq,size_t nn,std::string mode);


public:
	Tasks() {}
	Tasks(std::vector<int>& frames,std::vector<CartesianCoor3D>& qqq,size_t nn,std::string mode);

	void print();
	Tasks scope(int rank);
};

#endif
