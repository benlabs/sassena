/*
 *  log.hpp
 *
 *  Created on: July 13, 2009
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef LOG__LOG_HPP_
#define LOG__LOG_HPP_

// common header
#include "common.hpp"

// standard header
#include <iostream>
#include <map>
#include <string>

// info, warn and err classes wrap the stdio capabilities. they are singleton
class Info {
private:
	Info() {}
	Info(const Info&) {}
    Info& operator=(const Info&) { return *this; }
	
	static size_t counter;	
	
	std::string prefix;
public:

	static Info* Inst() { static Info instance; counter++; return &instance; }

	void set_prefix(std::string p) { prefix = p;}

	void write(std::string);
	size_t calls() { return counter; }
	~Info() {}
};

class Warn {
private:
	Warn() {}
	Warn(const Warn&) {}
    Warn& operator=(const Warn&) { return *this; }
	
	static size_t counter;
	std::string prefix;

public:

	static Warn* Inst() { static Warn instance; counter++; return &instance; }

	void set_prefix(std::string p) { prefix = p;}

	void write(std::string);
	size_t calls() { return counter; }
	~Warn() {}
};

class Err {
private:
	Err() {}
	Err(const Err&) {}
	Err& operator=(const Err&) { return *this; }
		
	static size_t counter;		

	std::string prefix;
public:

	static Err* Inst() { static Err instance; counter++; return &instance; }

	void set_prefix(std::string p) { prefix = p;}

	void write(std::string);
	size_t calls() { return counter; }
	~Err() {}	
};

#endif 

// end of file
