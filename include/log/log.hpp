/** \file
This file contains a set of singleton classes, which provide formatted console outputs.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/

#ifndef LOG__LOG_HPP_
#define LOG__LOG_HPP_

// common header
#include "common.hpp"

// standard header
#include <iostream>
#include <map>
#include <string>

/** 
Singelton class for formatted output to stdout, which should be used for information type messages.
*/
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

/** 
Singelton class for formatted output to stdout, which should be used for warning type messages.
*/
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

/** 
Singelton class for formatted output to stdout, which should be used for error type messages.
*/
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
