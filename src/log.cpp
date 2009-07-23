#include "log.hpp"

using namespace std;


void Info::write(std::string word) {
	
	clog << "INFO>> ";
	clog << word;
	clog << endl;
	clog.flush();
}

void Warn::write(std::string word) {
	
	clog << "WARNING>> ";
	clog << word;
	clog << endl;
	clog.flush();
}

void Err::write(std::string word) {
	
	cerr << "ERROR>> ";
	cerr << word;
	cerr << endl;
	cerr.flush();
}

size_t Info::counter=0;
size_t Warn::counter=0;
size_t Err::counter=0;

