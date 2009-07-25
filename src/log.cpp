#include "log.hpp"

using namespace std;


void Info::write(std::string word) {
	
	clog << prefix;
	clog << word;
	clog << endl;
	clog.flush();
}

void Warn::write(std::string word) {
	
	clog << prefix;
	clog << word;
	clog << endl;
	clog.flush();
}

void Err::write(std::string word) {
	
	cerr << prefix;
	cerr << word;
	cerr << endl;
	cerr.flush();
}

size_t Info::counter=0;
size_t Warn::counter=0;
size_t Err::counter=0;

