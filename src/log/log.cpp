#include "log/log.hpp"

#include <boost/format.hpp>

using namespace std;


void Info::write(std::string word) {
	boost::format linefmt("%|15|%|-65|");
	linefmt % prefix % word;
	clog << linefmt.str() << endl;
	clog.flush();
}

void Warn::write(std::string word) {
	boost::format linefmt("%|15|%|-65|");
	linefmt % prefix % word;
	clog << linefmt.str() << endl;
	clog.flush();
}

void Err::write(std::string word) {
	boost::format linefmt("%|15|%|-65|");
	linefmt % prefix % word;
	cerr << linefmt.str() << endl;
	cerr.flush();

}


size_t Info::counter=0;
size_t Warn::counter=0;
size_t Err::counter=0;

