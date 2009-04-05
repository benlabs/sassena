#include <iostream>
#include <libconfig.h++>

#include "settings.hpp"

// this is a declaration
std::map<std::string,libconfig::Setting&> settings;

using namespace std;

int main(int argc,char** argv) {
	
 	Settings::read(&argc,argv);
	if (argc !=0) { 
		cerr << "error during read out of the configuration file" << endl; exit(-1);
	}
	
	for (int i;i<setting["main"].getLength();i++) {
		cout << "Setting : " << setting[i].getName() << endl;
	}
	
	return 0;
}

// end of file
