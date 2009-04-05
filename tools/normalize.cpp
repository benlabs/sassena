#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

int main(int argc, char** argv) {
	
	ifstream infile(argv[1]);
	
	double normval = -1.0;
	
	bool first = true;
	double q,i,ierr;
	
	string line,token;
	while (getline(infile,line)) { 
		istringstream ss(line);
		ss >> q >> i >> ierr;
		if (first) { normval = i; first = false; }
		cout << q << "\t" << i/normval<< "\t" << ierr/normval << endl;
	}
	return 0;
}

// end of file
