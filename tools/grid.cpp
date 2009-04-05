#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <list>
#include <vector>

#include <coor3d.hpp>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//
// This routine takes the text output of a scatter file
// and sorts it acoording to isolines 
//
////////////////////////////////////////////////////////////////////////////////
typedef pair<CartesianCoor3D,pair<double,double> > keyval_t;
typedef list<keyval_t> scatlist_t;
	
bool cmpkeyval(keyval_t left,keyval_t right) {
	// then cluster by vectors
	// x being first iso line
	if (left.first.x < right.first.x)
		return true;

	if (left.first.x > right.first.x)
		return false;
		
	if (left.first.y < right.first.y)
		return true;

	if (left.first.y > right.first.y)
		return false;

	if (left.first.z < right.first.z)
		return true;

	if (left.first.z > right.first.z)
		return false;


	return false;
}

int main(int argc,char** argv) {

	if (argc!=4) {
		cerr << "ERROR>> " << "Need three arguments: sort-order scatter-file output-file" << endl;
		cerr << "ERROR>> " << "sort-order can be: x y z f" << endl;		
		return -1;
	}

	ifstream ifile(argv[2]);
	ofstream ofile(argv[3]);

	int skip =1;

	scatlist_t scat;


	clog << "INFO>> " << "Reading values from file " << argv[2] << endl;
	string line;
	while (true) {
		keyval_t keyval;
        float qx,qy,qz,fsquared,fsquared_dev;
		ifile >> qx >> qy >> qz >> fsquared >> fsquared_dev;
		if ( ifile.eof() ) break;			
					
		keyval.first = CartesianCoor3D(qx,qy,qz);
		keyval.second = make_pair(fsquared,fsquared_dev);
		scat.push_back(keyval);
	}
	
	clog << "INFO>> " << "Sorting values "<< endl;
	scat.sort(cmpkeyval);
	
	
	double lastx = scat.begin()->first.x;
	double lasty = scat.begin()->first.y;	
	double lastz = scat.begin()->first.z;		
	clog << "INFO>> " << "Writing values to file"<< argv[3] <<  endl;
	int xcount =0;
	int ycount =0;	
	int zcount =0;		
	bool nl = false;
	
		
		int mod = atoi(argv[1]);
	
		cout << "sort field is: '" << mod << "'" << endl;
	
	for (scatlist_t::iterator si=scat.begin();si!=scat.end();si++) {

		xcount++;
		
		if (lastx!=si->first.x)
		{ 
			ycount++;
		}	

		if (lasty!=si->first.y)
		{ 
			zcount++;
		}		

		if (mod == 1) {
			if ( (lastx!=si->first.x) &&  nl)
			{ 
				ofile << endl;		
				nl=false;
			}
		}
		else if (mod==2) {
			if ( (lasty!=si->first.y) &&  nl)
			{ 
				ofile << endl;		
				nl=false;
			}	
		}
		else if (mod==3) {
			if ( (lastz!=si->first.z) &&  nl)
			{ 
				ofile << endl;		
				nl=false;
			}	
		}
		else {
			cerr << "ERROR>> " << "sort token not understood" << endl;
			throw;
		}
		
		lastx = si->first.x;
 		lasty = si->first.y;
 		lastz = si->first.z;

		if ((xcount % skip)!=0) continue;
		if ((ycount % skip)!=0) continue;
		if ((zcount % skip)!=0) continue;
		
		nl = true;
		ofile << si->first.x << "\t";
		ofile << si->first.y << "\t";
		ofile << si->first.z << "\t";
		ofile << si->second.first << "\t"; 
		ofile << si->second.second;
		ofile << endl;
			
	}
	
	return 0;
};