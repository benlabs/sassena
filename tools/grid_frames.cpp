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
typedef pair<pair<CartesianCoor3D,size_t>,pair<double,double> > keyval_t;
typedef list<keyval_t> scatlist_t;
	
bool cmpkeyval(keyval_t left,keyval_t right) {
	// then cluster by vectors
	// x being first iso line
	if (left.first.first.x < right.first.first.x)
		return true;

	if (left.first.first.x > right.first.first.x)
		return false;
		
	if (left.first.first.y < right.first.first.y)
		return true;

	if (left.first.first.y > right.first.first.y)
		return false;

	if (left.first.first.z < right.first.first.z)
		return true;

	if (left.first.first.z > right.first.first.z)
		return false;

	if (left.first.second < right.first.second)
		return true;

	if (left.first.second > right.first.second)
		return false;


	return false;
}

int main(int argc,char** argv) {

	if (argc!=4) {
		cerr << "ERROR>> " << "Need three arguments: sort-field(numeric) scatter-file output-file" << endl;
		cerr << "ERROR>> " << "sort-order can be: 1 2 3 .." << endl;		
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
        float qx,qy,qz,fr,fsquared;
		ifile >> qx >> qy >> qz >> fr >> fsquared;
		if ( ifile.eof() ) break;			
					
		keyval.first = make_pair(CartesianCoor3D(qx,qy,qz),fr);
		keyval.second = make_pair(fsquared,0.0);
		scat.push_back(keyval);
	}
	
	clog << "INFO>> " << "Sorting values "<< endl;
	scat.sort(cmpkeyval);
	
	
	double lastx = scat.begin()->first.first.x;
	double lasty = scat.begin()->first.first.y;	
	double lastz = scat.begin()->first.first.z;	
	double lastfr = scat.begin()->first.second;
	clog << "INFO>> " << "Writing values to file"<< argv[3] <<  endl;
	int xcount =0;
	int ycount =0;	
	int zcount =0;		
	int frcount =0;	
	bool nl = false;
	
		
	int mod = atoi(argv[1]);
	
		cout << "sort field is: '" << mod << "'" << endl;
	
	for (scatlist_t::iterator si=scat.begin();si!=scat.end();si++) {

		xcount++;
		
		if (lastx!=si->first.first.x)
		{ 
			ycount++;
		}	

		if (lasty!=si->first.first.y)
		{ 
			zcount++;
		}	

		if (lastfr!=si->first.second)
		{ 
			frcount++;
		}			

		if (mod == 1) {
			if ( (lastx!=si->first.first.x) &&  nl)
			{ 
				ofile << endl;		
				nl=false;
			}
		}
		else if (mod==2) {
			if ( (lasty!=si->first.first.y) &&  nl)
			{ 
				ofile << endl;		
				nl=false;
			}	
		}
		else if (mod==3) {
			if ( (lastz!=si->first.first.z) &&  nl)
			{ 
				ofile << endl;		
				nl=false;
			}	
		}	
		else if (mod==4) {
			if ( (lastfr!=si->first.second) &&  nl)
			{ 
				ofile << endl;		
				nl=false;
			}	
		}			
		else {
			cerr << "ERROR>> " << "sort token not understood" << endl;
			throw;
		}
		
		lastx = si->first.first.x;
 		lasty = si->first.first.y;
 		lastz = si->first.first.z;
 		lastfr = si->first.second;
		if ((xcount % skip)!=0) continue;
		if ((ycount % skip)!=0) continue;
		if ((zcount % skip)!=0) continue;
		if ((frcount % skip)!=0) continue;
		
		nl = true;
		ofile << si->first.first.x << "\t";
		ofile << si->first.first.y << "\t";
		ofile << si->first.first.z << "\t";
		ofile << si->first.second << "\t";		
		ofile << si->second.first << "\t"; 
		ofile << si->second.second;
		ofile << endl;
			
	}
	
	return 0;
};