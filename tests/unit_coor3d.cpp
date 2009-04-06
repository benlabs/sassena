#include "coor3d.hpp"

#include <iostream>
#include <fstream>
#include <boost/math/special_functions.hpp>

#include "atoms.hpp"

using namespace std;
using namespace boost::math;
using namespace boost;
		
int main(int argc, char** argv) {


	ofstream ofile3(argv[1]);

	// let's try a Markov Chain!
	int num=0;
	double theta = 0;
	double phi = 0;	
	double theta2 = 0;
	double phi2 = 0;	
	
	double r =1.0;
	while (num<10000) {
		num++;
		if (num<9000) continue;
		double w = (1.0*rand()/RAND_MAX);
		double t = 2*M_PI*(1.0*rand()/RAND_MAX);
		double p = 2*M_PI*(1.0*rand()/RAND_MAX);
				
		if (w<0.5) {
			theta += t;
			phi += p;
		}
		else {
			phi += p;
			theta += t;
		}
						
		double w2 = (1.0*rand()/RAND_MAX);
		double t2 = 2*M_PI*(1.0*rand()/RAND_MAX);
		double p2 = 2*M_PI*(1.0*rand()/RAND_MAX);
		
		if (w2<0.5) {
			theta2 += t2;
			phi2 += p2;
		}
		else {
			phi2 += p2;
			theta2 += t2;
		}
		
		SphericalCoor3D o(0,phi2,theta2);
		CartesianCoor3D c(SphericalCoor3D(r,phi,theta));
		ofile3 << c.x << "\t" << c.y << "\t" << c.z << endl;
	}
	
	ofile3.close();

	return 0;

	ofstream ofile(argv[1]);
	
//	for (int k=1;k<=10;k++) {
		for (double t=-1;t<=1;t+=0.01) {
			for (double p=-1; p<=1;p+=0.01) {
				SphericalCoor3D s(1.0,0.5*M_PI+t*M_PI,0.5*M_PI+p*2*M_PI);
				CartesianCoor3D c = s;
				ofile << c.x << "\t"  << c.y << "\t" << c.z << endl;
//				ofile << 1.0 << "\t " << l << "\t "  << m << endl;
			}
		}
//	}
	ofile.close();
	
	//read a file in x,y,z format:
	ifstream infile(argv[1]);
	
	vector<CartesianCoor3D> testobject;
	
//	string line;
//	while (!(infile.eof())) {
//		double x,y,z;
//		infile >> x >> y >> z;
//		testobject.push_back(CartesianCoor3D(x,y,z));
//	}

	ofstream ofile2(argv[2]);
	
//  for(vector<CartesianCoor3D>::iterator ci=testobject.begin();ci!=testobject.end();ci++) {
//   SphericalCoor3D s = SphericalCoor3D(*ci);
//  s = SphericalCoor3D(s.r,s.phi,s.theta);
//   CartesianCoor3D c3b = s;
//	 c3b = rotate(*ci,"z",-s.phi);
//	 c3b = rotate(c3b,"y",0.01*M_PI);
//	 c3b = rotate(c3b,"z",s.phi);
//   *ci = c3b; 
//   ofile2 << c3b.x << "\t "  << c3b.y << "\t "  << c3b.z << endl;
//
//
//  }



	
	
	
	
	vector<CartesianCoor3D> o;
	// generate a 'random object' in 3d
	for(int i=0;i<10000; i++) {
		double x = (1.0*rand()/RAND_MAX);
		double y = (1.0*rand()/RAND_MAX);
		double z = (1.0*rand()/RAND_MAX);				
		o.push_back(CartesianCoor3D(x,y,z));
	}

	try {
	cout << "Testing class CartesianCoor3D... " << endl;
	
	CartesianCoor3D cart(0.0,0.0,0.0);

	cout << " CartesianCoor3D->CylinderCoor3D->CartesianCoor3D" << endl;

	for (vector<CartesianCoor3D>::iterator oi=o.begin();oi!=o.end();oi++) {
				CartesianCoor3D c1(oi->x,oi->y,oi->z);
				CylinderCoor3D c2(c1);
				CartesianCoor3D c3(c2);			
				if (abs(c1.x -c3.x)>1000*std::numeric_limits<float>::epsilon()) {
					cerr << "c1.x-c3.x : " << c1.x-c3.x << endl;
					cerr << "c1=" << c1 << endl;
					cerr << "c2=" << c2 << endl;
					cerr << "c3=" << c3 << endl;
					throw("CartesianCoor3D");
				}
				if (abs(c1.y-c3.y)>1000*std::numeric_limits<float>::epsilon()) {
					cerr << "c1.y-c3.y : " << c1.y-c3.y << endl;
					cerr << "c1=" << c1 << endl;
					cerr << "c2=" << c2 << endl;
					cerr << "c3=" << c3 << endl;
					throw("CartesianCoor3D");
				}
				if (abs(c1.z- c3.z)>1000*std::numeric_limits<float>::epsilon()) {
					cerr << "c1.z-c3.z : " << c1.z - c3.z << endl;
					cerr << "c1=" << c1 << endl;
					cerr << "c2=" << c2 << endl;
					cerr << "c3=" << c3 << endl;					
					throw("CartesianCoor3D");
				}					
	}

	cout << " CartesianCoor3D (conserve, rotate)" << endl;

	// Test for rotate and 'conserve'
	for (vector<CartesianCoor3D>::iterator oi=o.begin();oi!=o.end();oi++) {
		CartesianCoor3D c1(oi->x,oi->y,oi->z);
		SphericalCoor3D c2(c1);
		c2.phi = c2.phi + 0.1*M_PI;
		c2.theta = c2.theta + 0.1*M_PI;		
		CartesianCoor3D c3(c2);
		SphericalCoor3D c4(c3);
		c4.phi = c4.phi - 0.1*M_PI;
		c4.theta = c4.theta - 0.1*M_PI;		
		CartesianCoor3D c5(c4);
		if (abs(c1.x -c5.x)>1000*std::numeric_limits<float>::epsilon()) {
			cerr << "c1.x-c3.x : " << c1.x-c3.x << endl;
			cerr << "c1=" << c1 << endl;
			cerr << "c2=" << c2 << endl;
			cerr << "c3=" << c3 << endl;
			throw("CartesianCoor3D");
		}
		if (abs(c1.y-c5.y)>1000*std::numeric_limits<float>::epsilon()) {
			cerr << "c1.y-c3.y : " << c1.y-c3.y << endl;
			cerr << "c1=" << c1 << endl;
			cerr << "c2=" << c2 << endl;
			cerr << "c3=" << c3 << endl;
			throw("CartesianCoor3D");
		}
		if (abs(c1.z- c5.z)>1000*std::numeric_limits<float>::epsilon()) {
			cerr << "c1.z-c3.z : " << c1.z - c3.z << endl;
			cerr << "c1=" << c1 << endl;
			cerr << "c2=" << c2 << endl;
			cerr << "c3=" << c3 << endl;					
			throw("CartesianCoor3D");
		}
	}
	
	cout << " CartesianCoor3D (expl. rotate vs. rotate)" << endl;
	
	// Test for correct rotation properties
	// Test for rotate and 'conserve'
	for (vector<CartesianCoor3D>::iterator oi=o.begin();oi!=o.end();oi++) {
		CartesianCoor3D c1(oi->x,oi->y,oi->z);
		SphericalCoor3D c2(c1);
		c2.phi = c2.phi + 0.25*M_PI;
		c2.theta = c2.theta + 0.25*M_PI;		
		CartesianCoor3D c3(c2);
		CartesianCoor3D c1b(oi->x,oi->y,oi->z);
		SphericalCoor3D c2b(c1);
//		CartesianCoor3D c3b = rotate(rotate(rotate(c1b,"z",0.1*M_PI),"x",0.1*M_PI),"y",0.1*M_PI);
//		CartesianCoor3D c3b = rotate(rotate(rotate(c1b,"z",0.1*M_PI),"y",0.1*M_PI*cos(c2b.phi)),"x",0.1*M_PI*sin(c2b.phi));
		CartesianCoor3D c3b = rotate(c1b,"z",0.25*M_PI);
		 c3b = rotate(c3b,"z",-c2b.phi-0.25*M_PI);
		 c3b = rotate(c3b,"y",0.25*M_PI);
		 c3b = rotate(c3b,"z",c2b.phi+0.25*M_PI);	
		if (abs(c3.x -c3b.x)>1000*std::numeric_limits<float>::epsilon()) {
			cerr << "c3.x-c3b.x : " << c3b.x-c3.x << endl;
			cerr << "c1=" << c1 << endl;			
			cerr << "c3=" << c3 << endl;
			cerr << "c3b=" << c3b << endl;			
			throw("CartesianCoor3D");
		}
		if (abs(c3.y-c3b.y)>1000*std::numeric_limits<float>::epsilon()) {
			cerr << "c3.y-c3b.y : " << c3b.y-c3.y << endl;
			cerr << "c1=" << c1 << endl;			
			cerr << "c3=" << c3 << endl;
			cerr << "c3b=" << c3b << endl;
			throw("CartesianCoor3D");
		}
		if (abs(c3.z- c3b.z)>1000*std::numeric_limits<float>::epsilon()) {
			cerr << "c3.z-c3b.z : " << c3b.z - c3.z << endl;
			cerr << "c1=" << c1 << endl;			
			cerr << "c3=" << c3 << endl;
			cerr << "c3b=" << c3b << endl;
			throw("CartesianCoor3D");
		}
	}
	
	cout << " CartesianCoor3D->SphericalCoor3D->CartesianCoor3D" << endl;

	for (vector<CartesianCoor3D>::iterator oi=o.begin();oi!=o.end();oi++) {
				CartesianCoor3D c1(oi->x,oi->y,oi->z);
				SphericalCoor3D c2(c1);
				CartesianCoor3D c3(c2);			
				if (abs(c1.x -c3.x)>1000*std::numeric_limits<float>::epsilon()) {
					cerr << "c1.x-c3.x : " << c1.x-c3.x << endl;
					cerr << "c1=" << c1 << endl;
					cerr << "c2=" << c2 << endl;
					cerr << "c3=" << c3 << endl;					
					throw("CartesianCoor3D");
				}
				if (abs(c1.y-c3.y)>1000*std::numeric_limits<float>::epsilon()) {
					cerr << "c1.y-c3.y : " << c1.y-c3.y << ", " << c1.y - c3.y << endl;
					cerr << "c1=" << c1 << endl;
					cerr << "c2=" << c2 << endl;
					cerr << "c3=" << c3 << endl;					
					throw("CartesianCoor3D");
				}
				if (abs(c1.z- c3.z)>1000*std::numeric_limits<float>::epsilon()) {
					cerr << "c1.z-c3.z : " << c1.z - c3.z << endl;
					cerr << "c1=" << c1 << endl;
					cerr << "c2=" << c2 << endl;
					cerr << "c3=" << c3 << endl;					
					throw("CartesianCoor3D");
				}						
	}
	cout << " CartesianCoor3D->SphericalCoor3D(OUTOFBOUNDS)->CartesianCoor3D " << endl;
	for (vector<CartesianCoor3D>::iterator oi=o.begin();oi!=o.end();oi++) {
		CartesianCoor3D c1(oi->x,oi->y,oi->z);
			for (int p=-10;p<10;p++) {
				for (int t=-10;t<10;t++) {
						
				SphericalCoor3D c2(c1);
				c2.phi += p*2*M_PI;
				c2.theta += t*2*M_PI;
				
				CartesianCoor3D c3(c2);			
				if (abs(c1.x -c3.x)>1000*std::numeric_limits<float>::epsilon()) {
					cerr << "c1.x-c3.x : " << c1.x-c3.x << endl;
					cerr << "c1=" << c1 << endl;
					cerr << "c2=" << c2 << endl;
					cerr << "c3=" << c3 << endl;					
					throw("CartesianCoor3D");
				}
				if (abs(c1.y-c3.y)>1000*std::numeric_limits<float>::epsilon()) {
					cerr << "c1.y-c3.y : " << c1.y-c3.y << ", " << c1.y - c3.y << endl;
					cerr << "c1=" << c1 << endl;
					cerr << "c2=" << c2 << endl;
					cerr << "c3=" << c3 << endl;					
					throw("CartesianCoor3D");
				}
				if (abs(c1.z- c3.z)>1000*std::numeric_limits<float>::epsilon()) {
					cerr << "c1.z-c3.z : " << c1.z - c3.z << endl;
					cerr << "c1=" << c1 << endl;
					cerr << "c2=" << c2 << endl;
					cerr << "c3=" << c3 << endl;					
					throw("CartesianCoor3D");
				}	
				
				} }					
	}		
	
	cout << "Testing class CylinderCoor3D... " << endl;

	CylinderCoor3D cyl(0.0,0.0,0.0);

	
	for (int k=1;k<20;k++) {
		for (double l=0;l<2*M_PI;l+=M_PI/20) {
			for (int m=-10; m<10;m++) {
				CylinderCoor3D c1(k,l,m);
				CartesianCoor3D c2(c1);
				CylinderCoor3D c3(c2);			
				if (abs(c1.r -c3.r)>20*std::numeric_limits<double>::epsilon()) {
					cerr << "c1.r-c3.r : " << c1.r-c3.r << endl;
					cerr << "c1=" << c1 << endl;
					cerr << "c2=" << c2 << endl;
					cerr << "c3=" << c3 << endl;					
					throw("CylinderCoor3D");
				}
				if (abs(c1.phi-c3.phi)>20*std::numeric_limits<double>::epsilon()) {
					cerr << "c1.phi-c3.phi : " << c1.phi-c3.phi << endl;
					cerr << "c1=" << c1 << endl;
					cerr << "c2=" << c2 << endl;
					cerr << "c3=" << c3 << endl;					
					throw("CylinderCoor3D");
				}
				if (abs(c1.z- c3.z)>20*std::numeric_limits<double>::epsilon()) {
					cerr << "c1.z-c3.z : " << c1.z - c3.z << endl;
					throw("CylinderCoor3D");
				}					
			}
		}
	}

	ofstream tf("test");

	int k=10;
		for (double l=-4*M_PI;l<4*M_PI;l+=(4.0*M_PI/100)) {
			for (double m=-2*M_PI; m<2*M_PI;m+=(2.0*M_PI/100)) {
				SphericalCoor3D c1(k,l,m);
				CartesianCoor3D c2(c1);
				tf << c2.x << "\t" << c2.y << "\t" << c2.z << endl;
				SphericalCoor3D c3(c2);
			}
		}
	
	cout << "Testing class SphericalCoor3D... " << endl;
	
	
	
	SphericalCoor3D spher(0.0,0.0,0.0);
	
	for (int k=1;k<20;k++) {
		for (double l=0;l<2*M_PI;l+=M_PI/20) {
			for (int m=0; m<M_PI;m++) {
				if (m==0 && l!=0) continue;
				SphericalCoor3D c1(k,l,m);
				CartesianCoor3D c2(c1);
				SphericalCoor3D c3(c2);			
				if (abs(c1.r -c3.r)>20*std::numeric_limits<double>::epsilon()) {
					cerr << "c1.r-c3.r : " << c1.r-c3.r << endl;
					cerr << "c1=" << c1 << endl;
					cerr << "c2=" << c2 << endl;
					cerr << "c3=" << c3 << endl;					
					throw("SphericalCoor3D");
				}
				if (abs(c1.phi-c3.phi)>20*std::numeric_limits<double>::epsilon()) {
					cerr << "c1.phi-c3.phi : " << c1.phi-c3.phi << endl;
					cerr << "c1=" << c1 << endl;
					cerr << "c2=" << c2 << endl;
					cerr << "c3=" << c3 << endl;					
					throw("SphericalCoor3D");
				}
				if (abs(c1.theta- c3.theta)>20*std::numeric_limits<double>::epsilon()) {
					cerr << "c1.theta-c3.theta : " << c1.theta - c3.theta << endl;
					throw("SphericalCoor3D");
				}					
			}
		}
	}
	
	}
	catch(const char* e) {
		cout << "Error in Unit Test Coord3D: " << e << endl;
		return 1;
	}
	return 0;
}

// end of file
