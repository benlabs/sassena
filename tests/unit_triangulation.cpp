#include "coor3d.hpp"

#include <complex>
#include <cmath>
#include <iostream>
#include <fstream>
#include <geometry.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/math/special_functions.hpp>
#include <vector>

using namespace std;

using namespace boost::math;
//using namespace boost::numeric::ublas::detail; 
	
int main(int argc, char** argv) {

  double r = 1;
double resolution = 0.1;

  DrawSphereHelper d(CartesianCoor3D(0,0,0),r,resolution*r);
  d.draw();
  cout << "INFO>> VECTORS=" << d.vectors.size() << endl;

  CartesianCoor3D bt(0,0,0);
  CartesianCoor3D bl(0,0,0);
  CartesianCoor3D bs(0,0,0);
  CartesianCoor3D bmc(0,0,0);

	ofstream obt("test-bt.txt");

	for (std::vector<CartesianCoor3D>::iterator vi=d.vectors.begin();vi!=d.vectors.end();vi++) {
		bt = bt+ *vi;
		CartesianCoor3D c2 = *vi;
			obt << c2.x << "\t" << c2.y << "\t" << c2.z  << endl;
	}
	cout << bt << endl;

	const double M_2PI = 2*M_PI;	
	const double radincr = 2.0*(M_2PI)/360;	
	const double arcincr = radincr*r;

	ofstream obs("test-bs.txt");

	double thetal = (M_PI/2)*r; // this is the max distance we "walk on the sphere in theta direction"
	for (double t=-M_PI/2;((t==-M_PI/2) || (t<thetal));t+=arcincr) {	
		double dtheta = (t/thetal)*M_PI;
		double phil = 2*M_PI*sin(dtheta)*r; // this is the max distance we "walk on the sphere in phi direction"		
		if (phil>0) {
			for (double p=0;p<phil;p+=arcincr) {
				double dphi = (p/phil)*2*M_PI;
				SphericalCoor3D c(r,dphi,dtheta);
				bs = bs + c;
				CartesianCoor3D c2 = c;
				obs << c2.x << "\t" << c2.y << "\t" << c2.z  << endl;				
		    }
		}
		else {
			SphericalCoor3D c(r,0,dtheta);
			bs = bs + c;
			CartesianCoor3D c2 = c;
			obs << c2.x << "\t" << c2.y << "\t" << c2.z  << endl;			
		}
		obs << endl;
	}
	cout << bs << endl;

	for (double t=-M_PI/2;t<M_PI/2;t+=radincr) {	
		for (double p=0;p<M_2PI;p+=radincr) {
			SphericalCoor3D c(1,p,t);
			bl = bl +c;
	    }
	}
	cout << bl << endl;
	
	cout << CartesianCoor3D(SphericalCoor3D(1.0,0,0)) 		<< endl;
	cout << CartesianCoor3D(SphericalCoor3D(1.0,0,0.5*M_PI))<< endl;
	cout << CartesianCoor3D(SphericalCoor3D(1.0,0,M_PI))	<< endl;

	ofstream obl("test-bl.txt");
	srand(time(NULL));
		for (double t=0;t<M_PI;t+=radincr) {
	for (double p=0;p<M_2PI;p+=radincr) {
		SphericalCoor3D c(1,p,t);
		CartesianCoor3D c2 = c;
		obl << c2.x << "\t" << c2.y << "\t" << c2.z  << endl;
    }
	obl << endl;
	}
	cout << "mc:" << endl;
	ofstream omc("test-mc.txt");
	
	int num=0;
	while (num<200000) {
		double p = 2*M_PI*(rand()*1.0/RAND_MAX);
		double t = acos(2*(rand()*1.0/RAND_MAX)-1);
		SphericalCoor3D c(1,p,t);
		CartesianCoor3D c2 = c;
		omc << c2.x << "\t" << c2.y << "\t" << c2.z  << endl;
		num++;
		bmc = bmc + c;
    }
	cout << bmc << endl;



	// go up to the order of 5?
	const int lmax=60;
	
	std::vector<std::vector<complex<double> > > almv;

  	almv.resize(lmax+1);
  	for (int l=0;l<=lmax;l++) {
   		almv[l].resize(2*l+1,complex<double>(0,0));
  	}

	for (int l=0;l<=lmax;l++) {
		double esf = 1.0;
		complex<double> fmpiilesf = 4.0*M_PI*pow(complex<double>(0,1.0),l) * esf;
		double aabess = boost::math::sph_bessel(l,r*r);
		
		for (int m=-l;m<=l;m++) {
			complex<double> aa = conj(boost::math::spherical_harmonic(l,m,0.1*M_PI,0*M_PI)); 
			almv[l][m+l] += fmpiilesf * aabess* aa;
		}
	}

	double	almvq = 0;  // almvq = abs(conj(almv[0])*(almv[0])); loop(k): almvq += abs(conj(almv[k])*(almv[k]));
	for (int l=0;l<=lmax;l++) {
		for (int m=-l;m<=l;m++) { 
			almvq += (conj(almv[l][m+l])*(almv[l][m+l])).real(); 
		}
	}
	cout << "multipole:" << endl;
	cout << almvq/(4.0*M_PI) << endl;


  return 0;
}

// end of file
