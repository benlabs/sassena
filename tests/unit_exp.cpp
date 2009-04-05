

#include <iostream>
#include <complex>
#include <cmath>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/math/special_functions.hpp>
#include <vector>

#include "coor3d.hpp"

using namespace std;
using namespace boost::math;
using namespace boost::numeric::ublas::detail;

double mysinc_multipole(double r);
double mysinc_bruteforce(double r);
complex<double> myexp(CartesianCoor3D qc, CartesianCoor3D rc);
		
int main(int argc, char** argv) {

//	for (int i=0;i<20;i++) {
//	for (double x=0;x<30;x+=0.1) {
//		double mm= mysinc_multipole(x);
//		double mb= mysinc_bruteforce(x);
//		cout << x << "\t"<< "\t" << mm << "\t" << mb << "\t" << mm-mb  << endl;
		
//		cout << x << "\t"<< spherical_harmonic(atoi(argv[1]),atoi(argv[2]),1.7,x).real() << "\t" << spherical_harmonic(atoi(argv[1]),atoi(argv[2]),1.7,x).imag() << endl;		
//		cout << x << "\t"<< sph_bessel(i,x) << endl;		
		
//	}
//	cout << endl;
//  }

//	for (double x=0;x<50;x+=0.1) {
//		double mm= mysinc_multipole(x);
//		double mb= mysinc_bruteforce(x);
//		cout << x << "\t"<< "\t" << mm << "\t" << mb << "\t" << mm-mb  << endl;
//	}
	std::vector<CartesianCoor3D> vecs;

	vecs.push_back(CartesianCoor3D(1.9002,-6.83481,10.3639));
	vecs.push_back(CartesianCoor3D(2.47889,-7.16065,11.6203));
	vecs.push_back(CartesianCoor3D(1.83684,-7.05636,12.4825));

	CartesianCoor3D q(0.0,0.0,1.0);
	CartesianCoor3D r(0.0,0.0,1.0);

	for(std::vector<CartesianCoor3D>::iterator vi=vecs.begin();vi!=vecs.end();vi++) {
		complex<double> expi = exp(complex<double>(0.0,q*(*vi)));
		complex<double> mexpi = myexp(q,*vi);
		cout << "\t" << mexpi.real() <<"\t" << mexpi.imag() << "\t" << expi.real()<<"\t"<<expi.imag() << endl;
	}
	
	return  0;
};

double mysinc_multipole(double r) {
	
	using namespace boost::numeric::ublas::detail;

	// go up to the order of 5?
	const int lmax=50;

	std::vector<std::vector<complex<double> > > almv;

  almv.resize(lmax);
  for (int l=0;l<lmax;l++) {
   almv[l].resize(2*l+1,complex<double>(0,0));
  } 

CartesianCoor3D q(0.0,0.0,1.0);

	for (int i=0;i<10;i++) {
		for (int l=0;l<lmax;l++) {
			complex<double> fmpiilesf = 4.0*M_PI*pow(complex<double>(0,1.0),l);
			double aabess = sph_bessel(l,q.length()*r);

			for (int m=-l;m<=l;m++) {

//				complex<double> aa = conj(spherical_harmonic(l,m,0,0)); 
//				complex<double> fac = exp(complex<double>(0.0,q.z*r);
//				almv[l][m+l] += fmpiilesf * aabess* aa;
				almv[l][m+l] += fmpiilesf * aabess ;				
			}
		}

	}
	
	double	almvq = 0;  // almvq = abs(conj(almv[0])*(almv[0])); loop(k): almvq += abs(conj(almv[k])*(almv[k]));
	for (int l=0;l<lmax;l++) {
		for (int m=-l;m<=l;m++) { 
			almvq += (conj(almv[l][m+l])*(almv[l][m+l])).real(); 
		}
	}

	return almvq / (4.0*M_PI);
}


double mysinc_bruteforce(double r) {
	
	CartesianCoor3D q(0.0,0.0,1.0);
	
	const double M_2PI = 2*M_PI;	
	const double radincr = (M_2PI*1)/360;	
	int num=0;
	SphericalCoor3D qs(q);
	complex<double> A=0;
	double phi,theta;
	for (phi=0;phi<M_2PI;phi+=radincr) {
		for (theta=0;theta<M_PI;theta+=radincr) {
			complex<double> Aqs(0,0);
			SphericalCoor3D c(qs.r,qs.phi+phi,qs.theta+theta);
			for (int i=0;i<10;i++) {
				Aqs += 1.0;
			}
			A += conj(Aqs)*Aqs;
			num++;
		}
	}

	return A.real() / num;
}

complex<double> myexp(CartesianCoor3D qc, CartesianCoor3D rc) {
	
	SphericalCoor3D q(qc);
	SphericalCoor3D r(rc);
	
	int lmax=50;
	std::vector<std::vector<complex<double> > > almv;
	
	  almv.resize(lmax);
	  for (int l=0;l<lmax;l++) {
	   almv[l].resize(2*l+1,complex<double>(0,0));
	  } 

	for (int l=0;l<lmax;l++) {
		complex<double> fmpiilesf = 4.0*M_PI*pow(complex<double>(0,1.0),l);
		double aabess = sph_bessel(l,q.r*r.r);

		for (int m=-l;m<=l;m++) {
			complex<double> aa = conj(spherical_harmonic(l,m,r.theta,r.phi)); 
			complex<double> bb = spherical_harmonic(l,m,q.theta,q.phi); 			
			almv[l][m+l] += fmpiilesf * aabess * aa * bb;				
		}
	}

	complex<double>	almvq = 0;  // almvq = abs(conj(almv[0])*(almv[0])); loop(k): almvq += abs(conj(almv[k])*(almv[k]));
	for (int l=0;l<lmax;l++) {
		for (int m=-l;m<=l;m++) { 
			almvq += almv[l][m+l];
		}
	}

	return almvq;
}
