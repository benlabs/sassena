/*
 *  analysis.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

// direct header
#include "analysis.hpp"

// standard header
#include <complex>
#include <map>
#include <string>
#include <vector>

// special library headers
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/vector.hpp>


#include <boost/random/mersenne_twister.hpp>	
#include <boost/random/uniform_on_sphere.hpp>	
#include <boost/random/variate_generator.hpp>	

// other headers
#include "atom.hpp"
#include "atomselection.hpp"
#include "coor3d.hpp"
#include "geometry.hpp"
#include "grid.hpp"
#include "sample.hpp"
#include "settings.hpp"


using namespace std;
using namespace boost::accumulators;
using namespace boost::math;

namespace Analysis {

void qvectors_unfold_sphere(std::string avvectors, CartesianCoor3D q, uint32_t qseed,double resolution, vector<CartesianCoor3D>& qvectors) {
	
	SphericalCoor3D qs(q);
	
	// we can't unfold a point:
	if (qs.r==0.0) return;
	
	const double M_2PI = 2*M_PI;	
	const double radincr = resolution*(M_2PI)/360;		
	const double arcincr = radincr*qs.r;

	int num=0;

	srand(qseed);

	if (avvectors=="rasterlinear") {
		for (double t=0;t<M_PI;t+=radincr) {	
			for (double p=0;p<M_2PI;p+=radincr) {
				qvectors.push_back( SphericalCoor3D(qs.r,qs.phi+p,qs.theta+t) );
		    }
		}		
	}
	else if (avvectors=="rastersurface") {
		double thetal = M_PI*qs.r; // this is the max distance we "walk on the sphere in theta direction"
		for (double t=0;((t==0) || (t<thetal));t+=arcincr) {	
			double dtheta = (t/thetal)*M_PI;
			double phil = 2*M_PI*sin(dtheta)*qs.r; // this is the max distance we "walk on the sphere in phi direction"		
			if (phil>0) {
				for (double p=0;p<phil;p+=arcincr) {
					double dphi = (p/phil)*2*M_PI;
					qvectors.push_back( SphericalCoor3D(qs.r,qs.phi+dphi,qs.theta+dtheta) );
			    }
			}
			else {
				qvectors.push_back( SphericalCoor3D(qs.r,qs.phi,qs.theta+dtheta) );
			}
		}
	}
	else if (avvectors=="triangulation") {
		
		DrawSphereHelper d(CartesianCoor3D(0,0,0),qs.r,resolution*qs.r);
		d.draw();

		for(vector<CartesianCoor3D>::iterator vi=d.vectors.begin();vi!=d.vectors.end();vi++) {
			qvectors.push_back( SphericalCoor3D(qs.r,qs.phi+SphericalCoor3D(*vi).phi,qs.theta+SphericalCoor3D(*vi).theta) );
		}
	}
	else if (avvectors=="mcphiacos") {
		int num=0;
		while (num<resolution) {
			double p = 2*M_PI*(rand()*1.0/RAND_MAX);
			double t = acos(2*(rand()*1.0/RAND_MAX)-1);		
			qvectors.push_back( SphericalCoor3D(qs.r,qs.phi+p,qs.theta+t) );
			num++;
		}
	}
	else if (avvectors=="mcdoublesqrt") {
		double x1,x2;
		while (num<resolution) {
			x1 = (2.0*(rand()*1.0/RAND_MAX)) - 1.0;
			x2 = (2.0*(rand()*1.0/RAND_MAX)) - 1.0;
			double xl = powf(x1,2) + powf(x2,2);
			if ( xl >= 1.0 ) continue;

			double x = 2*x1* sqrt(1-xl);
			double y = 2*x2* sqrt(1-xl) / xl;
			double z = 1-2*xl;

			SphericalCoor3D s = CartesianCoor3D(x,y,z);
			qvectors.push_back( SphericalCoor3D(qs.r,qs.phi+s.phi,qs.theta+s.theta) );
			num++;
		}
	}
	else if (avvectors=="mcquaternion") {
		double x0,x1,x2,x3;
		
		while (num<resolution) {
			x0 = (2.0*(rand()*1.0/RAND_MAX)) - 1.0;
			x1 = (2.0*(rand()*1.0/RAND_MAX)) - 1.0;
			x2 = (2.0*(rand()*1.0/RAND_MAX)) - 1.0;
			x3 = (2.0*(rand()*1.0/RAND_MAX)) - 1.0;
			double xl = powf(x0,2) + powf(x1,2) + powf(x2,2) + powf(x3,2);
			if ( xl >= 1.0 ) continue;

			double x = 2* (x1*x3+x0*x2) / xl;
			double y = 2* (x2*x3-x0*x1) / xl;
			double z = 2* (powf(x0,2)+powf(x3,2)-powf(x1,2)+powf(x2,2)) / xl;

			SphericalCoor3D s = CartesianCoor3D(x,y,z);
			qvectors.push_back( SphericalCoor3D(qs.r,qs.phi+s.phi,qs.theta+s.theta) );
			num++;
		}
	}
	else if (avvectors=="mcboostunisphere") {
		
		boost::mt19937 rng; // that's my random number generator
		rng.seed(qseed);
		boost::uniform_on_sphere<double> s(3); // that's my distribution
		boost::variate_generator<boost::mt19937&, boost::uniform_on_sphere<double> > mysphere(rng,s);

		while (num<resolution) {

			vector<double> r = mysphere();

			double x = r[0];
			double y = r[1];
			double z = r[2];

			SphericalCoor3D s = CartesianCoor3D(x,y,z);
			qvectors.push_back( SphericalCoor3D(qs.r,qs.phi+s.phi,qs.theta+s.theta) );
			num++;
		}
	}


	
}

void qvectors_unfold_cylinder(std::string avvectors, CartesianCoor3D q, uint32_t qseed,double resolution, vector<CartesianCoor3D>& qvectors) {
	
	srand(qseed);
	
	const double M_2PI = 2*M_PI;
	const double radincr = (M_2PI*resolution)/360;			
	
	CylinderCoor3D qs(q);
	
	// we can't unfold a line:
	if (qs.r==0.0) return;
	
	// resolution is in degree, coordinates are in rad
	
	for (double phi=0;phi<M_2PI;phi+=radincr) {
		qvectors.push_back( CylinderCoor3D(qs.r,qs.phi+phi,qs.z) );	
	}
}

void set_scatteramp(Sample& sample,Atomselection as,CartesianCoor3D q,bool background) {
	
	Atoms& atoms = sample.atoms;
	double ql = q.length();
	map<string,double> sfquicklookup;
	map<pair<string,double>,double> esfquicklookup;

	string probetype =  (const char *) Settings::get("main")["scattering"]["probe"]["type"];
	string m = "constant";
	if (Settings::get("main")["scattering"]["probe"].exists("method")) {
		m = (const char *) Settings::get("main")["scattering"]["probe"]["method"];		
	}
				
	for (Atomselection::iterator asi=as.begin();asi!=as.end();asi++) {
		
		double sf=0;
		if (sfquicklookup.find(atoms[*asi].name)!=sfquicklookup.end()) {
			sf = sfquicklookup[atoms[*asi].name];	
		}
		else {
			if (m == "slater")
			{
				cerr<< "ERROR>> " << "slater xray not supported yet" << endl;
				throw;
			}
			else if (m == "table")
			{
				libconfig::Setting& coeff = Settings::get("scattering_factors")[probetype]["table"][atoms[*asi].name];
				double a1=coeff["a1"]; double b1=coeff["b1"];
				double a2=coeff["a2"]; double b2=coeff["b2"];
				double a3=coeff["a3"]; double b3=coeff["b3"];				
				double a4=coeff["a4"]; double b4=coeff["b4"];
				double c =coeff["c"];
				double arg2 = ql*ql;				
				sf = ( a1*exp(-b1*arg2)+a2*exp(-b2*arg2)+a3*exp(-b3*arg2)+a4*exp(-b4*arg2)+c ) ;
			}
			else {
				// ok , try to do a constant lookup by default
				try {
					sf = Settings::get("scattering_factors")[probetype][m][atoms[*asi].name];
				}
				catch(...) {
					cerr << "ERROR>> " << " scattering factor for atom '" << atoms[*asi].name << "' not found" << endl;
					throw;
				}
			}

			sfquicklookup[atoms[*asi].name]= sf;
		}			
		
		// calculate effective scattering length:
		if (background) {
			pair<string,double> key=make_pair(atoms[*asi].name,atoms[*asi].kappa);
			if (esfquicklookup.find(key)!=esfquicklookup.end()) {
				atoms[*asi].scatteramp = esfquicklookup[key];
			}
			else {
				double ev = atoms[*asi].volume;
				double k  = atoms[*asi].kappa;
				
				double esf = sf - sample.background*k*ev*exp(-1.0*powf(k*ev,2.0/3.0)*powf(ql,2)/(4*M_PI));
				atoms[*asi].scatteramp =  esf;
				esfquicklookup[key]=esf;				
			}
		}
		else {
			atoms[*asi].scatteramp = sf;
		}
	}
}


void scatter_none(Sample& sample,Atomselection as,CartesianCoor3D& q,std::vector<std::complex<double> >& aqs) {
	
	CoordinateSet& cs = sample.frames.current().coordinate_sets[as.name];

	for (size_t i=0;i<cs.x.size();i++) {
		double s = sample.atoms[cs.indexes[i]].scatteramp;

		CartesianCoor3D c(cs.x[i],cs.y[i],cs.z[i]);
		aqs.push_back( exp(-1.0*complex<double>(0,c*q)) * s );
	}
}

complex<double> scatter_none(Sample& sample,Atomselection as,CartesianCoor3D& q) {
	
	CoordinateSet& cs = sample.frames.current().coordinate_sets[as.name];
	
	complex<double> A = complex<double>(0,0);
	if (q.length()==0) {
		for (Atomselection::iterator asi=as.begin();asi!=as.end();asi++) {
			A += sample.atoms[*asi].scatteramp;
		}
	}
	
	else {


		for (size_t i=0;i<cs.x.size();i++) {
			double s = sample.atoms[cs.indexes[i]].scatteramp;

			CartesianCoor3D c(cs.x[i],cs.y[i],cs.z[i]);
			A += exp(-1.0*complex<double>(0,c*q)) * s;
		}
		
//	complex<double> A = complex<double>(0,0);		
//		const int N = as.size();
//		const int Nmod4 = N % 4;
//		const int Ndiv4 = N / 4;
//		
//	
//		for (int i=0;i<Ndiv4;i++) {
//			CartesianCoor3D c1 = sample.frames.current().coord3D(as[4*i  ]);
//			CartesianCoor3D c2 = sample.frames.current().coord3D(as[4*i+1]);
//			CartesianCoor3D c3 = sample.frames.current().coord3D(as[4*i+2]);
//			CartesianCoor3D c4 = sample.frames.current().coord3D(as[4*i+3]);
//				
//			A += exp(-1.0*complex<double>(0,c1*q)) * sample.atoms[as[4*i  ]].scatteramp;		
//			A += exp(-1.0*complex<double>(0,c2*q)) * sample.atoms[as[4*i+1]].scatteramp;		
//			A += exp(-1.0*complex<double>(0,c3*q)) * sample.atoms[as[4*i+2]].scatteramp;		
//			A += exp(-1.0*complex<double>(0,c4*q)) * sample.atoms[as[4*i+3]].scatteramp;		
//			
//		}
//		if (Nmod4==3) {
//			CartesianCoor3D c1 = sample.frames.current().coord3D(as[N - 3]);
//			CartesianCoor3D c2 = sample.frames.current().coord3D(as[N - 2]);
//			CartesianCoor3D c3 = sample.frames.current().coord3D(as[N - 1]);
//				
//			A += exp(-1.0*complex<double>(0,c1*q)) * sample.atoms[as[N - 3]].scatteramp;		
//			A += exp(-1.0*complex<double>(0,c2*q)) * sample.atoms[as[N - 2]].scatteramp;		
//			A += exp(-1.0*complex<double>(0,c3*q)) * sample.atoms[as[N - 1]].scatteramp;		
//		}
//		else if (Nmod4==2) {
//			CartesianCoor3D c1 = sample.frames.current().coord3D(as[N - 2]);
//			CartesianCoor3D c2 = sample.frames.current().coord3D(as[N - 1]);
//				
//			A += exp(-1.0*complex<double>(0,c1*q)) * sample.atoms[as[N - 2]].scatteramp;		
//			A += exp(-1.0*complex<double>(0,c2*q)) * sample.atoms[as[N - 1]].scatteramp;		
//		}
//		else if (Nmod4==1) {
//			CartesianCoor3D c1 = sample.frames.current().coord3D(as[N - 1  ]);
//				
//			A += exp(-1.0*complex<double>(0,c1*q)) * sample.atoms[as[N - 1  ]].scatteramp;		
//		}

//	for (Atomselection::iterator asi=as.begin();asi!=as.end();asi++) {
	}
	return A;

//	double qx = q.x;
//	double qy = q.y;
//	double qz = q.z;
//
//complex<double> A = complex<double>(0,0);	
//	vector<double> tmp;
//	tmp.resize(as.size(),0);
//
//	for (int i=0;i<as.size();i++) {
//		tmp[i] += sample.frames.current().x[as[i]] * qx;
//	}	
//	for (int i=0;i<as.size();i++) {
//		tmp[i] += sample.frames.current().y[as[i]] * qy;
//	}	
//	for (int i=0;i<as.size();i++) {
//		tmp[i] += sample.frames.current().z[as[i]] * qz;
//	}	
//
//	for (int i=0;i<as.size();i++) {
//		A += exp(-1.0*complex<double>(0,tmp[i])) * sample.atoms[as[i]].scatteramp;
//	}	
//	
//	return A;	


//	complex<double> A = complex<double>(0,0);
//	boost::numeric::ublas::vector<double> qv(3); qv(0) = q.x; qv(1) = q.y; qv(2) = q.z;
//	boost::numeric::ublas::vector<double> rvr(sample.frames.current().coord3Dmatrix.size1());
//	boost::numeric::ublas::vector<double> rvi(sample.frames.current().coord3Dmatrix.size1());	
//	rvr =  boost::numeric::ublas::prod(sample.frames.current().coord3Dmatrix,qv);
//	rvi = rvr;
//	
//	const int Na = sample.frames.current().coord3Dmatrix.size1();
//	
//	for (int j=0;j<Na;j++) {
//		rvr(j) = cos(rvr(j))*sample.atoms[as[j]].scatteramp;
//		rvi(j) = -1.0*sin(rvi(j))*sample.atoms[as[j]].scatteramp;	
//	}
//	int rvn = rvr.size();
//	double rvrs = 0;
//	double rvis = 0;
//	
//	for (int i=0;i<rvn;i++) {
//		rvrs += rvr(i);
//	}
//	for (int i=0;i<rvn;i++) {
//		rvis += rvi(i);
//	}
//	
//	return complex<double>(rvrs,rvis);
	

//	for (Atomselection::iterator asi=as.begin();asi!=as.end();asi++) {
//		CartesianCoor3D c = sample.frames.current().coord3D(*asi); 
//		A += exp(-1.0*complex<double>(0,c*q)) * sample.atoms[*asi].scatteramp;
//	}
	return A;
		
//	complex<double> A = complex<double>(0,0);
//	for (Atomselection::iterator asi=as.begin();asi!=as.end();asi++) {
//		CartesianCoor3D c = sample.frames.current().coord3D(*asi); 
//		A += exp(-1.0*complex<double>(0,c*q)) * sample.atoms[*asi].scatteramp;
//	}
//	return A;

//	double Aa = 0;		
//	double Ab = 0;
//	for (Atomselection::iterator asi=as.begin();asi!=as.end();asi++) {
//		CartesianCoor3D c = sample.frames.current().coord3D(*asi);
//		Aa += sample.atoms[*asi].scatteramp * cos(-1.0*(c*q));
//		Ab += sample.atoms[*asi].scatteramp * sin(-1.0*(c*q));
//	}	
//	return complex<double>(Aa,Ab);
}


void scatter_vectors (Sample& sample,Atomselection as,std::vector<CartesianCoor3D>& qvectors,std::vector<vector<std::complex<double> > >& scattering_amplitudes) {

	scattering_amplitudes.resize(qvectors.size());
	for(size_t i = 0; i < qvectors.size(); ++i)
	{
		vector<complex<double> >& scat = scattering_amplitudes[i];
		scatter_none(sample,as,qvectors[i],scat);
	}
}

void scatter_vectors (Sample& sample,Atomselection as,std::vector<CartesianCoor3D>& qvectors,std::vector<std::complex<double> >& scattering_amplitudes) {

	for(size_t i = 0; i < qvectors.size(); ++i)
	{
		scattering_amplitudes.push_back( scatter_none(sample,as,qvectors[i]) );
	}
}

void scatter_sphere_multipole (Sample& sample,Atomselection as,CartesianCoor3D q,double resolution,std::vector<std::complex<double> >& scattering_amplitudes) {

	using namespace boost::numeric::ublas::detail;
	
	const int lmax=int(resolution);
	
	std::vector<std::vector<complex<double> > > almv;

  	almv.resize(lmax+1);
  	for (int l=0;l<=lmax;l++) {
   		almv[l].resize(2*l+1,complex<double>(0,0));
  	}

  	for (Atomselection::iterator asi=as.begin();asi!=as.end();asi++) {
		SphericalCoor3D c = sample.frames.current().coord3D(*asi);
		double esf = sample.atoms[*asi].scatteramp;

		for (int l=0;l<=lmax;l++) {
			complex<double> fmpiilesf = 4.0*M_PI*pow(complex<double>(0,1.0),l) * esf;
			double aabess = sph_bessel(l,q.length()*c.r);
		
			for (int m=-l;m<=l;m++) {
		
			complex<double> aa = conj(spherical_harmonic(l,m,c.theta,c.phi)); 

			almv[l][m+l] += fmpiilesf * aabess* aa;
			}
		}
	
	}


	for (int l=0;l<=lmax;l++) {
		for (int m=-l;m<=l;m++) { 
			scattering_amplitudes.push_back( almv[l][m+l] / sqrt(4.0*M_PI) ); // maybe need to multiply w/ lmax^2
		}
	}
}

void scatter_cylinder_multipole(Sample& sample,Atomselection as,CartesianCoor3D q,double resolution,std::vector<std::complex<double> >& scattering_amplitudes) {

	Atoms& atoms = sample.atoms;
	
	CylinderCoor3D qs(q);

	const int lmax=resolution;
	
	vector<complex<double> > A,B,C,D;
	A.resize(lmax+1,complex<double>(0,0));
	B.resize(lmax+1,complex<double>(0,0));
	C.resize(lmax+1,complex<double>(0,0));
	D.resize(lmax+1,complex<double>(0,0));

	for (Atomselection::iterator asi=as.begin();asi!=as.end();asi++) {
		CylinderCoor3D c = sample.frames.current().coord3D(*asi);

		double esf = atoms[*asi].scatteramp;
		
		complex<double> expi = exp(complex<double>(0,c.z*qs.z));

//		A[0] += expi * (double)bessj(0,c.r*qs.r) *esf ;
		A[0] += expi * (double)cyl_bessel_j(0,c.r*qs.r) *esf ;

		for (int l=1;l<=lmax;l++) {
//			complex<double> fac1 = 2.0*powf(-1.0,l)*bessj(2*l,c.r*qs.r);
//			complex<double> fac2 = complex<double>(0,1.0)*double(2.0*powf(-1.0,l-1)*bessj(2*l-1,c.r*qs.r));

			complex<double> fac1 = 2.0*powf(-1.0,l)*cyl_bessel_j(2*l,c.r*qs.r);
			complex<double> fac2 = complex<double>(0,1.0)*double(2.0*powf(-1.0,l-1)*cyl_bessel_j(2*l-1,c.r*qs.r));

			A[l] += fac1*expi*cos(2*l*c.phi) *esf;
			B[l] += fac1*expi*sin(2*l*c.phi) *esf;
			C[l] += fac2*expi*cos((2*l-1)*c.phi) *esf;
			D[l] += fac2*expi*sin((2*l-1)*c.phi) *esf;
		}
	}

	scattering_amplitudes.push_back( A[0] ); // maybe need to multiply w/ lmax*4 + 1		
	for (int l=1;l<=lmax;l++) {
		scattering_amplitudes.push_back( sqrt(0.5)*A[l] ); // maybe need to multiply w/ lmax*4 + 1		
		scattering_amplitudes.push_back( sqrt(0.5)*B[l] ); // maybe need to multiply w/ lmax*4 + 1		
		scattering_amplitudes.push_back( sqrt(0.5)*C[l] ); // maybe need to multiply w/ lmax*4 + 1		
		scattering_amplitudes.push_back( sqrt(0.5)*D[l] ); // maybe need to multiply w/ lmax*4 + 1		
	}

}

complex<double> background(Sample& sample,int resolution, double hlayer, CartesianCoor3D q,map<string,vector<double> >& kappas) {
	
	vector<CartesianCoor3D> uc = sample.frames.current().unitcell;
	CartesianCoor3D origin = sample.frames.current().origin;
	
	string gi = Settings::get("main")["scattering"]["background"]["include"]; // name of background group
	string ge = Settings::get("main")["scattering"]["background"]["exclude"]; // name of particle group
		
	Atomselection as_system   = sample.atomselections["system"];
	Atomselection as_particle = sample.atomselections[ge];
	Atomselection as_solvent  = sample.atomselections[gi];
	
	vector<int> v; // an empty vector to initialize vGrid3D
	vGrid3D grid(resolution,uc,origin,v);
	vGrid3D solgrid(resolution,uc,origin,v);	
	vGrid3D partgrid(resolution,uc,origin,v);	
	Grid3D<bool> solventgrid(resolution,uc,origin,false);	
	Grid3D<bool> particlegrid(resolution,uc,origin,false);		
	Grid3D<bool> systemgrid(resolution,uc,origin,false);

	for(Atomselection::iterator asi=as_system.begin();asi!=as_system.end();asi++) {
		CartesianCoor3D c = sample.frames.current().coord3D(*asi);
		Gridkey3D gk = grid.get_cell(c);
		grid[gk].push_back(*asi);
	}	
		
	for(Atomselection::iterator asi=as_particle.begin();asi!=as_particle.end();asi++) {
		CartesianCoor3D c = sample.frames.current().coord3D(*asi);
		CartesianCoor3D hc(hlayer,hlayer,hlayer);
		vector<Gridkey3D> gks =particlegrid.get_cells(c,hc);
		
		for (vector<Gridkey3D>::iterator gki=gks.begin();gki!=gks.end();gki++) {
			particlegrid[*gki]=true;			
		}
	}	
	for(Atomselection::iterator asi=as_solvent.begin();asi!=as_solvent.end();asi++) {
		CartesianCoor3D c = sample.frames.current().coord3D(*asi);
		Gridkey3D gk = grid.get_cell(c);
		
		solventgrid[gk]=true;
	}
	
	// remove valid solvent grid cells...
	for(Grid3D<bool>::iterator gi=particlegrid.begin();gi!=particlegrid.end();gi++) {
		if (gi->second) {
			solventgrid[Gridkey3D(gi->first.ix,gi->first.iy,gi->first.iz)]=false;
		}
	}	
	// build an atomselection...
	Atomselection as_solvent_trunc(sample.atoms,false);
	int numsolatoms=0;
	double solventvolume=0;
	for(Grid3D<bool>::iterator gi=solventgrid.begin();gi!=solventgrid.end();gi++) {
		if (gi->second) {
			vector<int>& v = grid[Gridkey3D(gi->first.ix,gi->first.iy,gi->first.iz)];
			for (vector<int>::iterator vi=v.begin();vi!=v.end();vi++) {
				as_solvent_trunc.add(sample.atoms[*vi]);
				numsolatoms++;
			}
			solventvolume+= solventgrid.element_volume();
		}
	}	
	
	// calculate scattering amplitude...
	set_scatteramp(sample,as_solvent_trunc,q,false);		
//	complex<double> Ac = scatter(*(as_solvent_trunc.sample),as_solvent_trunc,q);

	complex<double> Ac = scatter_none(sample,as_solvent_trunc,q);	
//	cout << "scatter numatoms=" << numsolatoms << " , " << as_solvent_trunc.size() << endl;
//	cout << "scatter solvent_trunc=" << Ac << endl;
	
//	cout << "volume:" << totalvolume << endl;
//	cout << "scatter/volume:" << Ac/totalvolume << endl;
	
		
	// build a solvent exclusive grid..

	for(Atomselection::iterator asi=as_solvent.begin();asi!=as_solvent.end();asi++) {
		CartesianCoor3D c = sample.frames.current().coord3D(*asi);
		Gridkey3D gk = solgrid.get_cell(c);
		
		solgrid[gk].push_back(*asi);
	}
	

	// build a particle exclusive grid..
	for(Atomselection::iterator asi=as_particle.begin();asi!=as_particle.end();asi++) {
		CartesianCoor3D c = sample.frames.current().coord3D(*asi);
		Gridkey3D gk = partgrid.get_cell(c);
		
		partgrid[gk].push_back(*asi);
	}

	//scan solventgrid(bool) and grid(atoms) for any outer solvent molecules...
	double sum_shell_solvent=0;
	int    n_shell_solvent=0;
		int num_solcells=0;
		
		double dens = 0;
		
	for(Grid3D<bool>::iterator gi=solventgrid.begin();gi!=solventgrid.end();gi++) {
		if (gi->second) {
			num_solcells++;
			vector<int>& v = solgrid[Gridkey3D(gi->first.ix,gi->first.iy,gi->first.iz)];
			int ln =0;
			for (vector<int>::iterator vi=v.begin();vi!=v.end();vi++) {
				sum_shell_solvent+=sample.atoms[*vi].volume;
				n_shell_solvent++;
				ln++;
			}
			dens += 1e27*ln/(solventgrid.element_volume()) / 6.0221415e23 /3;
		}
	}	

	// scan particlegrid(bool) and solgrid(atoms) for inner solvent molecules...
	double hparticlevolume=0;
	double sum_all_solvent=0;
	double sum_inner_solvent=0;	
	for(Grid3D<bool>::iterator gi=particlegrid.begin();gi!=particlegrid.end();gi++) {
		vector<int>& v = solgrid[Gridkey3D(gi->first.ix,gi->first.iy,gi->first.iz)];		

		if (gi->second) {
			hparticlevolume+= solgrid.element_volume();
			for (vector<int>::iterator vi=v.begin();vi!=v.end();vi++) {
				sum_inner_solvent+=sample.atoms[*vi].volume;
			}
		}
		for (vector<int>::iterator vi=v.begin();vi!=v.end();vi++) {
			sum_all_solvent+=sample.atoms[*vi].volume;
		}		
	}	
	
	// scan particlegrid(bool) and partgrid(atoms) for particle molecules...
	double sum_particle = 0;
	for(Grid3D<bool>::iterator gi=particlegrid.begin();gi!=particlegrid.end();gi++) {
		if (gi->second) {
			vector<int>& v = partgrid[Gridkey3D(gi->first.ix,gi->first.iy,gi->first.iz)];
			for (vector<int>::iterator vi=v.begin();vi!=v.end();vi++) {
				sum_particle+=sample.atoms[*vi].volume;
			}
		}
	}	
	
	
	// kappa estimation:
	double totalvolume = grid.total_volume();
	double kappa_s = solventvolume/sum_shell_solvent;
	double particle_volume = totalvolume - sum_all_solvent * kappa_s;
//	double particle_volume = totalvolume - sum_shell_solvent * kappa_s;	
	double kappa_p = particle_volume / sum_particle;
//	double kappa_p = particle_volume / sum_particle;	

///	cout << resolution << "\t" << hlayer << "\t" << dens/num_solcells << endl;
	
//	cout << "number of solvent cells: " << num_solcells << endl;
//	cout << "n_shell_solvent: " << n_shell_solvent << endl;
//	cout << "numsolatoms: " << numsolatoms << endl;	
//	cout << "solventvolume: " << solventvolume << endl;
//	cout << "cell volume: " << xd*yd*zd << endl;
//	if (numsolatoms!=0)	cout << "solvent particles per cell: " << numsolatoms/num_solcells << endl;	
//	cout << "average density of bulk water (mol./L): " << 1e27*n_shell_solvent/solventvolume / 6.0221415e23 /3  << endl;

//	cout << "totalvolume=" << xd*yd*zd*resolution*resolution*resolution << endl;
//	cout << "solventvolume=" << solventvolume << endl;
//	cout << "hparticlevolume=" << hparticlevolume << endl;
//
//	cout << "sum_shell_solvent=" << sum_shell_solvent << endl;
//	cout << "sum_particle=" << sum_particle << endl;
//	cout << "particle_volume=" << particle_volume << endl;
//
//	cout << "sum_inner_solvent=" << sum_inner_solvent << endl;
//	cout << "kappa_s=" << kappa_s << endl;
//	cout << "kappa_p=" << kappa_p << endl;	
//	cout << resolution << "\t" << hlayer << "\t" << kappa_s <<  "\t" << kappa_p << endl;

	// currently we have only a 2 phase system (allowed)
	// return kappa values via argument
	kappas[gi].push_back(kappa_s);
	kappas[ge].push_back(kappa_p);
	
//   if (Settings::get("main")["scattering"]["background"].exists("volume")) {
//   	// store kappa values for atoms:
//   	for(Atomselection::iterator asi=as_solvent.begin();asi!=as_solvent.end();asi++) {
//   //		(*asi)->kappa = kappa_s;
//   		sample.atoms[*asi].kappa = 1.0;
//   	}
//   	for(Atomselection::iterator asi=as_particle.begin();asi!=as_particle.end();asi++) {
//   //		(*asi)->kappa = kappa_p;
//   		sample.atoms[*asi].kappa = 1.0;
//   	}		
//   }
	
	return Ac/solventvolume;
		
}

pair<double,double> background_avg(Sample& sample,int resolution,double hydration, CartesianCoor3D q) {
		
	accumulator_set<double,features<tag::mean,tag::variance> > acc_RE;
	accumulator_set<double,features<tag::mean,tag::variance> > acc_IM;

	map<string,vector<double> > kappas;
		
	int framestride = 1;
	if (Settings::get("main")["scattering"]["background"].exists("framestride")) {
		framestride = Settings::get("main")["scattering"]["background"]["framestride"];
	}
	
	int totalframes = (sample.frames.size()/framestride) + 1;
	if ((sample.frames.size() % framestride) ==0) totalframes--;
	clog << "INFO>> " << "Average background density over " << totalframes << " frames" << endl;

	clog << "INFO>> " << "Progress(%): ";
	int progress=0;
	for(size_t i=0;i<sample.frames.size();i+= framestride) {
		int percent = (i*100/sample.frames.size());
		if ( percent >= progress  ) {
			progress = 10*(percent/10) + 10;
			clog << percent << " ";
		}
				
		sample.frames.load(i,sample.atoms,sample.atomselections);
   	 	complex<double> b0v = Analysis::background(sample,resolution,hydration,q,kappas);
		acc_RE(b0v.real()); acc_IM(b0v.imag());
	}
	clog << "100" << endl;

	// write kappa values into the atom entries...
	// create temp. acc
	for (map<string,vector<double> >::iterator kmi=kappas.begin();kmi!=kappas.end();kmi++) {
		accumulator_set<double,features<tag::mean,tag::variance> > k;
		for (vector<double>::iterator ki=kmi->second.begin();ki!=kmi->second.end();ki++) {
			k(*ki);
		}
		double m = mean(k); double sv = sqrt(variance(k));
		clog << "INFO>> " << "atomic volume scaling factor for group " << kmi->first << ": " << m << " +- " << sv << endl;
		for (Atomselection::iterator asi=sample.atomselections[kmi->first].begin();asi!=sample.atomselections[kmi->first].end();asi++) {
			sample.atoms[*asi].kappa = m;
		}
	}

	return pair<double,double>(mean(acc_RE),sqrt(variance(acc_RE)));

}

} // end of namespace 

// some old code:


//CartesianCoor3D Analysis::center_of_mass(Sample& sample,int framenumber, string weight) {
//	double x,y,z,m,mi;
//	x = y = z = 0.0;
//	m = 0.0;
//
//	for (vector<Atom>::iterator atom_i=sample.atoms.begin();atom_i!=sample.atoms.end();atom_i++) {
//		if (weight=="mass") {
//			mi = atom_i->mass;
//		}
//		else if ((weight=="particle") || (weight=="none")) {
//			mi = 1.0;
//		}
//		else if (weight=="scatteramp") {
//			mi = atom_i->scatteramp;
//		}
//		else {
//			cerr << "Weighting scheme not defined" << endl; throw;
//		}
//
//		m += mi; x += atom_i->x*mi;	y += atom_i->y*mi; z += atom_i->z*mi;
//	}
//
//	return CartesianCoor3D(x/m,y/m,z/m);
//}

//CartesianCoor3D Analysis::center_of_geometry(Sample& sample) {
//	double xmax,ymax,zmax;
//	double xmin,ymin,zmin;
//
//	int counter =0;
//	bool first = true;
//	for (vector<Atom>::iterator atom_i=sample.atoms.begin();atom_i!=sample.atoms.end();atom_i++,counter++) {
//		if (first) {
//			xmin = xmax = sample.frames.current().x[counter];
//			ymin = ymax = sample.frames.current().y[counter];
//			zmin = zmax = sample.frames.current().z[counter];
//			first = false;
//		}
//		else {
//			if (xmin>sample.frames.current().x[counter]) xmin = sample.frames.current().x[counter];
//			if (ymin>sample.frames.current().y[counter]) ymin = sample.frames.current().y[counter];
//			if (zmin>sample.frames.current().z[counter]) zmin = sample.frames.current().z[counter];
//                     
//			if (xmax<sample.frames.current().x[counter]) xmax = sample.frames.current().x[counter];
//			if (ymax<sample.frames.current().y[counter]) ymax = sample.frames.current().y[counter];
//			if (zmax<sample.frames.current().z[counter]) zmax = sample.frames.current().z[counter];
//		}
//	}
//
//	return CartesianCoor3D((xmax+xmin)/2.0,(ymax+ymin)/2.0,(zmax+zmin)/2.0);
//}



//CartesianCoor3D Analysis::minimum_cell(Sample& sample) {
//	CartesianCoor3D min,tmp;
//
//	bool first = true;
//	for (int i=0;i<sample.dcdframes.number_of_frames();i++) {
//		// load currentframe
//		sample.read_frame(i);
//		if (sample.frames.current().block1.empty()) { cerr << "minimum_cell not working for unit-cell less dcds" << endl; throw;}
//		tmp.x = sample.frames.current().block1[0];
//		tmp.y = sqrt(sample.frames.current().block1[1]*sample.frames.current().block1[1]+ sample.frames.current().block1[2]*sample.frames.current().block1[2]);
//		tmp.z = sqrt(sample.frames.current().block1[3]*sample.frames.current().block1[3]+sample.frames.current().block1[4]*sample.frames.current().block1[4]+ sample.frames.current().block1[5]*sample.frames.current().block1[5]);
//
//		if (first) {
//			first = false;
//			min = tmp;
//		} else {
//			if (min.x>tmp.x) min.x = tmp.x;
//			if (min.y>tmp.y) min.y = tmp.y;
//			if (min.z>tmp.z) min.z = tmp.z;						
//		}
//	}
//
//	return min;
//}


//pair< cartrect,cartrect> Analysis::scan_borders(Sample& sample) {
//	CartesianCoor3D pmin,pmax; // particle borders
//	CartesianCoor3D smin,smax; // system borders
//
//	bool firstatom = true; bool firstparticle = true;
//	vector<Atom>::iterator ai = sample.atoms.begin();
//	for (size_t j=0;j<sample.frames.current().x.size();j++,ai++) {
//		if (ai->particle) {
//			if (firstparticle) {
//			pmin.x = pmax.x = sample.frames.current().x[j];
//			pmin.y = pmax.y = sample.frames.current().y[j];
//			pmin.z = pmax.z = sample.frames.current().z[j];
//			firstparticle = false;
//			}
//
//			if (pmin.x>sample.frames.current().x[j]) pmin.x = sample.frames.current().x[j];
//			if (pmin.y>sample.frames.current().y[j]) pmin.y = sample.frames.current().y[j];
//			if (pmin.z>sample.frames.current().z[j]) pmin.z = sample.frames.current().z[j];
//
//			if (pmax.x<sample.frames.current().x[j]) pmax.x = sample.frames.current().x[j];
//			if (pmax.y<sample.frames.current().y[j]) pmax.y = sample.frames.current().y[j];
//			if (pmax.z<sample.frames.current().z[j]) pmax.z = sample.frames.current().z[j];
//		}
//		if (firstatom) {
//		smin.x = smax.x = sample.frames.current().x[j];
//		smin.y = smax.y = sample.frames.current().y[j];
//		smin.z = smax.z = sample.frames.current().z[j];
//		firstatom = false;
//		}
//
//		if (smin.x>sample.frames.current().x[j]) smin.x = sample.frames.current().x[j];
//		if (smin.y>sample.frames.current().y[j]) smin.y = sample.frames.current().y[j];
//		if (smin.z>sample.frames.current().z[j]) smin.z = sample.frames.current().z[j];
//                                                      
//		if (smax.x<sample.frames.current().x[j]) smax.x = sample.frames.current().x[j];
//		if (smax.y<sample.frames.current().y[j]) smax.y = sample.frames.current().y[j];
//		if (smax.z<sample.frames.current().z[j]) smax.z = sample.frames.current().z[j];
//	}		
//	
//	return pair< cartrect,cartrect>(cartrect(smin,smax), cartrect(pmin,pmax));
//}

// broken
//pair<double,double> Analysis::kappa_sp(Sample& sample, int frame_number) {
//
//	// load currentframe
//	sample.read_frame(frame_number);
//	
//	// select shell
//	pair<double,double> rinrout = rinnerrouter(sample);
//
//	// loop over all atoms, test for: solvent? and in-shell?
//	double sum_shell_solvent=0;
//	double sum_inner_solvent=0;
//	double sum_particle=0;
//	double total=0;
//	int counter = 0;
//	for (vector<Atom>::iterator ai = sample.atoms.begin();ai!=sample.atoms.end();ai++,counter++) {
//		double xd = sample.frames.current().x[counter];
//		double yd = sample.frames.current().y[counter];
//		double zd = sample.frames.current().z[counter];	
//		double r = 0;
//		r = sqrt(xd*xd+yd*yd);
//		double ev = sqrt(powf(M_PI,3))*powf(Settings::get("excluded_volumes")[ai->name],3);
//		if (ai->solvent) {
//			if ( (r>rinrout.first) && (r<rinrout.second) )
//				sum_shell_solvent += ev;
//			if (r<rinrout.first)
//				sum_inner_solvent += ev;
//		}
//		if (ai->particle) {
//			sum_particle += ev;
//		}
//		total += ev;
//	}
//
//
//	double kappa_s, kappa_p;
//		// volume of the cylinder shell: V = PI * (router-rinner)^2 * z
//		double length;
//		if (sample.frames.current().has_unit_cell()) {
////			length = sample.frames.current().unit_cell().z;
//		}
//		else {
//			pair< cartrect,cartrect> borders = scan_borders(sample);
//			length = borders.first.second.z - borders.first.first.z;
//		}
//
//		kappa_s = (M_PI*(powf(rinrout.second,2)-powf(rinrout.first,2))*length)/sum_shell_solvent;
//		double particle_volume = M_PI*powf(rinrout.first,2)*length - sum_inner_solvent * kappa_s;
//		kappa_p = particle_volume / sum_particle;
//		
//	return pair<double,double>(kappa_s,kappa_p);
//
//}

//double Analysis::density(Atomselection& as_solvent,int resolution) {
//	Sample& sample = *(as_solvent.sample);
//	
//	vector<CartesianCoor3D> uc = sample.frames.current().unit_cell();
//	CartesianCoor3D origin = sample.frames.current().origin;
//	
//	vector<int> v; // an empty vector to initialize vGrid3D	
//	vGrid3D grid(resolution,uc,origin,v);
//		
//	for(Atomselection::iterator asi=as_solvent.begin();asi!=as_solvent.end();asi++) {
//		CartesianCoor3D c = as_solvent.sample->frames.current().coord3D((*asi)->index);
//		Gridkey3D gk = grid.get_cell(c);
//		grid[gk].push_back((*asi)->index);
//	}
//	
//	
//	//scan solventgrid(bool) and grid(atoms) for any outer solvent molecules...
//	double sum_shell_solvent=0;
//	int    n_shell_solvent=0;
//	int num_solcells=0;
//		
//	double dens = 0;
//		
//	ofstream densfile("dens.txt",ios_base::app);
//	for(vGrid3D::iterator gi=grid.begin();gi!=grid.end();gi++) {
//		num_solcells++;
//		vector<int>& v = grid[Gridkey3D(gi->first.ix,gi->first.iy,gi->first.iz)];
//		int ln =0;
//		for (vector<int>::iterator vi=v.begin();vi!=v.end();vi++) {
//			n_shell_solvent++;
//			ln++;
//		}
//		dens += 1e27*ln/(grid.element_volume()) / 6.0221415e23 /3;
//		densfile << 1e27*ln/(grid.element_volume()) / 6.0221415e23 /3 << endl;
//	}	
//
//	cout << resolution << "\t"  << dens/num_solcells << endl;
//	return dens/num_solcells;
//}

//pair<double,double> Analysis::densityavg(Sample& sample,Atomselection as_solvent,int resolution) {
////			accumulator_set<double,features<tag::mean,tag::variance> > acc_RE;
////			accumulator_set<double,features<tag::mean,tag::variance> > acc_IM;
//			
//	accumulator_set<double,features<tag::mean,tag::variance> > acc_RE;
//	for(int i=0;i<sample.dcdframes.number_of_frames();i++) {
////	for(int i=0;i<10;i++) {
////		sample.read_frame(i,as_particle);
//		sample.read_frame(i,as_solvent);		
//
//		double dens = Analysis::density(as_solvent,resolution);
////	}	
//		acc_RE(dens);
//	}
//
//	return pair<double,double>(mean(acc_RE),sqrt(variance(acc_RE)));
//
////	return make_pair(mean(acc_RE),sqrt(variance(acc_RE)));
//}




// end of file
