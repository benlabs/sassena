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
#include "log.hpp"
#include "density_grid.hpp"
#include "parameters.hpp"
#include "database.hpp"
#include "sample.hpp"
#include "timer.hpp"

using namespace std;
using namespace boost::accumulators;
using namespace boost::math;

namespace Analysis {


void scatter_cylinder_multipole(Sample& sample,Atomselection as,CartesianCoor3D q,double resolution,std::vector<std::complex<double> >& scattering_amplitudes) {

	Atoms& atoms = sample.atoms;
	
	CartesianCoor3D o = Params::Inst()->scattering.average.orientation.axis;
	o = o / o.length();
	
	// get the part of the scattering vector perpenticular to the o- orientation
	CartesianCoor3D qparallel = (o*q)*o; 
	CartesianCoor3D qperpenticular = q - qparallel; // define this as phi=0 
	double qperpenticular_l = qperpenticular.length();
	double qparallel_l = qparallel.length();
	
//	CylinderCoor3D qs(qptl,0.0,qprl);

	// resolution is in degree, coordinates are in rad
	const int lmax=resolution;
	
	vector<complex<double> > A,B,C,D;
	A.resize(lmax+1,complex<double>(0,0));
	B.resize(lmax+1,complex<double>(0,0));
	C.resize(lmax+1,complex<double>(0,0));
	D.resize(lmax+1,complex<double>(0,0));

	for (Atomselection::iterator asi=as.begin();asi!=as.end();asi++) {
		CartesianCoor3D ct = sample.frames.current().coord3D(*asi);
		
		CartesianCoor3D ctparallel = (o*ct)*o; 
		CartesianCoor3D ctperpenticular = ct - ctparallel; // this contains psi-phi
		double ctperpenticular_l = ctperpenticular.length();
		double ctparallel_l = ctparallel.length();

		double psiphi = 0;

		// if either qper_l or ctper_l is zero , then the bessel_terms vanish and delta_phipsi is irrelevant
		if ((qperpenticular_l!=0) && (ctperpenticular_l!=0)) {
			psiphi = acos( (ctperpenticular * qperpenticular) / (ctperpenticular_l*qperpenticular_l));
			CartesianCoor3D ctq = (ctperpenticular.cross_product(qperpenticular));
			ctq = ctq / ctq.length();				
			if ((ctq + o).length()>1.0) psiphi = 2*M_PI-psiphi; // if o || ctq -> 2; 0 otherwise			
		}

		double esf = atoms[*asi].scatteramp;

		double parallel_sign = 1.0;
		if ((ctparallel_l!=0) && (qparallel_l!=0)) {
			parallel_sign = (ctparallel*qparallel) / (ctparallel_l*qparallel_l);			
		}

//		complex<double> expi = exp(complex<double>(0,c.z*qs.z));		
		complex<double> expi = exp(complex<double>(0,parallel_sign*ctparallel_l*qparallel_l));

//		A[0] += expi * (double)bessj(0,c.r*qs.r) *esf ;
		A[0] += expi * (double)cyl_bessel_j(0,ctperpenticular_l*qperpenticular_l) *esf ;

		for (int l=1;l<=lmax;l++) {
//			complex<double> fac1 = 2.0*powf(-1.0,l)*bessj(2*l,c.r*qs.r);
//			complex<double> fac2 = complex<double>(0,1.0)*double(2.0*powf(-1.0,l-1)*bessj(2*l-1,c.r*qs.r));

			complex<double> fac1 = 2.0*powf(-1.0,l)*cyl_bessel_j(2*l,ctperpenticular_l*qperpenticular_l);
			complex<double> fac2 = complex<double>(0,1.0)*double(2.0*powf(-1.0,l-1)*cyl_bessel_j(2*l-1,ctperpenticular_l*qperpenticular_l));

			A[l] += fac1*expi*cos(2*l*psiphi) *esf;
			B[l] += fac1*expi*sin(2*l*psiphi) *esf;
			C[l] += fac2*expi*cos((2*l-1)*psiphi) *esf;
			D[l] += fac2*expi*sin((2*l-1)*psiphi) *esf;
		}
	}

	// we need to multiply w/ (lmax*4+1) here b/c single ampltiudes are average-summed
	scattering_amplitudes.push_back( sqrt(4*lmax+1)* A[0] ); 		
	double lmax41sqr = sqrt(4*lmax+1)*sqrt(0.5);	
	for (int l=1;l<=lmax;l++) {
		scattering_amplitudes.push_back( lmax41sqr*A[l] ); 		
		scattering_amplitudes.push_back( lmax41sqr*B[l] ); 		
		scattering_amplitudes.push_back( lmax41sqr*C[l] ); 		
		scattering_amplitudes.push_back( lmax41sqr*D[l] ); 		
	}

}

void compute_phase_factors(Sample& sample) {
	// first grid the coordinate system, then calculate the grid-volume for each phase
	// then calculate the unified kappa
	
	
	// for each phase calculate the summed volumes
	vector<ScatteringBackgroundPhaseParameters>& phases = Params::Inst()->scattering.background.phases;
	vector<double> total_volumes(phases.size());

	vector<CoordinateSet> coordinate_sets;
	// first prepare all coordinate set...
	for(size_t i = 0; i < phases.size(); ++i)
	{
		coordinate_sets.push_back( CoordinateSet(sample.frames.current(),sample.atoms.selections[phases[i].selection]) );
	}
	
	for(size_t i = 0; i < phases.size(); ++i)
	{
		Atomselection& phase_atoms = sample.atoms.selections[phases[i].selection];

//		double factor = 1.0;
//		if (phases[i].scaling=="manual") {
//			factor = phases[i].factor;
//		}
		
		double volume = 0.0;
		for(size_t j = 0; j < phase_atoms.size(); ++j)
		{
			volume += sample.atoms[phase_atoms[j]].volume;
		} 
		total_volumes[i] = volume;	

		DensityGrid densgrid(sample.frames.current().unitcell,0.1,sample.frames.current().origin,false);
		densgrid.set(coordinate_sets[i]);
		for(size_t j = 0; j < phases.size(); ++j)
		{
			if (i==j) continue;
			densgrid.unset(coordinate_sets[j],phases[i].nullrange);
		}			
	}
	

	vector<double> grid_volumes(phases.size());
	
}

complex<double> background(Sample& sample,Atomselection& as_system, Atomselection& as_solvent, Atomselection& as_particle, int resolution, double hlayer, CartesianCoor3D q,vector<double>& kappas_s,vector<double>& kappas_p) {
	
//	vector<CartesianCoor3D> uc = sample.frames.current().unitcell;
//	CartesianCoor3D origin = sample.frames.current().origin;
//	
//	vector<int> v; // an empty vector to initialize vGrid3D
//	vGrid3D grid(resolution,uc,origin,v);
//	vGrid3D solgrid(resolution,uc,origin,v);	
//	vGrid3D partgrid(resolution,uc,origin,v);	
//	Grid3D<bool> solventgrid(resolution,uc,origin,false);	
//	Grid3D<bool> particlegrid(resolution,uc,origin,false);		
//	Grid3D<bool> systemgrid(resolution,uc,origin,false);
//
//	for(Atomselection::iterator asi=as_system.begin();asi!=as_system.end();asi++) {
//		CartesianCoor3D c = sample.frames.current().coord3D(*asi);
//		Gridkey3D gk = grid.get_cell(c);
//		grid[gk].push_back(*asi);
//	}	
//		
//	for(Atomselection::iterator asi=as_particle.begin();asi!=as_particle.end();asi++) {
//		CartesianCoor3D c = sample.frames.current().coord3D(*asi);
//		CartesianCoor3D hc(hlayer,hlayer,hlayer);
//		vector<Gridkey3D> gks =particlegrid.get_cells(c,hc);
//		
//		for (vector<Gridkey3D>::iterator gki=gks.begin();gki!=gks.end();gki++) {
//			particlegrid[*gki]=true;			
//		}
//	}	
//	for(Atomselection::iterator asi=as_solvent.begin();asi!=as_solvent.end();asi++) {
//		CartesianCoor3D c = sample.frames.current().coord3D(*asi);
//		Gridkey3D gk = grid.get_cell(c);
//		
//		solventgrid[gk]=true;
//	}
//	
//	// remove valid solvent grid cells...
//	for(Grid3D<bool>::iterator gi=particlegrid.begin();gi!=particlegrid.end();gi++) {
//		if (gi->second) {
//			solventgrid[Gridkey3D(gi->first.ix,gi->first.iy,gi->first.iz)]=false;
//		}
//	}	
//	
//	// build an atomselection...
//	Atomselection as_solvent_trunc(sample.atoms,false);
//	int numsolatoms=0;
//	double solventvolume=0;
//	for(Grid3D<bool>::iterator gi=solventgrid.begin();gi!=solventgrid.end();gi++) {
//		if (gi->second) {
//			vector<int>& v = grid[Gridkey3D(gi->first.ix,gi->first.iy,gi->first.iz)];
//			for (vector<int>::iterator vi=v.begin();vi!=v.end();vi++) {
//				as_solvent_trunc.add(sample.atoms[*vi]);
//				numsolatoms++;
//			}
//			solventvolume+= solventgrid.element_volume();
//		}
//	}	
//	
//	// calculate scattering amplitude...
//	set_scatteramp(sample,as_solvent_trunc,q,false);		
////	complex<double> Ac = scatter(*(as_solvent_trunc.sample),as_solvent_trunc,q);
//
//	complex<double> Ac = scatter_none(sample,as_solvent_trunc,q);	
////	cout << "scatter numatoms=" << numsolatoms << " , " << as_solvent_trunc.size() << endl;
////	cout << "scatter solvent_trunc=" << Ac << endl;
//	
////	cout << "volume:" << totalvolume << endl;
////	cout << "scatter/volume:" << Ac/totalvolume << endl;
//	
//		
//	// build a solvent exclusive grid..
//
//	for(Atomselection::iterator asi=as_solvent.begin();asi!=as_solvent.end();asi++) {
//		CartesianCoor3D c = sample.frames.current().coord3D(*asi);
//		Gridkey3D gk = solgrid.get_cell(c);
//		
//		solgrid[gk].push_back(*asi);
//	}
//	
//
//	// build a particle exclusive grid..
//	for(Atomselection::iterator asi=as_particle.begin();asi!=as_particle.end();asi++) {
//		CartesianCoor3D c = sample.frames.current().coord3D(*asi);
//		Gridkey3D gk = partgrid.get_cell(c);
//		
//		partgrid[gk].push_back(*asi);
//	}
//
//	//scan solventgrid(bool) and grid(atoms) for any outer solvent molecules...
//	double sum_shell_solvent=0;
//	int    n_shell_solvent=0;
//		int num_solcells=0;
//		
//		double dens = 0;
//		
//	for(Grid3D<bool>::iterator gi=solventgrid.begin();gi!=solventgrid.end();gi++) {
//		if (gi->second) {
//			num_solcells++;
//			vector<int>& v = solgrid[Gridkey3D(gi->first.ix,gi->first.iy,gi->first.iz)];
//			int ln =0;
//			for (vector<int>::iterator vi=v.begin();vi!=v.end();vi++) {
//				sum_shell_solvent+=sample.atoms[*vi].volume;
//				n_shell_solvent++;
//				ln++;
//			}
//			dens += 1e27*ln/(solventgrid.element_volume()) / 6.0221415e23 /3;
//		}
//	}	
//
//	// scan particlegrid(bool) and solgrid(atoms) for inner solvent molecules...
//	double hparticlevolume=0;
//	double sum_all_solvent=0;
//	double sum_inner_solvent=0;	
//	for(Grid3D<bool>::iterator gi=particlegrid.begin();gi!=particlegrid.end();gi++) {
//		vector<int>& v = solgrid[Gridkey3D(gi->first.ix,gi->first.iy,gi->first.iz)];		
//
//		if (gi->second) {
//			hparticlevolume+= solgrid.element_volume();
//			for (vector<int>::iterator vi=v.begin();vi!=v.end();vi++) {
//				sum_inner_solvent+=sample.atoms[*vi].volume;
//			}
//		}
//		for (vector<int>::iterator vi=v.begin();vi!=v.end();vi++) {
//			sum_all_solvent+=sample.atoms[*vi].volume;
//		}		
//	}	
//	
//	// scan particlegrid(bool) and partgrid(atoms) for particle molecules...
//	double sum_particle = 0;
//	for(Grid3D<bool>::iterator gi=particlegrid.begin();gi!=particlegrid.end();gi++) {
//		if (gi->second) {
//			vector<int>& v = partgrid[Gridkey3D(gi->first.ix,gi->first.iy,gi->first.iz)];
//			for (vector<int>::iterator vi=v.begin();vi!=v.end();vi++) {
//				sum_particle+=sample.atoms[*vi].volume;
//			}
//		}
//	}	
//	
//	
//	// kappa estimation:
//	double totalvolume = grid.total_volume();
//	double kappa_s = solventvolume/sum_shell_solvent;
//	double particle_volume = totalvolume - sum_all_solvent * kappa_s;
////	double particle_volume = totalvolume - sum_shell_solvent * kappa_s;	
//	double kappa_p = particle_volume / sum_particle;
////	double kappa_p = particle_volume / sum_particle;	
//
////  cout << resolution << "\t" << hlayer << "\t" << dens/num_solcells << endl;
////  
////  cout << "number of solvent cells: " << num_solcells << endl;
////  cout << "n_shell_solvent: " << n_shell_solvent << endl;
////  cout << "numsolatoms: " << numsolatoms << endl;	
////  cout << "solventvolume: " << solventvolume << endl;
////	cout << "cell volume: " << xd*yd*zd << endl;
////  if (numsolatoms!=0)	cout << "solvent particles per cell: " << numsolatoms/num_solcells << endl;	
////  cout << "average density of bulk water (mol./L): " << 1e27*n_shell_solvent/solventvolume / 6.0221415e23 /3  << endl;
////
////  cout << "totalvolume=" << totalvolume << endl;
////  cout << "solventvolume=" << solventvolume << endl;
////  cout << "hparticlevolume=" << hparticlevolume << endl;
////
////  cout << "sum_shell_solvent=" << sum_shell_solvent << endl;
////  cout << "sum_particle=" << sum_particle << endl;
////  cout << "particle_volume=" << particle_volume << endl;
////
////  cout << "sum_inner_solvent=" << sum_inner_solvent << endl;
////  cout << "kappa_s=" << kappa_s << endl;
////  cout << "kappa_p=" << kappa_p << endl;	
////	cout << resolution << "\t" << hlayer << "\t" << kappa_s <<  "\t" << kappa_p << endl;
//
//	// currently we have only a 2 phase system (allowed)
//	// return kappa values via argument
//	kappas_s.push_back(kappa_s);
//	kappas_p.push_back(kappa_p);
//	
////   if (Settings::get("main")["scattering"]["background"].exists("volume")) {
////   	// store kappa values for atoms:
////   	for(Atomselection::iterator asi=as_solvent.begin();asi!=as_solvent.end();asi++) {
////   //		(*asi)->kappa = kappa_s;
////   		sample.atoms[*asi].kappa = 1.0;
////   	}
////   	for(Atomselection::iterator asi=as_particle.begin();asi!=as_particle.end();asi++) {
////   //		(*asi)->kappa = kappa_p;
////   		sample.atoms[*asi].kappa = 1.0;
////   	}		
////   }
//	
	return 0.0;
//	return Ac/solventvolume;
		
}

pair<double,double> background_avg(Sample& sample,int resolution,double hydration, CartesianCoor3D q) {
		
	accumulator_set<double,features<tag::mean,tag::variance> > acc_RE;
	accumulator_set<double,features<tag::mean,tag::variance> > acc_IM;

	vector<double> kappas_s;
	vector<double> kappas_p;
		
	
//	vector<string>& gi = Params::Inst()->scattering.background.include; // background groups
//	vector<string>& ge = Params::Inst()->scattering.background.exclude; // particle groups
////
//	clog << "INFO>> " << "Solvent defined by: ";
//		for(size_t i = 0; i < gi.size(); ++i)
//		{
//			clog << " " << gi[i];
//		}	
//	 clog << endl;
//
//	clog << "INFO>> " << "Solute defined by: ";
//		for(size_t i = 0; i < ge.size(); ++i)
//		{
//			clog << " " << ge[i];
//		}	
//	 clog << endl;
//		
//		
//	Atomselection as_particle(sample.atoms,false);
//	for(size_t i = 0; i < ge.size(); ++i)
//	{
//		sample.atoms.assert_selection(ge[i]);
//		as_particle += sample.atoms.selections[ge[i]];
//	}
//	Atomselection as_solvent(sample.atoms,false);
//	for(size_t i = 0; i < gi.size(); ++i)
//	{
//		sample.atoms.assert_selection(gi[i]);
//		as_solvent += sample.atoms.selections[gi[i]];
//	}
//	sample.atoms.assert_selection("system");
//	Atomselection as_system(sample.atoms,true);
//		
//	int framestride = Params::Inst()->scattering.background.framestride;
//	
//	int totalframes = (sample.frames.size()/framestride) + 1;
//	if ((sample.frames.size() % framestride) ==0) totalframes--;
//	clog << "INFO>> " << "Average background density over " << totalframes << " frames" << endl;
//
//	clog << "INFO>> " << "Progress(%): ";
//	int progress=0;
//	for(size_t i=0;i<sample.frames.size();i+= framestride) {
//		int percent = (i*100/sample.frames.size());
//		if ( percent >= progress  ) {
//			progress = 10*(percent/10) + 10;
//			clog << percent << " ";
//		}
//				
//		sample.frames.load(i,sample.atoms);
//   	 	complex<double> b0v = Analysis::background(sample,as_system,as_solvent,as_particle,resolution,hydration,q,kappas_s,kappas_p);
//		acc_RE(b0v.real()); acc_IM(b0v.imag());
//	}
//	clog << "100" << endl;
//
//	// write kappa values into the atom entries...
//	// create temp. acc
//	// first kappa will contain solvent group, others particles.
//	accumulator_set<double,features<tag::mean,tag::variance> > ks;
//	for (vector<double>::iterator ki=kappas_s.begin();ki!=kappas_s.end();ki++) {
//		ks(*ki);
//	}
//	accumulator_set<double,features<tag::mean,tag::variance> > kp;
//	for (vector<double>::iterator ki=kappas_p.begin();ki!=kappas_p.end();ki++) {
//		kp(*ki);
//	}
//	
//	double m; double sv;
//	m = mean(ks); sv = sqrt(variance(ks));
//	clog << "INFO>> " << "atomic volume scaling factor for solvent groups: " << m << " +- " << sv << endl;
//	for (Atomselection::iterator asi=as_solvent.begin();asi!=as_solvent.end();asi++) {
//		sample.atoms[*asi].kappa = m;
//	}
//
//	m = mean(kp); sv = sqrt(variance(kp));
//	clog << "INFO>> " << "atomic volume scaling factor for solute groups: " << m << " +- " << sv << endl;
//	for (Atomselection::iterator asi=as_particle.begin();asi!=as_particle.end();asi++) {
//		sample.atoms[*asi].kappa = m;
//	}

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
