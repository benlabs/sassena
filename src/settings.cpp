/*
 *  settings.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

// direct header
#include "settings.hpp"

// standard header
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

// special library headers
#include <libconfig.h++>
#include <boost/filesystem.hpp>

// other headers

using namespace std;
using namespace boost::filesystem;

std::vector<std::string> Settings::conf_files;
std::map<std::string,libconfig::Setting*> Settings::settings;

bool Settings::conf_files_set = false;
string Settings::config_rootpath = "";

inline double gaussvolume(double r) { return sqrt(powf(M_PI,3))*powf(r,3); }
inline double hardspherevolume(double r) { return (4.0/3.0)*M_PI*powf(r,3); }

libconfig::Setting* Settings::read_from_configuration(string filename,bool expandfilename = true) {
		// configuration file objects are not referenced/destroyed. mind the memory leak if using 'a lot'.
		libconfig::Config* pconf = new libconfig::Config;
		
		string truefilename = expandfilename ? get_filepath(filename) : filename;
		
		try { pconf->readFile(truefilename.c_str()); } 
		catch (libconfig::ParseException &e) { 
			cerr << "Parse exception (line " << e.getLine() << "): " << e.getError() << endl;
			delete pconf; 
			throw;
		}
//		catch (...) {
//			delete pconf; 
//			throw;
//		}

		return &(pconf->getRoot());
}


void Settings::unfold_settings(libconfig::Setting* setting) {
	
	if (!conf_files_set) {
		conf_files_set = true;
		conf_files.push_back("scattering_factors");
		conf_files.push_back("atom_names");
		conf_files.push_back("atom_masses");
		conf_files.push_back("excluded_volumes");	
	}
	
	for (int i=0;i<setting->getLength();i++) {
		string settings_name = (*setting)[i].getName();
		if ((*setting)[i].getType()==libconfig::Setting::TypeGroup) {
			unfold_settings(&((*setting)[i]));
		} else {
			for (vector<string>::iterator cfi=conf_files.begin();cfi!=conf_files.end();cfi++) {
				if (settings_name == *cfi) {
					string fn = (*setting)[i];
					settings[*cfi] = read_from_configuration(fn.c_str()); 
				}
			}
		}
	}
}

bool Settings::read(int& argc,char** argv) {
	bool status = true;
	
	string maindir;
	
	// single configuration file mode
	if (argc==1) {
		status = false;
	}
	else if (argc==2) {
		// make the directory of the main configuration file the root for all others
		if (path(argv[1]).is_complete()) 
			config_rootpath = path(argv[1]).parent_path().string();
		else 
			config_rootpath = ( initial_path() / path(argv[1]).parent_path() ).string();
		settings["main"] = read_from_configuration(argv[1],false);
	} else {
		cerr << "Command line mode currently not supported. Use a configuration file." << endl;
		status = false;
	}

	if (status) {
		// some entries are paths for yet another configuration file (powered by keyword list)
		// these entries have to be detected and loaded into the main map of settings
		unfold_settings(settings["main"]);	
	} else {
		cerr << "SASSIM usage:                                    " << endl;
		cerr << "configuration file mode: sassim sassim.conf      " << endl;
		cerr << "                                                 " << endl;
	}

	return status;
}

std::string Settings::get_filepath(string fname) {
	path fpath(fname);
	string fdir;
	if (fpath.parent_path().is_complete()) 
		fdir = fpath.parent_path().string(); 
	else if (!config_rootpath.empty())
		fdir = (path(config_rootpath) / fpath.parent_path()).string();
	else 
		fdir = (initial_path() / fpath.parent_path()).string();
	return (path(fdir) / fpath.filename()).string();
}

double Settings::get_volume(string atomname) {
	
	string volident = (const char *) Settings::get("main")["scattering"]["background"]["volumes"];			
	string volmethod = 	Settings::get("excluded_volumes")[volident]["method"];
	
	if (volmethod == "gaussradius") {
	 	return gaussvolume(Settings::get("excluded_volumes")[volident][atomname]);	
	}
	else if (volmethod == "hardradius") {
	 	return hardspherevolume(Settings::get("excluded_volumes")[volident][atomname]);	
	}
	else if (volmethod == "volume") {
	 	return Settings::get("excluded_volumes")[volident][atomname];	
	}
	else {
		cerr << "ERROR>> " << " Volume method not understood." << endl;
		throw;
	}
}

void Settings::get_qqqvectors(std::vector<CartesianCoor3D>& qqqvectors) {
	
	// lookup section to use:
	string qqqvectors_method = Settings::get("main")["scattering"]["vectors"]["method"];
	
	
	if (qqqvectors_method=="single") {	
		libconfig::Setting& s = Settings::get("main")["scattering"]["vectors"]["single"];
		double x = s["direction"][0];
		double y = s["direction"][1];
		double z = s["direction"][2];		
		CartesianCoor3D direction(x,y,z);
		
		double scal = s["scale"];				
		qqqvectors.push_back(scal*direction);
	}
	else if (qqqvectors_method=="linear") {	
		libconfig::Setting& s = Settings::get("main")["scattering"]["vectors"]["linear"];
		double x = s["direction"][0];
		double y = s["direction"][1];
		double z = s["direction"][2];		
		CartesianCoor3D direction(x,y,z);
				
		double from = s["from"];		
		double to = s["to"];				
		int points = s["points"];
		double exponent = 1.0;
		if (s.exists("exponent")) exponent = s["exponent"];
		
		for (int i=0;i<points;i++) {		
			double scal = from + powf((i*1.0/(points-1)),exponent)*(to-from);
			qqqvectors.push_back(scal*direction);
		}
	}
	else if (qqqvectors_method=="map") {
		libconfig::Setting& setting = Settings::get("main")["scattering"]["vectors"]["map"];

		libconfig::Setting& s1 = setting[0];
		libconfig::Setting& s2 = setting[1];
		
		int p1 = s1["points"];
		int p2 = s2["points"];
		double from1 = s1["from"];
		double from2 = s2["from"];
		double to1 = s1["to"];
		double to2 = s2["to"];
		double exponent1 = 1.0;
		double exponent2 = 1.0;		
		if (s1.exists("exponent")) exponent1 = s1["exponent"];
		if (s2.exists("exponent")) exponent2 = s2["exponent"];
		
		double x1 = s1["direction"][0];
		double y1 = s1["direction"][1];
		double z1 = s1["direction"][2];		
		CartesianCoor3D d1(x1,y1,z1);
		double x2 = s2["direction"][0];
		double y2 = s2["direction"][1];
		double z2 = s2["direction"][2];		
		CartesianCoor3D d2(x2,y2,z2);
		
		for (int i=0;i<p1;i++) {		
			double scal1 = from1 + powf((i*1.0/(p1-1)),exponent1)*(to1-from1);
			for (int j=0;j<p2;j++) {			
				double scal2 = from2 + powf((j*1.0/(p2-1)),exponent2)*(to2-from2);
				qqqvectors.push_back(scal1*d1+scal2*d2);
			}
		}
	}
	else if (qqqvectors_method=="cube") {
		libconfig::Setting& setting = Settings::get("main")["scattering"]["vectors"]["cube"];

		libconfig::Setting& s1 = setting[0];
		libconfig::Setting& s2 = setting[1];
		libconfig::Setting& s3 = setting[2];
		
		int p1 = s1["points"];
		int p2 = s2["points"];
		int p3 = s3["points"];		
		double from1 = s1["from"];
		double from2 = s2["from"];
		double from3 = s3["from"];		
		double to1 = s1["to"];
		double to2 = s2["to"];
		double to3 = s3["to"];
		double exponent1 = 1.0;
		double exponent2 = 1.0;		
		double exponent3 = 1.0;				
		if (s1.exists("exponent")) exponent1 = s1["exponent"];
		if (s2.exists("exponent")) exponent2 = s2["exponent"];
		if (s3.exists("exponent")) exponent3 = s3["exponent"];
				
		double x1 = s1["direction"][0];
		double y1 = s1["direction"][1];
		double z1 = s1["direction"][2];		
		CartesianCoor3D d1(x1,y1,z1);
		double x2 = s2["direction"][0];
		double y2 = s2["direction"][1];
		double z2 = s2["direction"][2];		
		CartesianCoor3D d2(x2,y2,z2);
		double x3 = s3["direction"][0];
		double y3 = s3["direction"][1];
		double z3 = s3["direction"][2];		
		CartesianCoor3D d3(x3,y3,z3);
				
		for (int i=0;i<p1;i++) {		
			double scal1 = from1 + powf((i*1.0/(p1-1)),exponent1)*(to1-from1);
			for (int j=0;j<p2;j++) {			
				double scal2 = from2 + powf((j*1.0/(p2-1)),exponent2)*(to2-from2);				
				for (int k=0;k<p3;k++) {			
					double scal3 = from3 + powf((k*1.0/(p3-1)),exponent3)*(to3-from3);									
					qqqvectors.push_back(scal1*d1+scal2*d2+scal3*d3);
				}
			}
		}
	}	
	else if (qqqvectors_method=="file") {
		string qqqfilename = Settings::get("main")["scattering"]["vectors"]["file"];
		ifstream qqqfile(Settings::get_filepath(qqqfilename).c_str());
		
		double x,y,z; 
		while (qqqfile >> x >> y >> z) {
			qqqvectors.push_back(CartesianCoor3D(x,y,z));
		}
	}
}

// end of file
