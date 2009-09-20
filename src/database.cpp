// direct header
#include "database.hpp"

// standard header
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

// special library headers
#include <boost/filesystem.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/regex.hpp>

// other headers
#include "log.hpp"
#include "parameters.hpp"

using namespace std;

void Database::read_conf(std::string filename) {
		
	// configuration file objects are not referenced/destroyed. mind the memory leak if using 'a lot'.
	libconfig::Config* pconf = new libconfig::Config;

	// store the database line by line into an internal buffer, 
	// this is for keeping history
	ifstream dbfile(filename.c_str());
	string line;
	while (getline(dbfile,line)) {
		carboncopy.push_back(line);
	}
	dbfile.close();

	// construct the libconfig interface
	try { pconf->readFile(filename.c_str()); } 
		catch (libconfig::ParseException &e) { 
		Err::Inst()->write(string("Parse exception (line ")+to_s(e.getLine())+string("): ")+string(e.getError()));
		delete pconf; 
		throw;
	}

	libconfig::Setting&	rootsetting = pconf->getRoot();

	// START OF database section //

	if (rootsetting["names"].exists("pdb")) {
		for (int i=0;i<rootsetting["names"]["pdb"].getLength();i++) {	
			string label = static_cast<string>(rootsetting["names"]["pdb"][i].getName());
			string regexp = static_cast<string>((const char*) rootsetting["names"]["pdb"][i]);
			names.pdb.reg(label,regexp);
		}		
	}

	for (int i=0;i<rootsetting["masses"].getLength();i++) {	
		string label = static_cast<string>(rootsetting["masses"][i].getName());
		double mass = static_cast<double>(rootsetting["masses"][i]);
		size_t ID = atomIDs.get(label);
		masses.reg(ID,mass);
	}		

	for (int i=0;i<rootsetting["sizes"].getLength();i++) {	
		
		string label = static_cast<string>(rootsetting["sizes"][i].getName());
		size_t sizes_function_type = static_cast<long>(rootsetting["sizes"][i][0]);
		
		vector<double> values;
		for(size_t k = 1; k < rootsetting["sizes"][i].getLength(); ++k)
		{
			values.push_back(static_cast<double>(rootsetting["sizes"][i][k]));
		}					
		
		size_t ID = atomIDs.get(label);
		volumes.reg(ID,values,sizes_function_type);
	}	
	
	// determine which probe-type to load: 
	bool probe_found=false;
	for (int i=0;i<rootsetting["sfactors"].getLength();i++) {	
		if (static_cast<string>(rootsetting["sfactors"][i].getName())== Params::Inst()->scattering.scatterfactors) {
		
			for (int j=0;j<rootsetting["sfactors"][i].getLength();j++) {	
				
				string label = static_cast<string>(rootsetting["sfactors"][i][j].getName());
				size_t sfactors_function_type = static_cast<long>(rootsetting["sfactors"][i][j][0]);
				size_t ID = atomIDs.get(label);
				
				vector<double> values;
				for(size_t k = 1; k < rootsetting["sfactors"][i][j].getLength(); ++k)
				{
					values.push_back(static_cast<double>(rootsetting["sfactors"][i][j][k]));
				}					

				sfactors.reg(ID,values,sfactors_function_type);
			}			
			
			probe_found = true;
			break;
		}
	}
	if (!probe_found) throw;
	
	
	// END OF database section //

	delete pconf;
};

string Database::guessformat(string filename) {
	// do the best you can to guess the format
	return "conf";
}

bool Database::check() {
	// implement a sanity check for the database
}

void Database::init(std::string filename,std::string format) {
	if (format=="") {
		Info::Inst()->write(string("Detecting format of database file: ") + filename);
		format = guessformat(filename);
		Info::Inst()->write(string("Detected format: ") + format);
	}
	
	if (format=="conf") {
		read_conf(filename);
	}
	if (format=="xml") {
		Err::Inst()->write("XML Database not yet implemented");
		throw;
	}

	Info::Inst()->write("Checking parameters...");	
	
	if (check()) {
		Info::Inst()->write("Check succeeded. Parameters seem OK.");			
	}
	else {
		Err::Inst()->write("Check failed. Please check your input file.");
		throw;		
	}
	
}

void Database::write(std::string filename, std::string format) {
	if (format=="") {
		Info::Inst()->write("No database format specified. Will use conf.");
		format = "conf";
	}
	
	if (format=="conf") {
		Err::Inst()->write("conf database format (write) not yet implemented.");
		throw;
	}
	if (format=="xml") {
		Err::Inst()->write("XML database format (write) not yet implemented.");
		throw;
	}
}


void DatabaseNamesPDBParameters::reg(string label, string regexp) {
	m_label2regexp[label]=regexp;
}

void DatabaseNamesPDBParameters::add_quicklookup(string testlabel,string label) {
	if (quicklookup_counter>1000000) clear_quicklookup();
	m_quicklookup[testlabel]=label;
	quicklookup_counter++;	
}

void DatabaseNamesPDBParameters::clear_quicklookup() {
	m_quicklookup.clear();
	quicklookup_counter=0;
}

string DatabaseNamesPDBParameters::get(string testlabel) {
		
   int matches = 0;
   string label = "<Unknown>";

	if (m_quicklookup.find(testlabel)!=m_quicklookup.end()) {
		return m_quicklookup[testlabel];
	}
	
	for(map<string,string>::iterator lri = m_label2regexp.begin();lri!=m_label2regexp.end();lri++ ) {
	   	boost::regex e(lri->second);
	
	   	if (regex_match(testlabel,e)) {
	   		if (matches>0) {
				Err::Inst()->write(string("PDB atom name matches")+string(label)+string(" and ")+string(lri->first));
				Err::Inst()->write(string("expression to match: ")+string(testlabel));
	   			throw;
	   		} else 
			{ 
				label = lri->first; 
				matches++; 
			}
	   	}
	}

   if (matches==0) {
	Err::Inst()->write(string("PDB atom name not recognized: ")+string(testlabel));
   	throw;
   }
   
	m_quicklookup[testlabel]=label;
	return label;	
}



void DatabaseAtomIDsParameters::reg(string label) {
	if (m_IDs.find(label)==m_IDs.end()) {
		m_IDs[label]=nextID;
		m_rIDs[nextID]=label;
		nextID++;
	}
}

size_t DatabaseAtomIDsParameters::get(string label) {
	if (m_IDs.find(label)==m_IDs.end()) {
		reg(label);
	}
	return m_IDs[label];
}

string DatabaseAtomIDsParameters::rget(size_t ID) {
	return m_rIDs[ID];
}


void DatabaseMassesParameters::reg(size_t ID, double value) {
	m_masses[ID]=value;
}

double DatabaseMassesParameters::get(size_t ID) {
	return m_masses[ID];
}




void DatabaseVolumesParameters::reg(size_t ID, vector<double> constants, size_t function_type) {
	m_constants[ID]=constants;
	m_functiontypes[ID] = function_type;
}

void DatabaseVolumesParameters::add_quicklookup(size_t ID,double value) {
	if (quicklookup_counter>1000000) clear_quicklookup();
	m_quicklookup[ID]=value;
	quicklookup_counter++;	
}

void DatabaseVolumesParameters::clear_quicklookup() {
	m_quicklookup.clear();
	quicklookup_counter=0;
}

double DatabaseVolumesParameters::get(size_t ID) {
		
	size_t ft = m_functiontypes[ID];
	vector<double> v = m_constants[ID];
	if (ft==0) { // 0 = constants // volume direkt
		return v[0];
	}
	
	if (ft==1) { // hardsphere radius
		return (4.0/3.0)*M_PI*powf(v[0],3);
	}
	
	if (ft==2) { // gausssphere radius
		return sqrt(powf(M_PI,3))*powf(v[0],3);
	}

	// don't need a quicklookup for volume calculations, b/c pretty fast
}




void DatabaseSFactorsParameters::reg(size_t ID, vector<double> constants, size_t function_type) {
	m_constants[ID]=constants;
	m_functiontypes[ID] = function_type;
}

void DatabaseSFactorsParameters::add_quicklookup(size_t ID, double q,double value) {
	if (quicklookup_counter>1000000) clear_quicklookup();
	m_quicklookup[q][ID]=value;
	quicklookup_counter++;	
}

void DatabaseSFactorsParameters::clear_quicklookup() {
	m_quicklookup.clear();
	quicklookup_counter=0;
}

double DatabaseSFactorsParameters::get(size_t ID,double q) {
		
	size_t& ft = m_functiontypes[ID];
	vector<double>& v = m_constants[ID];
	if (ft==0) { // 0 = constants
		return v[0];
	}
	
	// a shortcut for expensive calculations, this thing grows. 
	if ((ft==1) || (ft==2)) {
		if ((m_quicklookup.find(q)!=m_quicklookup.end()) && (m_quicklookup[q].find(ID)!=m_quicklookup[q].end())) {
			return m_quicklookup[q][ID];
		}
	}
	
	if (ft==1) {
		// slater:
		//     den=0.0
		//      do j=1,5
		//        if (coef(z,1,j).ne.0.0) then
		//          c = ((2.*coef(z,1,j)/coef(z,3,j))**(coef(z,3,j)+0.5))/
		//     1        sqrt(fak(int(2.*coef(z,3,j))))
		//          den = den+coef(z,2,j)*(      c  *  r**(coef(z,3,j)-1.0)  *  exp(-coef(z,1,j)  *  r/coef(z,3,j)/a)       )**2
		//     1          exp(-coef(z,1,j)*r/coef(z,3,j)/a))**2
		//        end if
		//      end do
		//      den = den/4./pi
		double den = 0.0;
		for(size_t j = 0; j < 5; ++j)
		{
			size_t j2 = 5+j;
			size_t j3 = 10+j;
			if (v[j]!=0.0) {
				double c = 0.0;
				c += powf( 2.0*v[j]/v[j3], v[j3]+0.5 );
				c /= sqrt( boost::math::factorial<double>(2.0*v[j3]) );
				den += v[j2]* powf(  c*powf(q, v[j3]-1.0)*exp(-1.0*v[j]*q/v[j3]/1.0) , 2);
			}
		}
		
		add_quicklookup(ID,q,den / (4.0*M_PI));
		return (den / (4.0*M_PI));
	}
	
	if (ft==2) {
		double den = 0.0;
		den += v[0]*exp(-v[1]*q*q);
		den += v[2]*exp(-v[3]*q*q);
		den += v[4]*exp(-v[5]*q*q);
		den += v[6]*exp(-v[7]*q*q);
		den += v[8];

		add_quicklookup(ID,q,den);
		return den;
	}
}