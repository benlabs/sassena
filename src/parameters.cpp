// direct header
#include "parameters.hpp"

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

using namespace std;



// provide a set of overloaded function to wrap the type-dependency of the configuration file params
std::string getstring(const char * value) { stringstream ss; ss << value; return ss.str(); }
std::string getstring(std::string value)  { stringstream ss; ss << value; return ss.str(); }
std::string getstring(double value)       { stringstream ss; ss << value; return ss.str(); }
std::string getstring(long value)         { stringstream ss; ss << value; return ss.str(); }
std::string getstring(bool value)         { stringstream ss; ss << value; return ss.str(); }
std::string getstring(libconfig::Setting& s) { 
	if (s.getType()==libconfig::Setting::TypeBoolean) {
		stringstream ss; ss << bool(s); return ss.str();		
	}
	if (s.getType()==libconfig::Setting::TypeFloat) {
		stringstream ss; ss << double(s); return ss.str();		
	}
	if (s.getType()==libconfig::Setting::TypeInt) {
		stringstream ss; ss << long(s); return ss.str();		
	}	
	if (s.getType()==libconfig::Setting::TypeInt64) {
		stringstream ss; ss << long(s); return ss.str();		
	}
	if (s.getType()==libconfig::Setting::TypeString) {
		string str = s;
		stringstream ss; ss << str; return ss.str();		
	}
	Err::Inst()->write("error in getstring(Setting&)");
	throw;
}

double getdouble(const char * value) { return atof(value); }
double getdouble(std::string value)  { return atof(value.c_str()); }
double getdouble(double value)       { return value; }
double getdouble(long value)         { return double(value); }
double getdouble(bool value)         { if (value) return 1.0; else return 0.0; }
double getdouble(libconfig::Setting& s) { 
	if (s.getType()==libconfig::Setting::TypeBoolean) {
		if (s) return 1.0; else return 0.0; 		
	}
	if (s.getType()==libconfig::Setting::TypeFloat) {
		return s;
	}
	if (s.getType()==libconfig::Setting::TypeInt) {
		return double(s);
	}
	if (s.getType()==libconfig::Setting::TypeInt64) {
		return double(s);
	}
	if (s.getType()==libconfig::Setting::TypeString) {
		return atof(s);
	}
	Err::Inst()->write("error in getdouble(Setting&)");
	throw;
}

long getlong(const char * value) { return atol(value); }
long getlong(std::string value)  { return atol(value.c_str()); }
long getlong(double value)       { return long(value); }
long getlong(long value)         { return value; }
long getlong(bool value)         { if (value) return 1; else return 0; }
long getlong(libconfig::Setting& s) { 
	if (s.getType()==libconfig::Setting::TypeBoolean) {
		if (s) return 1; else return 0; 		
	}
	if (s.getType()==libconfig::Setting::TypeFloat) {
		return long(s);
	}
	if (s.getType()==libconfig::Setting::TypeInt) {
		return s;
	}	
	if (s.getType()==libconfig::Setting::TypeInt64) {
		return s;
	}
	if (s.getType()==libconfig::Setting::TypeString) {
		return atol(s);
	}
	if ( (s.getType()==libconfig::Setting::TypeGroup) ||
		 (s.getType()==libconfig::Setting::TypeArray) ||
 		 (s.getType()==libconfig::Setting::TypeList) ) {
			Err::Inst()->write("error in getlong(Setting&): section to basic");
	}
	
	Err::Inst()->write("error in getlong(Setting&)");
	throw;
}

bool getbool(const char * value) { stringstream ss; ss << value; if (ss.str()=="true" || ss.str()=="yes" || ss.str()=="on" || ss.str()=="1") return true; else return false; }
bool getbool(std::string value)  { stringstream ss; ss << value; if (ss.str()=="true" || ss.str()=="yes" || ss.str()=="on" || ss.str()=="1") return true; else return false; }
bool getbool(double value)       { if (value==1.0) return true; else return false; }
bool getbool(long value)         { if (value==1) return true; else return false; }
bool getbool(bool value)         { if (value) return true; else return false; }
bool getbool(libconfig::Setting& s) { 
	if (s.getType()==libconfig::Setting::TypeBoolean) {
		if (s) return true; else return false; 		
	}
	if (s.getType()==libconfig::Setting::TypeFloat) {
		if (double(s)==1.0) return true; else return false;
	}
	if (s.getType()==libconfig::Setting::TypeInt) {
		if (long(s)==1) return true; else return false;
	}
	if (s.getType()==libconfig::Setting::TypeInt64) {
		if (long(s)==1) return true; else return false;
	}
	if (s.getType()==libconfig::Setting::TypeString) {
		string str = s;
		if ((str=="true") || (str=="yes") || (str=="on") || (str=="1")) return true; else return false; 
	}
	Err::Inst()->write("error in getbool(Setting&)");
	throw;
}


//////////////////////////////////////////////////////////////////////////////////////////////////
/////// Class Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////


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
		
	size_t ft = m_functiontypes[ID];
	vector<double> v = m_constants[ID];
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
				c /= sqrt( boost::math::factorial<int>(int(2.0*v[j3])) );
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
		den += v[5]*exp(-v[6]*q*q);
		den += v[7];

		add_quicklookup(ID,q,den);
		return den;
	}
}



string Params::get_filepath(string fname) {
	using namespace boost::filesystem;

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


void Params::read_conf(std::string filename) {
		
	// configuration file objects are not referenced/destroyed. mind the memory leak if using 'a lot'.
	libconfig::Config* pconf = new libconfig::Config;

	string truefilename = get_filepath(filename);

	// store the configuration line by line into an internal buffer, 
	// this is for keeping history
	ifstream conffile(truefilename.c_str());
	string line;
	while (getline(conffile,line)) {
		configuration.push_back(line);
	}
	conffile.close();

	// construct the libconfig interface
	try { pconf->readFile(truefilename.c_str()); } 
		catch (libconfig::ParseException &e) { 
		Err::Inst()->write(string("Parse exception (line ")+to_s(e.getLine())+string("): ")+string(e.getError()));
		delete pconf; 
		throw;
	}

	libconfig::Setting&	rootsetting = pconf->getRoot();
	// some flags , which we will need
	bool flag_nogroups = false;
	
	// START OF sample section //	
	
	// now read the parameters
	sample.structure.file   = get_filepath(getstring(rootsetting["sample"]["structure"]["file"]));
	sample.structure.format = getstring(rootsetting["sample"]["structure"]["format"]);

	// read in group definitions
	if (rootsetting["sample"].exists("groups") && (rootsetting["sample"]["groups"].getLength()>0)) {
		for(size_t i = 0; i < rootsetting["sample"]["groups"].getLength(); ++i)
		{
			libconfig::Setting& sgroup = rootsetting["sample"]["groups"][i];					
			string gn = "";
			if (sgroup.exists("name")) gn = getstring(sgroup["name"]); // groupname
			string fn = get_filepath(getstring(sgroup["file"])); // filename
			string ft = getstring(sgroup["type"]); // filetype
			string sn = "";
			if (sgroup.exists("select")) sn = getstring(sgroup["select"]); // selection field name
			double sv = 0.0;
			if (sgroup.exists("select_value")) sv = getdouble(sgroup["select_value"]); // selection field value, positive match								
			sample.groups[gn] = SampleGroupParameters(gn,fn,ft,sn,sv);
		}
	} else {
		// no group definition
		Warn::Inst()->write("No groups defined. Switching to full system scattering");
		flag_nogroups = true;
	}	

	// apply deuteration
	if (rootsetting["sample"].exists("deuter")) {
		if (rootsetting["sample"]["deuter"].getType()==libconfig::Setting::TypeString) {
			Info::Inst()->write(string("Deuteration of group ") + getstring(rootsetting["sample"]["deuter"]));		
			sample.deuter.push_back(rootsetting["sample"]["deuter"]);
		}
		else if (rootsetting["sample"]["deuter"].getType()==libconfig::Setting::TypeList) {
	    	for (int i=0;i<rootsetting["sample"]["deuter"].getLength();i++) {
				Info::Inst()->write(string("Deuteration of group ") + getstring(rootsetting["sample"]["deuter"][i]));		
				sample.deuter.push_back(rootsetting["sample"]["deuter"][i]);				
			}
		}
	}
	else {
		Info::Inst()->write("No explicit deuteration");		
	}
	
	
	// read in frame information
	for (int i=0;i<rootsetting["sample"]["frames"].getLength();i++) {
		string fn = get_filepath(getstring(rootsetting["sample"]["frames"][i]["file"]));
		string ft = getstring(rootsetting["sample"]["frames"][i]["type"]);
		sample.frames.push_back(make_pair(fn,ft));
		Info::Inst()->write(string("Added frames from ")+fn+string(" using format: ")+ft);				
	}
	
	// periodic boundary behavior and/or postprocessing
	
	if (rootsetting["sample"]["pbc"]["wrapping"]) {
		sample.pbc.wrapping = true;
		sample.pbc.center = "";
		if (rootsetting["sample"]["pbc"].exists("center")) sample.pbc.center = getstring(rootsetting["sample"]["pbc"]["center"]);
		
		Info::Inst()->write(string("Turned wrapping ON with center group ")+sample.pbc.center);
	}
	else {
		sample.pbc.wrapping = false;		
		Info::Inst()->write("Turned wrapping OFF");
	}	
	
	// END OF sample section //
	// START OF scattering section //

	if (rootsetting["scattering"].exists("background")) {
		string m = getstring(rootsetting["scattering"]["background"]["method"]);
		bool rescale = getbool(rootsetting["scattering"]["background"]["rescale"]);
		// rescale automatically triggers background analysis routine. Maybe independent in the future...
		scattering.background.method = m;

		if (m=="auto" || ((m == "fixed") && rescale)) {

			if (m=="fixed") Info::Inst()->write(" rescale=true triggered analysis of the background scattering factor!");

			if (rootsetting["scattering"]["background"].exists("include")) {
				if ( (rootsetting["scattering"]["background"]["include"].getType()==libconfig::Setting::TypeList) ||
					 (rootsetting["scattering"]["background"]["include"].getType()==libconfig::Setting::TypeArray) ) {
						for(size_t i = 0; i < rootsetting["scattering"]["background"]["include"].getLength(); ++i)
						{
							scattering.background.include.push_back(rootsetting["scattering"]["background"]["include"][i]);
						}
					}
			}
			if (rootsetting["scattering"]["background"].exists("exclude")) {
				if ( (rootsetting["scattering"]["background"]["exclude"].getType()==libconfig::Setting::TypeList) ||
					 (rootsetting["scattering"]["background"]["exclude"].getType()==libconfig::Setting::TypeArray) ) {
						for(size_t i = 0; i < rootsetting["scattering"]["background"]["exclude"].getLength(); ++i)
						{
							scattering.background.exclude.push_back(rootsetting["scattering"]["background"]["exclude"][i]);
						}
					}
			}
			// determine background scattering density from grid based analysis
			// necessary:
			// group = name of the group which accounts for the background
			// resolution = grid resolution
			// hydration = hydration layer around each particle, 
			//    ie. solvent molecules within this grid distance are not counted towards bulk
			scattering.background.resolution = getlong(rootsetting["scattering"]["background"]["resolution"]); // resolution of grid
			
			scattering.background.hydration  = getdouble(rootsetting["scattering"]["background"]["hydration"]); // hydration layer of each particle (not background)

			Info::Inst()->write("Background: <auto>");			

			scattering.background.framestride = 1;

			if (rootsetting["scattering"]["background"].exists("framestride")) scattering.background.framestride = getlong(rootsetting["scattering"]["background"]["framestride"]);
		}
		else if (m=="fixed") {
			Info::Inst()->write("Background: <fixed>");
			scattering.background.value = getdouble(rootsetting["scattering"]["background"]["value"]);
		}
		else if (m=="none" || m=="vaccuum") {
			Info::Inst()->write("Background: <vaccuum>");
			scattering.background.method = "fixed";
			scattering.background.value = 0.0;
		}
	}
	else {
		Info::Inst()->write("Background: <vaccuum>");
		scattering.background.method = "fixed";
		scattering.background.value = 0.0;
	}


	// generating qqqvectors, i.e. the spectrum
	string vm = getstring(rootsetting["scattering"]["vectors"]["method"]);
	if (vm=="single") {	
		libconfig::Setting& s = rootsetting["scattering"]["vectors"]["single"];
		double x = getdouble(s["direction"][0]);
		double y = getdouble(s["direction"][1]);
		double z = getdouble(s["direction"][2]);		
		CartesianCoor3D direction(x,y,z);
		
		double scal = getdouble(s["scale"]);				
		scattering.qqqvectors.push_back(scal*direction);
	}
	else if (vm=="linear") {	
		libconfig::Setting& s = rootsetting["scattering"]["vectors"]["linear"];
		double x = getdouble(s["direction"][0]);
		double y = getdouble(s["direction"][1]);
		double z = getdouble(s["direction"][2]);		
		CartesianCoor3D direction(x,y,z);
				
		double from = getdouble(s["from"]  );		
		double to =   getdouble(s["to"]    );				
		int points =  getlong(s["points"]);
		double exponent = 1.0;
		if (s.exists("exponent")) exponent = getdouble(s["exponent"]);
		
		for (int i=0;i<points;i++) {		
			double scal = from + powf((i*1.0/(points-1)),exponent)*(to-from);
			scattering.qqqvectors.push_back(scal*direction);
		}
	}
	else if (vm=="map") {
		libconfig::Setting& setting = rootsetting["scattering"]["vectors"]["map"];
		libconfig::Setting& s1 = setting[0];
		libconfig::Setting& s2 = setting[1];
		
		int p1 = getlong(s1["points"]);
		int p2 = getlong(s2["points"]);
		double from1 = getdouble(s1["from"]);
		double from2 = getdouble(s2["from"]);
		double to1 = getdouble(s1["to"]);
		double to2 = getdouble(s2["to"]);
		double exponent1 = 1.0;
		double exponent2 = 1.0;		
		if (s1.exists("exponent")) exponent1 = getdouble(s1["exponent"]);
		if (s2.exists("exponent")) exponent2 = getdouble(s2["exponent"]);
		
		double x1 = getdouble(s1["direction"][0]);
		double y1 = getdouble(s1["direction"][1]);
		double z1 = getdouble(s1["direction"][2]);		
		CartesianCoor3D d1(x1,y1,z1);
		double x2 = getdouble(s2["direction"][0]);
		double y2 = getdouble(s2["direction"][1]);
		double z2 = getdouble(s2["direction"][2]);		
		CartesianCoor3D d2(x2,y2,z2);
		
		for (int i=0;i<p1;i++) {		
			double scal1 = from1 + powf((i*1.0/(p1-1)),exponent1)*(to1-from1);
			for (int j=0;j<p2;j++) {			
				double scal2 = from2 + powf((j*1.0/(p2-1)),exponent2)*(to2-from2);
				scattering.qqqvectors.push_back(scal1*d1+scal2*d2);
			}
		}
	}
	else if (vm=="cube") {
		libconfig::Setting& setting = rootsetting["scattering"]["vectors"]["cube"];
		libconfig::Setting& s1 = setting[0];
		libconfig::Setting& s2 = setting[1];
		libconfig::Setting& s3 = setting[2];
		
		int p1 = getlong(s1["points"]);
		int p2 = getlong(s2["points"]);
		int p3 = getlong(s3["points"]);		
		double from1 = getdouble(s1["from"]);
		double from2 = getdouble(s2["from"]);
		double from3 = getdouble(s3["from"]);		
		double to1 = getdouble(s1["to"]);
		double to2 = getdouble(s2["to"]);
		double to3 = getdouble(s3["to"]);
		double exponent1 = 1.0;
		double exponent2 = 1.0;		
		double exponent3 = 1.0;				
		if (s1.exists("exponent")) exponent1 = getdouble(s1["exponent"]);
		if (s2.exists("exponent")) exponent2 = getdouble(s2["exponent"]);
		if (s3.exists("exponent")) exponent3 = getdouble(s3["exponent"]);
				
		double x1 = getdouble(s1["direction"][0]);
		double y1 = getdouble(s1["direction"][1]);
		double z1 = getdouble(s1["direction"][2]);		
		CartesianCoor3D d1(x1,y1,z1);
		double x2 = getdouble(s2["direction"][0]);
		double y2 = getdouble(s2["direction"][1]);
		double z2 = getdouble(s2["direction"][2]);		
		CartesianCoor3D d2(x2,y2,z2);
		double x3 = getdouble(s3["direction"][0]);
		double y3 = getdouble(s3["direction"][1]);
		double z3 = getdouble(s3["direction"][2]);		
		CartesianCoor3D d3(x3,y3,z3);
				
		for (int i=0;i<p1;i++) {		
			double scal1 = from1 + powf((i*1.0/(p1-1)),exponent1)*(to1-from1);
			for (int j=0;j<p2;j++) {			
				double scal2 = from2 + powf((j*1.0/(p2-1)),exponent2)*(to2-from2);				
				for (int k=0;k<p3;k++) {			
					double scal3 = from3 + powf((k*1.0/(p3-1)),exponent3)*(to3-from3);									
					scattering.qqqvectors.push_back(scal1*d1+scal2*d2+scal3*d3);
				}
			}
		}
	}	
	else if (vm=="file") {
		string qqqfilename = get_filepath(getstring(rootsetting["scattering"]["vectors"]["file"]));
		ifstream qqqfile(qqqfilename.c_str());
		
		double x,y,z; 
		while (qqqfile >> x >> y >> z) {
			scattering.qqqvectors.push_back(CartesianCoor3D(x,y,z));
		}
	}

	scattering.framestride = 1;
	if (rootsetting["scattering"].exists("framestride")) scattering.framestride = getlong(rootsetting["scattering"]["framestride"]);

	if (rootsetting["scattering"].exists("correlation")) {
		scattering.correlation.type = getstring(rootsetting["scattering"]["correlation"]["type"]);
	}
	else {
		scattering.correlation.type = "none";
	}

	scattering.target = getstring(rootsetting["scattering"]["target"]);
	
	if (rootsetting["scattering"].exists("average")) {
		scattering.average.type   = getstring(rootsetting["scattering"]["average"]["type"]);
		if (scattering.average.type!="none") {
			scattering.average.method = getstring(rootsetting["scattering"]["average"]["method"]);
			if (scattering.average.method=="bruteforce") {		
				scattering.average.vectors = getstring(rootsetting["scattering"]["average"]["vectors"]);
			}
			scattering.average.resolution =  getdouble(rootsetting["scattering"]["average"]["resolution"]);
		}
	}
	else {
		scattering.average.type = "none";
	}
	
	if (rootsetting["scattering"].exists("interference")) {
		scattering.interference.type = getstring(rootsetting["scattering"]["interference"]["type"]);
	}
	else {
		scattering.interference.type = "all";
	}
		
	scattering.probe = getstring(rootsetting["scattering"]["probe"]);
	
	// END OF scattering section //
	// START OF output section //
	
	if (rootsetting["output"].exists("prefix")) {
		output.prefix=getstring(rootsetting["output"]["prefix"]);
	}
	else {
		output.prefix="scattering-data-";
	}
	
	for (int i=0;i<rootsetting["output"].getLength();i++) {
		

		if (rootsetting["output"][i].getType()!=libconfig::Setting::TypeGroup) continue;
		string settings_name = getstring(rootsetting[i].getName());
		
		string otype   =getstring(rootsetting["output"][i]["type"]);
		string oformat =getstring(rootsetting["output"][i]["format"]);

		stringstream ffname; ffname << output.prefix << settings_name << "." << oformat;
		OutputFileParameters op(otype,oformat,get_filepath(ffname.str()));
		output.files.push_back(op);
	}
	
	if (output.files.size()==0) {
		Err::Inst()->write("Aborting. No output files defined.");
		throw;
	}

	
	// END OF output section //
	// START OF database section //

	if (rootsetting["database"]["names"].exists("pdb")) {
		for (int i=0;i<rootsetting["database"]["names"]["pdb"].getLength();i++) {	
			string label = getstring(rootsetting["database"]["names"]["pdb"][i].getName());
			string regexp = getstring(rootsetting["database"]["names"]["pdb"][i]);
			database.names.pdb.reg(label,regexp);
		}		
	}

	for (int i=0;i<rootsetting["database"]["masses"].getLength();i++) {	
		string label = getstring(rootsetting["database"]["masses"][i].getName());
		double mass = getdouble(rootsetting["database"]["masses"][i]);
		size_t ID = database.atomIDs.get(label);
		database.masses.reg(ID,mass);
	}		

	for (int i=0;i<rootsetting["database"]["sizes"].getLength();i++) {	
		
		string label = getstring(rootsetting["database"]["sizes"][i].getName());
		size_t sizes_function_type = getlong(rootsetting["database"]["sizes"][i][0]);
		
		vector<double> values;
		for(size_t k = 1; k < rootsetting["database"]["sizes"][i].getLength(); ++k)
		{
			values.push_back(getdouble(rootsetting["database"]["sizes"][i][k]));
		}					
		
		size_t ID = database.atomIDs.get(label);
		database.volumes.reg(ID,values,sizes_function_type);
	}	
	
	// determine which probe-type to load: 
	bool probe_found=false;
	for (int i=0;i<rootsetting["database"]["sfactors"].getLength();i++) {	
		if (getstring(rootsetting["database"]["sfactors"][i].getName())== scattering.probe) {
		
			for (int j=0;j<rootsetting["database"]["sfactors"][i].getLength();j++) {	
				
				string label = getstring(rootsetting["database"]["sfactors"][i][j].getName());
				size_t sfactors_function_type = getlong(rootsetting["database"]["sfactors"][i][j][0]);
				size_t ID = database.atomIDs.get(label);
				
				vector<double> values;
				for(size_t k = 1; k < rootsetting["database"]["sfactors"][i][j].getLength(); ++k)
				{
					values.push_back(getdouble(rootsetting["database"]["sfactors"][i][j][k]));
				}					

				database.sfactors.reg(ID,values,sfactors_function_type);
			}			
			
			probe_found = true;
			break;
		}
	}
	if (!probe_found) throw;
	
	
	// END OF database section //
	// START OF limits section //
	
	// some limits which are used internally: 
	// frame cache limit

	limits.framecache_max = 2;

	if (rootsetting.exists("limits")) {
		if (rootsetting["limits"].exists("framecache_max")) limits.framecache_max = getlong(rootsetting["limits"]["framecache_max"]);
	}

	delete pconf;
};

string Params::guessformat(string filename) {
	// do the best you can to guess the format
	return "conf";
}

bool Params::check() {
	// implement a sanity check for the parameters
}

void Params::init(std::string filename,std::string format) {
	if (format=="") {
		Info::Inst()->write(string("Detecting format of configuration file: ") + filename);
		format = guessformat(filename);
		Info::Inst()->write(string("Detected format: ") + format);
	}
	
	// make the directory of the main configuration file the root for all others
	if (boost::filesystem::path(filename).is_complete()) 
		config_rootpath = boost::filesystem::path(filename).parent_path().string();
	else 
		config_rootpath = ( boost::filesystem::initial_path() / boost::filesystem::path(filename).parent_path() ).string();		
	
	
	if (format=="conf") {
		read_conf(filename);
	}
	if (format=="xml") {
		Err::Inst()->write("XML Configuration not yet implemented");
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


// end of file
