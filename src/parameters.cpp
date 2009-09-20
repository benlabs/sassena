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
#include "xml_interface.hpp"

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


////////////////////////////////////////////////////////////////////////////////////////////////////
// Params section
////////////////////////////////////////////////////////////////////////////////////////////////////


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
	throw;		
//	// configuration file objects are not referenced/destroyed. mind the memory leak if using 'a lot'.
//	libconfig::Config* pconf = new libconfig::Config;
//
//	// store the configuration line by line into an internal buffer, 
//	// this is for keeping history
//	ifstream conffile(filename.c_str());
//	string line;
//	while (getline(conffile,line)) {
//		carboncopy.push_back(line);
//	}
//	conffile.close();
//
//	// construct the libconfig interface
//	try { pconf->readFile(filename.c_str()); } 
//		catch (libconfig::ParseException &e) { 
//		Err::Inst()->write(string("Parse exception (line ")+to_s(e.getLine())+string("): ")+string(e.getError()));
//		delete pconf; 
//		throw;
//	}
//
//	libconfig::Setting&	rootsetting = pconf->getRoot();
//	// some flags , which we will need
//	bool flag_nogroups = false;
//	
//	// START OF sample section //	
//	
//	// now read the parameters
//	sample.structure.file   = get_filepath(getstring(rootsetting["sample"]["structure"]["file"]));
//	sample.structure.format = getstring(rootsetting["sample"]["structure"]["format"]);
//
//	// read in group definitions
//	if (rootsetting["sample"].exists("groups") && (rootsetting["sample"]["groups"].getLength()>0)) {
//		for(size_t i = 0; i < rootsetting["sample"]["groups"].getLength(); ++i)
//		{
//			libconfig::Setting& sgroup = rootsetting["sample"]["groups"][i];					
//			string gn = "";
//			if (sgroup.exists("name")) gn = getstring(sgroup["name"]); // groupname
//			string fn = get_filepath(getstring(sgroup["file"])); // filename
//			string ft = getstring(sgroup["type"]); // filetype
//			string sn = "";
//			if (sgroup.exists("select")) sn = getstring(sgroup["select"]); // selection field name
//			double sv = 0.0;
//			if (sgroup.exists("select_value")) sv = getdouble(sgroup["select_value"]); // selection field value, positive match								
//			sample.groups[gn] = SampleGroupParameters(gn,fn,ft,sn,sv);
//		}
//	} else {
//		// no group definition
//		Warn::Inst()->write("No groups defined. Switching to full system scattering");
//		flag_nogroups = true;
//	}	
//
//	// apply deuteration
//	if (rootsetting["sample"].exists("deuter")) {
//		if (rootsetting["sample"]["deuter"].getType()==libconfig::Setting::TypeString) {
//			Info::Inst()->write(string("Deuteration of group ") + getstring(rootsetting["sample"]["deuter"]));		
//			sample.deuter.push_back(rootsetting["sample"]["deuter"]);
//		}
//		else if (rootsetting["sample"]["deuter"].getType()==libconfig::Setting::TypeList) {
//	    	for (int i=0;i<rootsetting["sample"]["deuter"].getLength();i++) {
//				Info::Inst()->write(string("Deuteration of group ") + getstring(rootsetting["sample"]["deuter"][i]));		
//				sample.deuter.push_back(rootsetting["sample"]["deuter"][i]);				
//			}
//		}
//	}
//	else {
//		Info::Inst()->write("No explicit deuteration");		
//	}
//	
//	
//	// read in frame information
//	size_t def_first=0;	size_t def_last=0;	bool def_last_set=false; size_t def_stride = 1;
//	if (rootsetting["sample"]["frames"].exists("first"))   def_first  = getlong(rootsetting["sample"]["frames"]["first"]);
//	if (rootsetting["sample"]["frames"].exists("last"))  { def_last   = getlong(rootsetting["sample"]["frames"]["last"]); def_last_set = true; }
//	if (rootsetting["sample"]["frames"].exists("stride"))  def_stride = getlong(rootsetting["sample"]["frames"]["stride"]);
//	
//	for (int i=0;i<rootsetting["sample"]["frames"].getLength();i++) {
//		if (rootsetting["sample"]["frames"][i].getType()!=libconfig::Setting::TypeGroup) continue;
//		SampleFramesetParameters fset;
//		fset.first = def_first;	fset.last = def_last; fset.last_set = def_last_set; fset.stride = def_stride;
//		fset.filename = get_filepath(getstring(rootsetting["sample"]["frames"][i]["file"]));
//		fset.type = getstring(rootsetting["sample"]["frames"][i]["type"]);
//		if (rootsetting["sample"]["frames"][i].exists("first"))   fset.first  = getlong(rootsetting["sample"]["frames"][i]["first"]);
//		if (rootsetting["sample"]["frames"][i].exists("last"))  { fset.last   = getlong(rootsetting["sample"]["frames"][i]["last"]); fset.last_set = true; }
//		if (rootsetting["sample"]["frames"][i].exists("stride"))  fset.stride = getlong(rootsetting["sample"]["frames"][i]["stride"]);
//		
//		sample.frames.push_back(fset);
//		Info::Inst()->write(string("Added frames from ")+fset.filename+string(" using format: ")+fset.type);				
//	}
//	
//	// periodic boundary behavior and/or postprocessing
//	
//	if (rootsetting["sample"]["pbc"]["wrapping"]) {
//		sample.pbc.wrapping = true;
//		sample.pbc.center = "";
//		if (rootsetting["sample"]["pbc"].exists("center")) sample.pbc.center = getstring(rootsetting["sample"]["pbc"]["center"]);
//		
//		Info::Inst()->write(string("Turned wrapping ON with center group ")+sample.pbc.center);
//	}
//	else {
//		sample.pbc.wrapping = false;		
//		Info::Inst()->write("Turned wrapping OFF");
//	}	
//	
//	// END OF sample section //
//	// START OF scattering section //
//
//	if (rootsetting["scattering"].exists("background")) {
//		string m = getstring(rootsetting["scattering"]["background"]["type"]);
//		// rescale automatically triggers background analysis routine. Maybe independent in the future...
//		scattering.background.type = m;
//
//		if (rootsetting["scattering"]["background"].exists("factor")) {
//			scattering.background.factor = getdouble(rootsetting["scattering"]["background"]["factor"]);
//		}
//
//		for(size_t i = 0; i < rootsetting["scattering"]["background"]["phases"].getLength(); ++i)
//		{
//			libconfig::Setting& phase = rootsetting["scattering"]["background"]["phases"][i];
//			if (phase.getType()!=libconfig::Setting::TypeGroup) continue;
//			ScatteringBackgroundPhaseParameters phasestruct;
//			phasestruct.selection = getstring(phase["selection"]);
//			phasestruct.scaling = getstring(phase["scaling"]);
//			phasestruct.factor = getstring(phase["factor"]);
//			phasestruct.nullrange = 0.0;
//			if (phase.exists("nullrange")) phasestruct.nullrange = getstring(phase["nullrange"]);
//			scattering.background.phases.push_back(phasestruct);
//		}
//		
//		scattering.background.gridspacing = 3;
//		if (rootsetting["scattering"]["background"].exists("gridspacing")) {
//			scattering.background.gridspacing = getdouble(rootsetting["scattering"]["background"]["gridspacing"]);
//		}
//
//		scattering.background.gridspacing = 3;
//		if (rootsetting["scattering"]["background"].exists("gridspacing")) {
//			scattering.background.gridspacing = getdouble(rootsetting["scattering"]["background"]["gridspacing"]);
//		}
//	}
//
//	// generating qqqvectors, i.e. the spectrum
//	string vm = getstring(rootsetting["scattering"]["vectors"]["method"]);
//	if (vm=="single") {	
//		libconfig::Setting& s = rootsetting["scattering"]["vectors"]["single"];
//		double x = getdouble(s["direction"][0]);
//		double y = getdouble(s["direction"][1]);
//		double z = getdouble(s["direction"][2]);		
//		CartesianCoor3D direction(x,y,z);
//		
//		double scal = getdouble(s["scale"]);				
//		scattering.qqqvectors.push_back(scal*direction);
//	}
//	else if (vm=="linear") {	
//		libconfig::Setting& s = rootsetting["scattering"]["vectors"]["linear"];
//		double x = getdouble(s["direction"][0]);
//		double y = getdouble(s["direction"][1]);
//		double z = getdouble(s["direction"][2]);		
//		CartesianCoor3D direction(x,y,z);
//				
//		double from = getdouble(s["from"]  );		
//		double to =   getdouble(s["to"]    );				
//		int points =  getlong(s["points"]);
//		double exponent = 1.0;
//		if (s.exists("exponent")) exponent = getdouble(s["exponent"]);
//		
//		for (int i=0;i<points;i++) {		
//			double scal = from + powf((i*1.0/(points-1)),exponent)*(to-from);
//			scattering.qqqvectors.push_back(scal*direction);
//		}
//	}
//	else if (vm=="map") {
//		libconfig::Setting& setting = rootsetting["scattering"]["vectors"]["map"];
//		libconfig::Setting& s1 = setting[0];
//		libconfig::Setting& s2 = setting[1];
//		
//		int p1 = getlong(s1["points"]);
//		int p2 = getlong(s2["points"]);
//		double from1 = getdouble(s1["from"]);
//		double from2 = getdouble(s2["from"]);
//		double to1 = getdouble(s1["to"]);
//		double to2 = getdouble(s2["to"]);
//		double exponent1 = 1.0;
//		double exponent2 = 1.0;		
//		if (s1.exists("exponent")) exponent1 = getdouble(s1["exponent"]);
//		if (s2.exists("exponent")) exponent2 = getdouble(s2["exponent"]);
//		
//		double x1 = getdouble(s1["direction"][0]);
//		double y1 = getdouble(s1["direction"][1]);
//		double z1 = getdouble(s1["direction"][2]);		
//		CartesianCoor3D d1(x1,y1,z1);
//		double x2 = getdouble(s2["direction"][0]);
//		double y2 = getdouble(s2["direction"][1]);
//		double z2 = getdouble(s2["direction"][2]);		
//		CartesianCoor3D d2(x2,y2,z2);
//		
//		for (int i=0;i<p1;i++) {		
//			double scal1 = from1 + powf((i*1.0/(p1-1)),exponent1)*(to1-from1);
//			for (int j=0;j<p2;j++) {			
//				double scal2 = from2 + powf((j*1.0/(p2-1)),exponent2)*(to2-from2);
//				scattering.qqqvectors.push_back(scal1*d1+scal2*d2);
//			}
//		}
//	}
//	else if (vm=="cube") {
//		libconfig::Setting& setting = rootsetting["scattering"]["vectors"]["cube"];
//		libconfig::Setting& s1 = setting[0];
//		libconfig::Setting& s2 = setting[1];
//		libconfig::Setting& s3 = setting[2];
//		
//		int p1 = getlong(s1["points"]);
//		int p2 = getlong(s2["points"]);
//		int p3 = getlong(s3["points"]);		
//		double from1 = getdouble(s1["from"]);
//		double from2 = getdouble(s2["from"]);
//		double from3 = getdouble(s3["from"]);		
//		double to1 = getdouble(s1["to"]);
//		double to2 = getdouble(s2["to"]);
//		double to3 = getdouble(s3["to"]);
//		double exponent1 = 1.0;
//		double exponent2 = 1.0;		
//		double exponent3 = 1.0;				
//		if (s1.exists("exponent")) exponent1 = getdouble(s1["exponent"]);
//		if (s2.exists("exponent")) exponent2 = getdouble(s2["exponent"]);
//		if (s3.exists("exponent")) exponent3 = getdouble(s3["exponent"]);
//				
//		double x1 = getdouble(s1["direction"][0]);
//		double y1 = getdouble(s1["direction"][1]);
//		double z1 = getdouble(s1["direction"][2]);		
//		CartesianCoor3D d1(x1,y1,z1);
//		double x2 = getdouble(s2["direction"][0]);
//		double y2 = getdouble(s2["direction"][1]);
//		double z2 = getdouble(s2["direction"][2]);		
//		CartesianCoor3D d2(x2,y2,z2);
//		double x3 = getdouble(s3["direction"][0]);
//		double y3 = getdouble(s3["direction"][1]);
//		double z3 = getdouble(s3["direction"][2]);		
//		CartesianCoor3D d3(x3,y3,z3);
//				
//		for (int i=0;i<p1;i++) {		
//			double scal1 = from1 + powf((i*1.0/(p1-1)),exponent1)*(to1-from1);
//			for (int j=0;j<p2;j++) {			
//				double scal2 = from2 + powf((j*1.0/(p2-1)),exponent2)*(to2-from2);				
//				for (int k=0;k<p3;k++) {			
//					double scal3 = from3 + powf((k*1.0/(p3-1)),exponent3)*(to3-from3);									
//					scattering.qqqvectors.push_back(scal1*d1+scal2*d2+scal3*d3);
//				}
//			}
//		}
//	}	
//	else if (vm=="file") {
//		string qqqfilename = get_filepath(getstring(rootsetting["scattering"]["vectors"]["file"]));
//		ifstream qqqfile(qqqfilename.c_str());
//		
//		double x,y,z; 
//		while (qqqfile >> x >> y >> z) {
//			scattering.qqqvectors.push_back(CartesianCoor3D(x,y,z));
//		}
//	}
//
//	scattering.framestride = 1;
//	if (rootsetting["scattering"].exists("framestride")) scattering.framestride = getlong(rootsetting["scattering"]["framestride"]);
//
//	if (rootsetting["scattering"].exists("correlation")) {
//		scattering.correlation.type = getstring(rootsetting["scattering"]["correlation"]["type"]);
//	}
//	else {
//		scattering.correlation.type = "none";
//	}
//
//	scattering.target = getstring(rootsetting["scattering"]["target"]);
//	
//	if (rootsetting["scattering"].exists("average")) {
//		scattering.average.type   = getstring(rootsetting["scattering"]["average"]["type"]);
//		if (scattering.average.type!="none") {
//			scattering.average.method = getstring(rootsetting["scattering"]["average"]["method"]);
//			if (scattering.average.method=="bruteforce") {		
//				scattering.average.vectors = getstring(rootsetting["scattering"]["average"]["vectors"]);
//			}
//			scattering.average.resolution =  getdouble(rootsetting["scattering"]["average"]["resolution"]);
//		}
//	}
//	else {
//		scattering.average.type = "none";
//	}
//	
//	if (rootsetting["scattering"].exists("interference")) {
//		scattering.interference.type = getstring(rootsetting["scattering"]["interference"]["type"]);
//	}
//	else {
//		scattering.interference.type = "all";
//	}
//		
//	scattering.probe = getstring(rootsetting["scattering"]["probe"]);
//	
//	// END OF scattering section //
//	// START OF output section //
//	
//	if (rootsetting["output"].exists("prefix")) {
//		output.prefix=getstring(rootsetting["output"]["prefix"]);
//	}
//	else {
//		output.prefix="scattering-data-";
//	}
//	
//	for (int i=0;i<rootsetting["output"].getLength();i++) {
//		
//
//		if (rootsetting["output"][i].getType()!=libconfig::Setting::TypeGroup) continue;
//		string settings_name = getstring(rootsetting[i].getName());
//		
//		string otype   =getstring(rootsetting["output"][i]["type"]);
//		string oformat =getstring(rootsetting["output"][i]["format"]);
//
//		stringstream ffname; ffname << output.prefix << settings_name << "." << oformat;
//		OutputFileParameters op(otype,oformat,get_filepath(ffname.str()));
//		output.files.push_back(op);
//	}
//	
//	if (output.files.size()==0) {
//		Err::Inst()->write("Aborting. No output files defined.");
//		throw;
//	}
//
//	
//	// END OF output section //
//	// START OF limits section //
//	
//	// some limits which are used internally: 
//	// frame cache limit
//
//	limits.framecache_max = 2;
//	limits.static_load_imbalance_max = 0.05; // default is 5% loss due to bad partioning.
//	limits.buffers.allgather_max = 200*1000*1000; // don't use more than: 200 Megabyte for this buffer! calculated: nn*(nf/nn)*2 , e.g. 200*(1000010/200)*2 * 8 byte= 160 Mbyte
//
//	if (rootsetting.exists("limits")) {
//		if (rootsetting["limits"].exists("framecache_max")) limits.framecache_max = getlong(rootsetting["limits"]["framecache_max"]);
//		if (rootsetting["limits"].exists("static_load_imbalance_max")) limits.static_load_imbalance_max = getdouble(rootsetting["limits"]["static_load_imbalance_max"]);		
//		if (rootsetting.exists("buffers")) {
//			if (rootsetting["limits"].exists("framecache_max")) limits.framecache_max = getlong(rootsetting["limits"]["framecache_max"]);
//			if (rootsetting["limits"].exists("static_load_imbalance_max")) limits.static_load_imbalance_max = getdouble(rootsetting["limits"]["static_load_imbalance_max"]);		
//		}
//		
//	}
//	// END OF limits section //
//	// START OF debug section //
//	
//	debug.timer = false; // this adds a log message when a timer is started/stopped
//	debug.barriers = false; // this de-/activates collective barriers before each collective operation, this way all nodes are synchronized before the communication takes place. This is an important step towards analysis of timing.
//	debug.scatter_from_frame = false;
//	
//	if (rootsetting.exists("debug")) {
//		if (rootsetting["debug"].exists("timer")) debug.timer = getbool(rootsetting["debug"]["timer"]);
//		if (rootsetting["debug"].exists("barriers")) debug.barriers = getbool(rootsetting["debug"]["barriers"]);		
//		if (rootsetting["debug"].exists("scatter_from_frame")) debug.scatter_from_frame = getbool(rootsetting["debug"]["scatter_from_frame"]);				
//	}
//
//	delete pconf;
};

void Params::read_xml(std::string filename) {
		
	// store the configuration line by line into an internal buffer, 
	// this is for keeping history
	ifstream conffile(filename.c_str());
	string line;
	while (getline(conffile,line)) {
		carboncopy.push_back(line);
	}
	conffile.close();
		
	// START OF sample section //	
	XMLInterface xmli(filename);

	// now read the parameters
	
	sample.structure.file   = get_filepath(xmli.get_value<std::string>("//sample/structure/file"));
	sample.structure.format = xmli.get_value<std::string>("//sample/structure/format");
	vector<XMLElement> selections = xmli.get("//sample/selections/selection");
	for(size_t i = 0; i < selections.size(); ++i)
	{
		xmli.set_current(selections[i]);
		string sname, fn, ff, sn;
		double sv;
		if (xmli.exists("./name")) {
			sname = xmli.get_value<string>("./name") ;
		}
		if (xmli.exists("./file")) {		
			fn = get_filepath( xmli.get_value<string>("./file") );
		}
		if (xmli.exists("./format")) {		
			ff = xmli.get_value<string>("./format") ;
		}
		if (xmli.exists("./select")) {		
			sn = xmli.get_value<string>("./select") ;
		}
		if (xmli.exists("./select_value")) {
			sv = xmli.get_value<double>("./select_value") ;
		}
		sample.groups[sname] = SampleGroupParameters(sname,fn,ff,sn,sv);
	}

	// apply deuteration
	if (xmli.exists("//sample/deuterations")) {
		vector<XMLElement> deuterations = xmli.get("//sample/deuterations/deuteration");
		for(size_t i = 0; i < deuterations.size(); ++i)
		{
			xmli.set_current(deuterations[i]);
			string deuter = xmli.get_value<string>(".");
			Info::Inst()->write(string("Deuteration of group ") + deuter);	
			sample.deuter.push_back(deuter);
		}
	}
	else {
		Info::Inst()->write("No explicit deuteration");		
	}
	
	// END OF sample section //	

	// START OF scattering section //	

	// END OF scattering section //	

	// START OF output section //	

	// END OF output section //	

	// START OF limits section //	

	// END OF limits section //	

	// START OF debug section //	

	// END OF debug section //	
	
	// START OF sample section //	

	// read in frame information
	size_t def_first=0;	size_t def_last=0;	bool def_last_set=false; size_t def_stride = 1;
	if (xmli.exists("//sample/frames/first"))   def_first  = xmli.get_value<size_t>("//sample/frames/first");
	if (xmli.exists("//sample/frames/last"))  { def_last   = xmli.get_value<size_t>("//sample/frames/last"); def_last_set = true; }
	if (xmli.exists("//sample/frames/stride"))  def_stride = xmli.get_value<size_t>("//sample/frames/stride");
	
	vector<XMLElement> framesets = xmli.get("//sample/frames/frameset");
	for(size_t i = 0; i < framesets.size(); ++i)
	{
		xmli.set_current(framesets[i]);
		SampleFramesetParameters fset;	
		fset.first = def_first;	fset.last = def_last; fset.last_set = def_last_set; fset.stride = def_stride;
		fset.filename = get_filepath(xmli.get_value<string>("./file"));
		fset.type = xmli.get_value<string>("./format");
		if (xmli.exists("./first"))   fset.first  = xmli.get_value<size_t>("./first");
		if (xmli.exists("./last"))  { fset.last   = xmli.get_value<size_t>("./last"); fset.last_set = true; }
		if (xmli.exists("./stride"))  fset.stride = xmli.get_value<size_t>("./stride");
		
		sample.frames.push_back(fset);
		Info::Inst()->write(string("Added frames from ")+fset.filename+string(" using format: ")+fset.type);
		Info::Inst()->write(string("Options: first=")+to_s(fset.first)+string(", last=")+to_s(fset.last)+string(", lastset=")+to_s(fset.last_set)+string(", stride=")+to_s(fset.stride));		
	}
		
	if (xmli.exists("//sample/motions")) {
	
		vector<XMLElement> motions = xmli.get("//sample/motions/motion");
		for(size_t i = 0; i < motions.size(); ++i)
		{
			xmli.set_current(motions[i]);
			SampleMotionParameters motion;	
			motion.type = "linear";	
			motion.displace = 0.0; 
			motion.direction=CartesianCoor3D(1,0,0);
			motion.frequency=2*M_PI/1000.0; // corresponds to one full cycle per 1000 frames, used for linear oscillation and rotation
			if (xmli.exists("./type"))   motion.type  = xmli.get_value<string>("./type");
			if (xmli.exists("./displace"))  motion.displace   = xmli.get_value<double>("./displace");
			if (xmli.exists("./frequency"))  motion.frequency   = xmli.get_value<double>("./frequency");			
			if (xmli.exists("./direction")) {
				motion.direction.x   = xmli.get_value<double>("./direction/x");
				motion.direction.y   = xmli.get_value<double>("./direction/y");
				motion.direction.z   = xmli.get_value<double>("./direction/z");				
			} 

			sample.motions.push_back(motion);
			Info::Inst()->write(string("Adding additional motion to sample: type=")+motion.type+string(", displacement=")+to_s(motion.displace));
		}
	}	
		
	// periodic boundary behavior and/or postprocessing
	sample.pbc.wrapping = false;	
	sample.pbc.center = "";
	if (xmli.exists("//sample/pbc")) {
		if (xmli.exists("//sample/pbc/wrapping")) sample.pbc.wrapping = xmli.get_value<bool>("//sample/pbc/wrapping");
		if (xmli.exists("//sample/pbc/center"))   sample.pbc.center   = xmli.get_value<string>("//sample/pbc/center");
		
		if (sample.pbc.wrapping) {
			Info::Inst()->write(string("Turned wrapping ON with center group ")+sample.pbc.center);
		} else {
			Info::Inst()->write("Turned wrapping OFF");
		}
			
	} else {
		Info::Inst()->write("No periodic boundary treatment");
	}
	
	// END OF sample section //
	// START OF scattering section //

	if (xmli.exists("//scattering/background")) {
		scattering.background.type = xmli.get_value<string>("//scattering/background/type");
		if (xmli.exists("//scattering/background/factor")) {
			scattering.background.factor = xmli.get_value<double>("//scattering/background/factor");			
		}

		if (xmli.exists("//scattering/background/phases")) {
			vector<XMLElement> phases = xmli.get("//scattering/background/phases");
				
			for(size_t i = 0; i < phases.size(); ++i)
			{
				xmli.set_current(phases[i]);
				ScatteringBackgroundPhaseParameters struc;
				struc.selection = xmli.get_value<string>("./selection");					
				if (xmli.exists("./scaling")) struc.scaling = xmli.get_value<string>("./scaling");	
				if (xmli.exists("./factor")) struc.factor = xmli.get_value<double>("./factor");		
				scattering.background.phases.push_back(struc);
			}
		}
		if (xmli.exists("//scattering/background/gridspacing")) {
			scattering.background.gridspacing = xmli.get_value<double>("//scattering/background/gridspacing");  				
		}
		if (xmli.exists("//scattering/background/stride")) {
			scattering.background.stride  = xmli.get_value<double>("//scattering/background/stride"); 				
		}		
	}	
	// generating qqqvectors, i.e. the spectrum

	string vt = xmli.get_value<string>("//scattering/vectors/type");
	if (vt=="single") {	

		double x = xmli.get_value<double>("//scattering/vectors/single/x");
		double y = xmli.get_value<double>("//scattering/vectors/single/y");
		double z = xmli.get_value<double>("//scattering/vectors/single/z");

		scattering.qvectors.push_back(CartesianCoor3D(x,y,z));
	}
	else if (vt=="scan") {	
		
		vector<XMLElement> scans = xmli.get("//scattering/vectors/scan");

		for(size_t i = 0; i < scans.size(); ++i)
		{
			xmli.set_current(scans[i]);
			ScatteringVectorsScanParameters sc;
			
			sc.basevector.x = xmli.get_value<double>("./basevector/x");
			sc.basevector.y = xmli.get_value<double>("./basevector/y");
			sc.basevector.z = xmli.get_value<double>("./basevector/z");

			sc.from     = xmli.get_value<double>("./from");
			sc.to       = xmli.get_value<double>("./to");
			sc.points   = xmli.get_value<size_t>("./points");
			if (xmli.exists("./exponent")) sc.exponent = xmli.get_value<double>("./exponent");
			scattering.qvectors.scans.push_back(sc);
		}
		
		scattering.qvectors.create_from_scans();
	}
	else if (vt=="file") {
		string qqqfilename = get_filepath(xmli.get_value<string>("//scattering/vectors/file"));
		ifstream qqqfile(qqqfilename.c_str());
		
		double x,y,z; 
		while (qqqfile >> x >> y >> z) {
			scattering.qvectors.push_back(CartesianCoor3D(x,y,z));
		}
	}
	
	if (xmli.exists("//scattering/correlation")) {
		if (xmli.exists("//scattering/correlation/type")) {
			scattering.correlation.type = xmli.get_value<string>("//scattering/correlation/type");
		}
		if (xmli.exists("//scattering/correlation/method")) {
			scattering.correlation.method = xmli.get_value<string>("//scattering/correlation/method");
		}
	}

	if (xmli.exists("//scattering/average")) {
		if (xmli.exists("//scattering/average/orientation")) {
			if (xmli.exists("//scattering/average/orientation/type")) { // sphere cylinder none
				scattering.average.orientation.type = xmli.get_value<string>("//scattering/average/orientation/type");
			}
			if (xmli.exists("//scattering/average/orientation/method")) { // bruteforce multipole debye
				scattering.average.orientation.method = xmli.get_value<string>("//scattering/average/orientation/method");
			}
			if (xmli.exists("//scattering/average/orientation/vectors")) { // mcboostuniform linearraster...
				scattering.average.orientation.vectors = xmli.get_value<string>("//scattering/average/orientation/vectors");
			}
			if (xmli.exists("//scattering/average/orientation/resolution")) { // count vectors ... , or order for multipole...
				scattering.average.orientation.resolution = xmli.get_value<double>("//scattering/average/orientation/resolution");
			}
			if (xmli.exists("//scattering/average/orientation/axis")) { // necessary for cylinder
				scattering.average.orientation.axis.x = xmli.get_value<double>("//scattering/average/orientation/axis/x");
				scattering.average.orientation.axis.y = xmli.get_value<double>("//scattering/average/orientation/axis/y");
				scattering.average.orientation.axis.z = xmli.get_value<double>("//scattering/average/orientation/axis/z");
			}
			if (xmli.exists("//scattering/average/orientation/origin")) { // necessary for cylinder
				scattering.average.orientation.origin = xmli.get_value<string>("//scattering/average/orientation/origin");
			}
		}
		if (xmli.exists("//scattering/average/motion")) {
			if (xmli.exists("//scattering/average/motion/type")) { // sphere cylinder none
				scattering.average.motion.type = xmli.get_value<string>("//scattering/average/motion/type");
			}
			if (xmli.exists("//scattering/average/motion/method")) { // bruteforce multipole debye
				scattering.average.motion.method = xmli.get_value<string>("//scattering/average/motion/method");
			}
			if (xmli.exists("//scattering/average/motion/vectors")) { // mcboostuniform linearraster...
				scattering.average.motion.vectors = xmli.get_value<string>("//scattering/average/motion/vectors");
			}
			if (xmli.exists("//scattering/average/motion/resolution")) { // count vectors ... , or order for multipole...
				scattering.average.motion.resolution = xmli.get_value<double>("//scattering/average/motion/resolution");
			}
			if (xmli.exists("//scattering/average/motion/axis")) { // necessary for cylinder
				scattering.average.motion.axis.x = xmli.get_value<double>("//scattering/average/motion/axis/x");
				scattering.average.motion.axis.y = xmli.get_value<double>("//scattering/average/motion/axis/y");
				scattering.average.motion.axis.z = xmli.get_value<double>("//scattering/average/motion/axis/z");
			}
			if (xmli.exists("//scattering/average/motion/origin")) { // necessary for cylinder
				scattering.average.motion.origin = xmli.get_value<string>("//scattering/average/motion/origin");
			}
		}
		
	}
	
	if (xmli.exists("//scattering/interference")) {
		if (xmli.exists("//scattering/interference/type")) {
			scattering.interference.type = xmli.get_value<string>("//scattering/interference/type");
		}
	}
		
	if (xmli.exists("//scattering/target")) {
		scattering.target = xmli.get_value<string>("//scattering/target");
	} else {
		Err::Inst()->write("You have to specify a target (any valid selection)");
		throw;
	}

	if (xmli.exists("//scattering/scatterfactors")) {
		scattering.scatterfactors = xmli.get_value<string>("//scattering/scatterfactors");
	} else {
		Err::Inst()->write("You have to specify a scatterfactor type (e.g. neutron-coherent)");
		throw;
	}
	
	// END OF scattering section //
	// START OF output section //
	
	if (!xmli.exists("//output")) {
		Err::Inst()->write("You have to specify an output section");
		throw;
	}
	
	if (xmli.exists("//output/prefix")) {
		output.prefix = xmli.get_value<string>("//output/prefix");
	}
	
	vector<XMLElement> outputfiles = xmli.get("//output/file");
	
	for(size_t i = 0; i < outputfiles.size(); ++i)
	{
		xmli.set_current(outputfiles[i]);
		OutputFileParameters ofp;
		ofp.name = xmli.get_value<string>("./name");
		ofp.method = xmli.get_value<string>("./method");
		ofp.format = xmli.get_value<string>("./format");
		ofp.filename = get_filepath(output.prefix + ofp.name + string(".") + ofp.format  );
		output.files.push_back(ofp);
	}
	
	if (output.files.size()==0) {
		Err::Inst()->write("Aborting. No output files defined.");
		throw;
	}

	// END OF output section //
	// START OF limits section //
	
	// some limits which are used internally: 
	// frame cache limit

	limits.framecache_max = 1;
	limits.coordinatesets_cache_max = 0;
	limits.static_load_imbalance_max = 0.05; // default is 5% loss due to bad partioning.
	limits.buffers.allgather_max = 200*1000*1000; // don't use more than: 200 Megabyte for this buffer! calculated: nn*(nf/nn)*2 , e.g. 200*(1000010/200)*2 * 8 byte= 160 Mbyte

	if (xmli.exists("//limits")) {
		if (xmli.exists("//limits/framecache_max")) {
			limits.framecache_max = xmli.get_value<size_t>("//limits/framecache_max");
		}
		if (xmli.exists("//limits/coordinatesets_cache_max")) {
			limits.coordinatesets_cache_max = xmli.get_value<size_t>("//limits/coordinatesets_cache_max");
		}		
		if (xmli.exists("//limits/static_load_imbalance_max")) {
			limits.static_load_imbalance_max = xmli.get_value<double>("//limits/static_load_imbalance_max");
		}
		
	}

	// END OF limits section //
	// START OF debug section //

	
	debug.timer = false; // this adds a log message when a timer is started/stopped
	debug.barriers = false; // this de-/activates collective barriers before each collective operation, this way all nodes are synchronized before the communication takes place. This is an important step towards analysis of timing.
	debug.scatter_from_frame = false;
	debug.decomposition_split=0; // 0 stands for automatic
	if (xmli.exists("//debug")) {
		if (xmli.exists("//debug/timer")) {
			debug.timer = xmli.get_value<bool>("//debug/timer");
		}
		if (xmli.exists("//debug/barriers")) {
			debug.barriers = xmli.get_value<bool>("//debug/barriers");
		}
		if (xmli.exists("//debug/scatter_from_frame")) {
			debug.scatter_from_frame = xmli.get_value<bool>("//debug/scatter_from_frame");
		}	
		if (xmli.exists("//debug/decomposition_split")) {
			debug.decomposition_split = xmli.get_value<size_t>("//debug/decomposition_split");
		}		
	}
	
};

string Params::guessformat(string filename) {
	// do the best you can to guess the format
	// guess by filename
	boost::regex e_conf(".*\\.conf$");
	boost::regex e_xml(".*\\.xml$");
	
	if (boost::regex_match(filename,e_conf)) return "conf";
	if (boost::regex_match(filename,e_xml)) return "xml";
	
	// else: problem
	Err::Inst()->write("Parameter file format could not be detected (ending with: conf , xml ?)");
	throw;
}

bool Params::check() {
	return true;
	// implement a sanity check for the parameters
}

void Params::init(std::string filename) {
	Info::Inst()->write(string("Detecting format of configuration file: ") + filename);
	string format = guessformat(filename);
	Info::Inst()->write(string("Detected format: ") + format);
	
	// make the directory of the main configuration file the root for all others
	if (boost::filesystem::path(filename).is_complete()) 
		config_rootpath = boost::filesystem::path(filename).parent_path().string();
	else 
		config_rootpath = ( boost::filesystem::initial_path() / boost::filesystem::path(filename).parent_path() ).string();		

	if (format=="conf") {
		read_conf(filename);
	}
	if (format=="xml") {
		read_xml(filename);
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

void Params::write(std::string filename, std::string format) {
	if (format=="") {
		Info::Inst()->write("No parameters format specified. Will use conf.");
		format = "conf";
	}
	
	if (format=="conf") {
		Err::Inst()->write("conf parameters format (write) not yet implemented.");
		throw;
	}
	if (format=="xml") {
		Err::Inst()->write("XML parameters format (write) not yet implemented.");
		throw;
	}
}

void ScatteringVectorsParameters::create_from_scans() {
	// read out the scans and push the results onto the internal vector

	if (scans.size()>3) {
		Err::Inst()->write("More than 3 scan definitions are not supported.");
		throw;
	}
	
	// local vector unfolds each scan first, then we'll do element-wise vector addition
	vector< vector<CartesianCoor3D> > qvectors(scans.size());
	for(size_t i = 0; i < scans.size(); ++i)
	{
		for (int j=0;j<scans[i].points;j++) {		
			double scal = scans[i].from + powf((j*1.0/(scans[i].points-1)),scans[i].exponent)*(scans[i].to-scans[i].from);
			qvectors[i].push_back(scal*scans[i].basevector);
		}		
	}

	// trivial case: only one scan!
	if (scans.size()==1) {
		for(size_t i = 0; i < qvectors[0].size(); ++i)
		{
			this->push_back(qvectors[0][i]);
		}
		return;
	}
	// if scans.size()>1 , not trivial..
	
	if (scans.size()==2) {
		for(size_t i = 0; i < qvectors[0].size(); ++i)
		{
			for(size_t j = 0; j < qvectors[1].size(); ++j)
			{
				this->push_back(qvectors[0][i]+qvectors[1][j]);
			}
		}	
	}
	
	if (scans.size()==3) {
		for(size_t i = 0; i < qvectors[0].size(); ++i)
		{
			for(size_t j = 0; j < qvectors[1].size(); ++j)
			{
				for(size_t k = 0; k < qvectors[2].size(); ++k)
				{
					this->push_back(qvectors[0][i]+qvectors[1][j]+qvectors[2][k]);
				}
			}
		}	
	}	
}

// end of file
