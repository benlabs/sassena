#include "coor3d.hpp"

#include <complex>
#include <cmath>
#include <iostream>
#include <fstream>
#include <geometry.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/math/special_functions.hpp>
#include <vector>

#include <boost/random/mersenne_twister.hpp>	
#include <boost/random/uniform_on_sphere.hpp>	
#include <boost/random/variate_generator.hpp>	

using namespace std;

using namespace boost::math;
using namespace boost;

//using namespace boost::numeric::ublas::detail; 
	
int main(int argc, char** argv) {

	boost::mt19937 rng; // that's my random number generator

	boost::uniform_on_sphere<double> s(3); // that's my distribution

	boost::variate_generator<boost::mt19937&, boost::uniform_on_sphere<double> > mysphere(rng,s);

	vector<double> r = mysphere();
	
	for(vector<double>::iterator ri=r.begin();ri!=r.end();ri++) {
		cout << *ri << endl;		
	}
	
	return 0;
}