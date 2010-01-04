#include "common.hpp"

std::vector<double> flatten(std::vector<std::complex<double> >& cvalues) {
	std::vector<double> rvalues(cvalues.size()*2);
	for(size_t i = 0; i < cvalues.size(); ++i)
	{
		rvalues[2*i]=cvalues[i].real();
		rvalues[2*i+1]=cvalues[i].imag();
	}
	return rvalues;
}

std::vector<std::complex<double> > compress(std::vector<double>& rvalues) {
	std::vector<std::complex<double> > cvalues(rvalues.size()/2);
	for(size_t i = 0; i < cvalues.size(); ++i)
	{
		cvalues[i]=std::complex<double>(rvalues[2*i],rvalues[2*i+1]);
	}
	return cvalues;
}

// end of file
