#include <iostream>
#include <libconfig.h++>

using namespace libconfig;
using namespace std;

int main(int argc,char** argv) {
	
	Config conf;
	
	conf.readFile(argv[1]);
	
	Setting& setting = conf.getRoot();
	
	int len = setting.getLength();
	for (int i=0;i<len;i++) {
		string settings_name = setting[i].getName();
		cout << settings_name << endl;
	}
	
	return 0;
}

// end of file
