#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>

#include <coor3d.hpp>

using namespace std;

int main(int argc,char** argv) {
        const int resX =101;
        const int resY =101;
        double valmap[resX][resY];
        int histo[255];
        int i,j;

        // format : <q> <sumA>  <sumB>  <fsquared>                                                                                                             
        ifstream ifile(argv[1]);

        for (i=0;i<resX;i++) {                                                                                                                         
                for (j=0;j<resY;j++) {
                        valmap[i][j] = 0;
                }
        }

		typedef vector<pair<pair<CartesianCoor3D,int>, double> > scatvec;
		scatvec scat;

		clog << "INFO>> " << "Reading values from file " << argv[1] << endl;
		string line;
		while (true) {
			pair<pair<CartesianCoor3D,int>,double> keyval;
	        float qx,qy,qz,fsquared;
			int frame;
			ifile >> qx >> qy >> qz >> frame >> fsquared;
			if ( ifile.eof() ) break;			
						
			keyval.first.first = CartesianCoor3D(qx,qy,qz);
			keyval.first.second = frame;
			keyval.second = fsquared;
			scat.push_back(keyval);
		}

		// get limits first
		CartesianCoor3D firstcoor = (scat.begin()->first).first;
		double left = firstcoor.x;
		double right = firstcoor.x;		
		double top = firstcoor.y;		
		double bottom = firstcoor.y;	

		clog << "INFO>> " << "Scanning limits..., centering " << endl;
			
		for(scatvec::iterator si=scat.begin();si!=scat.end();si++) {
			CartesianCoor3D c = (si->first).first;
			left = c.x < left ?  c.x : left;
			right = c.x > right ?  c.x : right;
			top = c.x > top ?  c.y : top;
			bottom = c.x < bottom ?  c.y : bottom;
		}
		clog << "INFO>> " << "limit left "   << left << endl;
		clog << "INFO>> " << "limit right "  << right << endl;
		clog << "INFO>> " << "limit top "    << top << endl;
		clog << "INFO>> " << "limit bottom " << bottom << endl;

		
		clog << "INFO>> " << "Values to write: " << scat.size() << endl;
		
		double fsmax=0;
		double rescorX = 0.5 * (right-left) / resX;
		double rescorY = 0.5 * (top-bottom) / resY;		
		for(scatvec::iterator si=scat.begin();si!=scat.end();si++) {
			CartesianCoor3D q = (si->first).first;
			int xi = int (rescorX + resX * ( q.x - left ) / (right - left) );
			int yi = int (rescorY + resY * ( q.y - bottom ) / (top - bottom) );
			if (xi==resX) xi--;
			if (xi==resY) yi--;			
			// overwrite values!
			// need new version which interpolates...
			valmap[xi][yi] = log(si->second);
			if (fsmax<(si->second)) fsmax = log(si->second);
		}

		clog << "INFO>> " << "Maximum log intensity: " << fsmax << endl;
		
        for (i=0;i<256;i++) histo[i]=0;

        for (i=0;i<resX;i++) {
//              strtok(line,)                                                                                                                          
           for (j=0;j<resY;j++) {
                   int val =  int(255*valmap[i][j] / fsmax);
                   histo[val]++;
           }
        }

        int histomax=0;
        for (i=0;i<256;i++) if (histo[i]>histo[histomax]) histomax=i;

        ofstream ofile(argv[2]);

        ofile<< "/* XPM */" << endl;
        ofile<< "static char * example_xpm[] = {" << endl;
        ofile<< "\"" <<  resX << " " << resY << " " << 256 << " " << 2 << " " << "\"," << endl;
        // generate color map;                                                                                                                                 
        for (i=0;i<256;i++) {
                ofile<< "\"";
                ofile.fill('0');
                ofile <<  hex << uppercase << setw (2) << i ;
                ofile << " c #";
//              ofile << hex << uppercase << setw (2) <<  (i>histomax?(i-histomax)*(254/(254-histomax)):0) << setw (2) << (i>(histomax/2)?(i-(histomax/2))*(254/(254-(histomax/2))):0) << setw (2) << (255-i) ;                                                                                                              
                ofile << hex << uppercase << setw (2) <<  (i>(256/2)?((i-(256/2))*2):0) << setw (2) << (i>(256/4)?((i-(256/4))*4/3):0) << setw (2) << (255-i);
                ofile << "\"," << endl;
        }

        for (i=0;i<resX;i++) {
                ofile<< "\"";
                for (j=0;j<resY;j++) {
                        int val =  int(255*valmap[i][j] / fsmax);
                        ofile.fill('0');
                        ofile << setw (2) << hex << uppercase << setw (2) << val;
                }
                if (i==(resX-1)) 
                        ofile<< "\"};" << endl;
                        else ofile<< "\"," << endl;
        }
		
        ofile.close();

        return 0;
}
				
// end of file
