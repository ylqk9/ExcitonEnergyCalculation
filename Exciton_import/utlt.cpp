#include <string>
#include <iostream>
#include <fstream>
#include "utlt.h"

using std::string;
using std::endl;
using std::cerr;
using std::ofstream;

void band::renormalization(double dx)
{
	double norm1 = (wf.adjoint()*wf)(0,0).real()*dx;
	double sqrt_norm = sqrt(norm1);
	density = density/norm1;
	wf = wf/sqrt_norm;
	return;
}

system_spec::system_spec (VectorXd dwidth, VectorXd dheight, string name1):
	region_num(dwidth.rows()), name(name1) 
{
	region_width.resize(region_num); region_width = dwidth;
	region_height.resize(region_num); region_height = dheight;
};

system_spec::system_spec (int num, string express) {}
	
int system_spec::Write_system_spec ()
{
	double dx = region_width.sum()/100;
	string file_name = name + "_system.dat";
	ofstream file_system_spec(file_name.c_str(), ofstream::out);
	file_system_spec.precision(16);
	if (!file_system_spec) {
		cerr<< "error: unable to open input file:"<<endl;
		return -1;
	}
	for ( int j=0; j!=100; ++j) {
		double x = (j+0.5)*dx;
		int i=0;
		double diff=x-region_width(i);
		while (diff>0) {
			++i;
			diff=diff-region_width(i);
		}
		file_system_spec<<x<<" "<<region_height(i)<<endl;
	}
	file_system_spec.clear();
	file_system_spec.close();
	return 0;
}