#include <complex>
#include <string>
#include <Eigen/Dense>

//header guard at start of header file
#ifndef UTLT_H
#define UTLT_H

#define isnan(x) _isnan(x)

using Eigen::VectorXd;
using Eigen::VectorXcd;
using std::string;

//complex vector and complex 2D matrix
typedef std::complex<double> dcmplx;

class band {
public:
	int band_index;
	double energy,phase,k_vec;
	VectorXcd wf;
	VectorXd density;
	void renormalization(double dx);
};

class system_spec {
public:
	int region_num;
	string name;
	VectorXd region_width;
	VectorXd region_height;
	int Write_system_spec ();
	system_spec (VectorXd dwidth, VectorXd dheight, string name1);
	system_spec (int num, string express);
};

//End guard at bottom of header file
#endif 