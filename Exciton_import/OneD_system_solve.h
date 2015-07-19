#include <string>
#include <vector>
#include <Eigen/Dense>
#include "utlt.h"

//header guard at start of header file
#ifndef ONED_SYSTEM_SOLVE_H
#define ONED_SYSTEM_SOLVE_H

using Eigen::Matrix;
//using Eigen::Matrix2d;
using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::Dynamic;
using std::vector;

class OneD_system : public system_spec {
public:
	// band----k_vec,energy,phase
	// band_structure.resize(used_bandnum, k_pointnum);
	Matrix<band, Dynamic, Dynamic> band_structure;
	vector<vector<band> >band_structure_sketch;
	int bound_bandnum, used_bandnum, k_pointnum, gridnum;
	int Solve_Energy_sketch ();
	int Solve_Energy ();
	int Solve_wf_dens ();
	int Write_Energy_Sketch();
	int Write_Energy ();
	int Write_wf_dens ();
	OneD_system (VectorXd w1, VectorXd h1, string name1, int grids, int ks, int bands);
private:
	double root_precison, unit_cell_width, grid_size, dk, k_init;
	bool energy_sketched, energy_solved, wf_solved;
	VectorXi x_region_index;
	VectorXd x;
	Matrix<MatrixXd, Dynamic, 1> transfer_mat;
	MatrixXd transfer_mat_prod;
	double energy2half_trace(double energy);
	double bi_section(double k, double upper, double lower);
	double bi_section_edge(double k_check, double energy_check, int band_index, double k_edge);
	void PrepareTransferMat(double energy);
	void LoopProdTransferMat(double energy);
};
//End guard at bottom of header file
#endif 