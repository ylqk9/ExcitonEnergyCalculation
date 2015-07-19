#include <iostream>
#include <vector>
#include <string>
#include "utlt.h"
#include "OneD_system_Solve.h"

//header guard at start of header file
#ifndef ONED_EXCITON_TDDFT_LINEAR_RESP_H
#define ONED_EXCITON_TDDFT_LINEAR_RESP_H

using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::VectorXd;
using Eigen::VectorXcd;
using std::string;

class exciton_base {
public:
	VectorXd excitation_unshift, excitation_shift;
	VectorXd oscillator_strength;
	MatrixXd excitation_mix;
	MatrixXcd excitation_evectors;
	int excit_from,excit_to;
	MatrixXcd f_kq;
	exciton_base ( const OneD_system &sys, int e_from, int e_to );
	virtual void Solve_f_kq () {}
	void Solve_excitation_shift ();
	int Write_excitation_shift ();
	int Write_excitation_mix ();
	int Write_oscillator_strength ();
	int Write_matrix_fkq();
	int Write_eff_fkq();
	int Write_sum_fkq();
	string name;
protected:
	int gridnum, k_pointnum;
	double dk,dx,unit_cell_length,amplitude;
	const OneD_system &csys;
private:
	bool shift_solved;
};

class exciton_non_interaction : public exciton_base {
public:
	exciton_non_interaction ( const OneD_system &sys, int e_from, int e_to );
	void Solve_f_kq();
	void Solve_non_interaction_f_kq ();
	void Cal_oscillator_strength ();
};

class exciton_dummy_f_kq : public exciton_base {
public:
	exciton_dummy_f_kq ( const OneD_system &sys, int e_from, int e_to );
	void Solve_f_kq ( double amp );
	void Solve_dummy_f_kq ( double amp );
};

class exciton_TDDFT_Tamm_Dancoff_delta : public exciton_base {
public:
	exciton_TDDFT_Tamm_Dancoff_delta ( const OneD_system &sys, int e_from, int e_to );
	void Solve_f_kq ( double amp );
	void Solve_delta_kernel_f_kq ( double amp );
};

class exciton_TDDFT_Tamm_Dancoff_Soft_Coulumb : public exciton_base {
public:
	exciton_TDDFT_Tamm_Dancoff_Soft_Coulumb ( const OneD_system &sys, int e_from, int e_to );
	void Solve_f_kq (double amp );
	void Solve_Soft_Coulumb_kernel_f_kq ( double amp );
};

class exciton_TDDFT_full : public exciton_base {
};

class exciton_TDDFT_Tamm_Dancoff_any_kernel : public exciton_base {
public:
	exciton_TDDFT_Tamm_Dancoff_any_kernel ( const OneD_system &sys, int e_from, int e_to );
	void Solve_f_kq ( double amp );
	void Solve_any_kernel_f_kq ( double amp );
	double Kernel ( double r, double rp );
};
//End guard at bottom of header file
#endif 