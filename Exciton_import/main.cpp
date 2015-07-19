#include <iostream>
#include <string>
#include <windows.h>
#include <time.h>
#include <Eigen/Dense>
#include "utlt.h"
#include "OneD_system_solve.h"
#include "OneD_exciton_TDDFT_Linear_resp.h"

using namespace std;
using Eigen::VectorXd;

int main()
{
	clock_t t_begin,t_end;
	t_begin=clock();
	int band_considered(3), gridnum(150), kpointnum(151);
	double amplitude1(3.5645),amplitude2(1.0);
	VectorXd w1(2), h1(w1.rows());
	// w1<<2.4, 0.6; h1<<0.0, 8.0;
	w1<<0.5, 0.5; h1<<0.0, 20.0;
	VectorXd w2(4), h2(w2.rows());
	w2<<2.4, 0.6, 2.8, 0.2;
	h2<<0.0, 8.0, 0.0, 10.0;
	// --------------------------------------------------------------------------------
	OneD_system KP_model(w1, h1, "My_kp_model", gridnum, kpointnum, band_considered);
	cout<<KP_model.name<<" is running..."<<endl;
	KP_model.Solve_wf_dens();
	KP_model.Write_Energy();
	KP_model.Write_system_spec ();
	exciton_TDDFT_Tamm_Dancoff_Soft_Coulumb KP_exciton( KP_model,1,2 );
	KP_exciton.Solve_f_kq (0.01);
	KP_exciton.Write_excitation_shift ();
	KP_exciton.Write_excitation_mix ();
	KP_exciton.Write_matrix_fkq ();
	KP_exciton.Write_eff_fkq ();
	KP_exciton.Write_sum_fkq();
	//--------------------------------------------------------------------------------
	// OneD_system Double_model(w2, h2, "dmodel", gridnum, kpointnum, band_considered);
	// cout<<Double_model.name<<" is running..."<<endl;
	// Double_model.Write_system_spec ();
	// Double_model.Solve_wf_dens ();
	// exciton_TDDFT_Tamm_Dancoff_delta Dw_exciton( Double_model,1,2 );
	// Dw_exciton.Solve_f_kq (amplitude1);
	// Dw_exciton.Write_excitation_shift ();
	// Dw_exciton.Write_excitation_mix ();
	// Dw_exciton.Write_matrix_fkq ();
	// Dw_exciton.Write_eff_fkq ();
	// Dw_exciton.Write_sum_fkq();
	//--------------------------------------------------------------------------------
	// OneD_system scmodel(w1, h1, "scmodel", gridnum, kpointnum, band_considered);
	// cout<<scmodel.name<<" is running..."<<endl;
	// scmodel.Write_system_spec ();
	// scmodel.Solve_wf_dens ();
	// exciton_TDDFT_Tamm_Dancoff_Soft_Coulumb sc_exciton( scmodel,1,2 );
	// sc_exciton.Solve_f_kq (amplitude2);
	// sc_exciton.Write_excitation_shift ();
	// sc_exciton.Write_excitation_mix ();
	// sc_exciton.Write_matrix_fkq ();
	// sc_exciton.Write_eff_fkq ();
	// sc_exciton.Write_sum_fkq();
	//---------------------------------------------------------------------------------
	// OneD_system anymodel(w1, h1, "anymodel", gridnum, kpointnum, band_considered);
	// cout<<anymodel.name<<" is running..."<<endl;
	// anymodel.Write_system_spec ();
	// anymodel.Solve_wf_dens ();
	// exciton_TDDFT_Tamm_Dancoff_any_kernel any_exciton( anymodel,1,2 );
	// any_exciton.Solve_f_kq (amplitude1);
	// any_exciton.Write_excitation_shift ();
	// any_exciton.Write_excitation_mix ();
	// any_exciton.Write_matrix_fkq ();
	// any_exciton.Write_eff_fkq ();
	// any_exciton.Write_sum_fkq();
	//---------------------------------------------------------------------------------
	Eigen::Matrix2cd a;
	a<<dcmplx(1,2),dcmplx(2,3),dcmplx(1,4),dcmplx(2,4);
	cout<<a<<endl;
	cout<<a.cwiseAbs()<<endl;
	t_end=clock();
	cout<<"program running time "<<(double(t_end)-double(t_begin))/CLOCKS_PER_SEC<<endl;
	MessageBox(NULL, L"Work finished!", L"Notice", MB_OK);
	return 0;
}