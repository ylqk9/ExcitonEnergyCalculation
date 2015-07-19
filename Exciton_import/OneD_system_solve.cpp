#include <iostream>
#include <iomanip>
#include <fstream>
#include "utlt.h"
#include "OneD_system_solve.h"

using namespace std;
using Eigen::Matrix;
using Eigen::Matrix2d;
using Eigen::MatrixXcd;
using Eigen::Vector2cd;

OneD_system::OneD_system (VectorXd w1, VectorXd h1, string name1, int grids, int ks, int bands):
	system_spec (w1, h1, name1), k_pointnum(ks), gridnum(grids), bound_bandnum(0), root_precison(1e-12), 
	used_bandnum(bands),energy_sketched(false), energy_solved(false), wf_solved(false)
{
	unit_cell_width = region_width.sum();
	grid_size = unit_cell_width/gridnum;
	x_region_index.resize(gridnum);
	x.resize(gridnum);
	transfer_mat.resize(region_num, 1);
	for ( int i = 0; i != region_num; ++i ) transfer_mat(i).resize(2, 2);
	for ( int i = 0; i != gridnum; ++i )  {
		x(i) = (i+0.5)*grid_size;
		int j = 0;
		double dx = x(i) - region_width(j);
		while (dx>0) {
			++j;
			dx = dx - region_width(j);
		}
		x_region_index(i)= j ;
	}
	dk=2*M_PI/unit_cell_width/(k_pointnum/2*2);
	k_init = 0.5*((k_pointnum + 1)%2);
	transfer_mat_prod.resize(2, 2);
};

int OneD_system::Solve_Energy_sketch ()
{
	int searchsteps(20000),band_elementnum;
	double step_energy(4.0*region_height.maxCoeff()/searchsteps),search_energy;
	double half_trace(0);
	band new_band_element;
	bool band_leaving(false),band_entering(false);
	band_structure_sketch.push_back(vector<band>());
	for (int i=0; i!=searchsteps-1; ++i){
		search_energy=i*step_energy;
		half_trace=energy2half_trace(search_energy);
		if (half_trace>=-1 && half_trace<=1) {
			new_band_element.energy=search_energy;
			new_band_element.phase=acos(half_trace);
			new_band_element.k_vec=new_band_element.phase/unit_cell_width;
			band_structure_sketch[bound_bandnum].push_back(new_band_element);
			band_entering=band_leaving;
			band_leaving=true;
		} else {
			band_entering=band_leaving;
			band_leaving=false;
		}
		if (band_entering==true && band_leaving==false) {
			++bound_bandnum;
			band_structure_sketch.push_back(vector<band>());
		}
	}
	//order the energy by ascending k vector
	band swap;
	for (int i=0; i!=bound_bandnum; ++i) {
		if (i%2==1) {
			band_elementnum=band_structure_sketch[i].size();
			for (int j=0; j!=(int)ceil(static_cast<float>(band_elementnum)); ++j) {
				swap=band_structure_sketch[i][j];
				band_structure_sketch[i][j]=band_structure_sketch[i][band_elementnum-j-1];
				band_structure_sketch[i][band_elementnum-j-1]=swap;
			}
		}
	}
	if (used_bandnum>bound_bandnum) cout<<"there is less band number than required."<<endl;
	return 0;
	energy_sketched = true;
}

int OneD_system::Solve_Energy ()
{
	if(!energy_sketched) Solve_Energy_sketch();
	band_structure.resize(used_bandnum, k_pointnum);
	double energy_0_end, energy_pi_end;
	band band_structure_sketch_begin, band_structure_sketch_end, band_new_element;
	int band_size, band_curvature_indicator(-1);
	band band_begin, band_end;
	for ( int j=0; j!=used_bandnum; ++j ){
		band_size=band_structure_sketch[j].size();
		band_structure_sketch_begin=band_structure_sketch[j][0];
		band_structure_sketch_end=band_structure_sketch[j][band_size-1];
		energy_0_end=bi_section_edge(band_structure_sketch_begin.k_vec, band_structure_sketch_begin.energy, j, 0);
		energy_pi_end=bi_section_edge(band_structure_sketch_end.k_vec, band_structure_sketch_end.energy, j, M_PI/unit_cell_width);
		if(energy_0_end < energy_pi_end) band_curvature_indicator = 1;
		else band_curvature_indicator = -1;
		energy_0_end*=(1+band_curvature_indicator*root_precison);                     //improve last digit (need to be improved)
		energy_pi_end*=(1-band_curvature_indicator*root_precison);                    //improve last digit (need to be improved)
		if (k_pointnum%2) {
			band_new_element.k_vec = -M_PI/unit_cell_width;
			band_new_element.energy = energy_pi_end;
			band_new_element.phase = -M_PI;
		} else {
			band_new_element.k_vec = (-k_pointnum/2+k_init)*dk;
			band_new_element.energy = bi_section(band_new_element.k_vec,energy_0_end,energy_pi_end);
			band_new_element.phase = band_new_element.k_vec*unit_cell_width;
		}
		band_structure(j,0) = band_new_element;
		for ( int i=1; i!=k_pointnum/2+(k_pointnum+1)/2; ++i ){
			band_new_element.k_vec = (i-k_pointnum/2+k_init)*dk;
			band_new_element.energy = bi_section(band_new_element.k_vec,energy_0_end,energy_pi_end);
			band_new_element.phase = band_new_element.k_vec*unit_cell_width;
			band_structure(j,i) = band_new_element;
		}
		for ( int i=k_pointnum-1; i!=k_pointnum/2; --i ) {
			band_structure(j,i) = band_structure(j,k_pointnum-1-i);
			band_structure(j,i).k_vec = -band_structure(j,i).k_vec;
			band_structure(j,i).phase = -band_structure(j,i).phase;
		}
		if (k_pointnum%2) {
			band_new_element.energy = energy_0_end;
			band_new_element.k_vec = 0;
			band_new_element.phase = 0;
			band_structure(j,k_pointnum/2) = band_new_element;
		}
	}
	energy_solved = true;
	return 0;
}

double OneD_system::energy2half_trace(double energy)
{
	double half_trace(0);
	MatrixXd transfer_mat_mul(2, 2);
	transfer_mat_mul(0,0)=1.0;transfer_mat_mul(1,1) = 1.0;
	transfer_mat_mul(1,0)=0.0;transfer_mat_mul(0,1) = 0.0;
	PrepareTransferMat(energy);
	for (int j=0; j!=region_num; ++j) transfer_mat_mul = transfer_mat_mul*transfer_mat(j);
	half_trace = transfer_mat_mul.trace()/2;
	return half_trace;
}

double OneD_system::bi_section(double k, double upper, double lower)
{
	double error(fabs(upper-lower)/(upper+lower)),middle,middle_f,upper_f,lower_f;
	double sign(1);
	if(k<0)sign=-sign; 
	upper_f=sign*acos(energy2half_trace(upper))/unit_cell_width-k;
	lower_f=sign*acos(energy2half_trace(lower))/unit_cell_width-k;
	while (error>root_precison)
	{
		middle=(upper+lower)/2;
		middle_f=sign*acos(energy2half_trace(middle))/unit_cell_width-k;
		if (upper_f*middle_f<0)
		{
			lower_f=middle_f;
			lower=middle;
		} else {
			upper_f=middle_f;
			upper=middle;
		}
		error=fabs((upper-lower)/(upper+lower));
	}
	return (upper+lower)/2;
}

double OneD_system::bi_section_edge(double k_check, double energy_check, int band_index, double k_edge)
{
	int side((k_edge!=0)?1:-1);
	side=(band_index%2!=1)?side:-side;
	double error(fabs(k_check-k_edge)/(k_check+k_edge)),energy_edge(0),half_trace;
	double energy_shift,energy_middle;
	bool acos_nan(false);
	energy_shift=energy_check;
	while (!acos_nan)
	{
		energy_shift=energy_shift*(1+side*0.01);
		half_trace=energy2half_trace(energy_shift);
		acos_nan=(isnan(acos(half_trace)))?true:false;
	}
	while (error>root_precison/100)
	{
		energy_middle=(energy_shift+energy_check)/2;
		half_trace=energy2half_trace(energy_middle);
		acos_nan=(isnan(acos(half_trace)))?true:false;
		if (acos_nan) energy_shift=energy_middle;
		else energy_check=energy_middle;
		error=fabs((energy_check-energy_shift)/(energy_check+energy_shift));
	}
	return energy_check;
}

int OneD_system::Write_Energy_Sketch()
{
	if(!energy_sketched) Solve_Energy_sketch();
	string file_name = name + "_sketch.dat";
	ofstream file_SketchEnergy(file_name.c_str(), ofstream::out);
	file_SketchEnergy.precision(16);
	if (!file_SketchEnergy) {
		cerr<< "error: unable to open input file:"<<endl;
		return -1;
	}
	for ( int j=0; j!=band_structure_sketch.size(); ++j) {
		for ( int i=0; i!=band_structure_sketch[j].size(); ++i) {
			file_SketchEnergy<<setw(25)<<band_structure_sketch[j][i].k_vec<<" "<<setw(25)<<band_structure_sketch[j][i].energy<<endl;
		}
		file_SketchEnergy<<endl;
	}
	file_SketchEnergy.clear();
	file_SketchEnergy.close();
	return 0;
}

int OneD_system::Write_Energy()
{
	if(!energy_solved) Solve_Energy();
	string file_name = name + "_energy_band.dat";
	ofstream file_Energy(file_name.c_str(), ofstream::out);
	file_Energy.precision(16);
	if (!file_Energy) {
		cerr<< "error: unable to open input file:"<<endl;
		return -1;
	}
	for ( int j=0; j!=used_bandnum; ++j ) {
		for ( int i=0; i!=k_pointnum; ++i ) {
			file_Energy<<setw(25)<<band_structure(j,i).k_vec<<" "<<setw(25)<<band_structure(j,i).energy<<" "<<setw(25)<<band_structure(j,i).phase<<endl;
		}
		file_Energy<<endl;
	}
	file_Energy.clear();
	file_Energy.close();
	return 0;
}

void OneD_system::PrepareTransferMat(double energy)
{
	MatrixXd transfer_next(2, 2);
	double k_in_region=0;
	for (int i=0; i!=region_num; ++i){
		if (energy>region_height(i)) {
			k_in_region=sqrt(2*(energy-region_height(i)));
			transfer_next(0,0)=cos(k_in_region*region_width(i));
			transfer_next(0,1)=-sin(k_in_region*region_width(i))/k_in_region;
			transfer_next(1,0)=sin(k_in_region*region_width(i))*k_in_region;
			transfer_next(1,1)=cos(k_in_region*region_width(i));
		} else if (energy<region_height(i)) {
			k_in_region=sqrt(2*(region_height(i)-energy));
			transfer_next(0,0)=cosh(k_in_region*region_width(i));
			transfer_next(0,1)=-sinh(k_in_region*region_width(i))/k_in_region;
			transfer_next(1,0)=-sinh(k_in_region*region_width(i))*k_in_region;
			transfer_next(1,1)=cosh(k_in_region*region_width(i));
		} else {
			k_in_region=0;
			transfer_next(0,0)=1;
			transfer_next(0,1)=-region_width(i);
			transfer_next(1,0)=0;
			transfer_next(1,1)=1;
		}
		transfer_mat(i)=transfer_next;
	}
	return;
}

void OneD_system::LoopProdTransferMat(double energy)
{
	PrepareTransferMat(energy);
	transfer_mat_prod(0,0)=1;transfer_mat_prod(1,1)=1;
	transfer_mat_prod(0,1)=0;transfer_mat_prod(1,0)=0;
	for ( int i=1; i!=region_num; ++i ) {
		transfer_mat_prod=transfer_mat_prod*transfer_mat(i);
	}
	transfer_mat_prod=transfer_mat_prod*transfer_mat(0);
}

int OneD_system::Solve_wf_dens()
{
	if(!energy_solved) Solve_Energy();
	Matrix<VectorXcd, Dynamic, 1>  coeffs;
	coeffs.resize(region_num, 1);
	for ( int i=0; i != region_num; ++i ) coeffs(i).resize(2, 2);
	for ( int i=0; i != used_bandnum; ++i ) {
		for ( int j=0; j!=k_pointnum; ++j ) {
			band_structure(i,j).wf.resize(gridnum);
			band_structure(i,j).density.resize(gridnum);
			LoopProdTransferMat(band_structure(i,j).energy);
			double t11(transfer_mat_prod(0,0));
			double t22(transfer_mat_prod(1,1));
			double t21(transfer_mat_prod(1,0));
			coeffs(0)(0) = (t11-t22)/2+dcmplx(0,1)*sin(band_structure(i,j).phase);
			coeffs(0)(1) = t21;
			for ( int k=1; k!=region_num; ++k ) {
				transfer_mat(k)(0,1) = -transfer_mat(k)(0,1);
				transfer_mat(k)(1,0) = -transfer_mat(k)(1,0);				
				coeffs(k) = transfer_mat(k)*coeffs(k-1);
			}
			for ( int k=0; k!=gridnum; ++k ) {
				double k1(sqrt(2*abs(region_height(x_region_index(k)) - band_structure(i,j).energy)));
				dcmplx wf_x_cos, wf_x_sin;
				double region_sum(0);
				for ( int l=0; l!=x_region_index(k)+1; ++l ) {
					region_sum += region_width(l);
				}
				if (band_structure(i,j).energy>region_height(x_region_index(k))) {
					wf_x_cos = coeffs(x_region_index(k))(0)*cos(k1*(x(k)-region_sum));
					wf_x_sin = coeffs(x_region_index(k))(1)*sin(k1*(x(k)-region_sum))/k1;
				} else if (band_structure(i,j).energy<region_height(x_region_index(k))) {
					wf_x_cos = coeffs(x_region_index(k))(0)*cosh(k1*(x(k)-region_sum));
					wf_x_sin = coeffs(x_region_index(k))(1)*sinh(k1*(x(k)-region_sum))/k1;
				} else {
					wf_x_cos = coeffs(x_region_index(k))(0);
					wf_x_sin = coeffs(x_region_index(k))(1)*(x(k)-region_sum);
				}
				band_structure(i,j).wf(k) = wf_x_cos+wf_x_sin;
			}
			band_structure(i,j).renormalization(grid_size);
			band_structure(i,j).density = (band_structure(i,j).wf.conjugate().array()*band_structure(i,j).wf.array()).real();
			
		}
	}
	wf_solved = true;
	return 0;
}

int OneD_system::Write_wf_dens()
{
	if(!wf_solved) Solve_wf_dens();
	string file_name = name + "_wf.dat";
	ofstream file_wf(file_name.c_str(), ofstream::out);
	file_wf.precision(16);
	if (!file_wf) {
		cerr<< "error: unable to open create function file to write:"<<endl;
		return -1;
	}
	for ( int j=0; j!=used_bandnum; ++j ) {
		for ( int i=0; i!=k_pointnum; ++i ) {
			for ( int k=0; k!=gridnum; ++k ) {
				file_wf.precision(6);
				file_wf<<setw(10)<<x(k);
				file_wf.precision(16);
				file_wf<<setw(30)<<band_structure(j,i).wf(k).real();
				file_wf<<setw(30)<<band_structure(j,i).wf(k).imag();
				file_wf<<" "<<i<<" "<<j<<endl;
			}
			file_wf<<endl;
		}
	}
	file_wf.clear();
	file_wf.close();
	file_name = name + "_density.dat";
	ofstream file_density(file_name.c_str(), ofstream::out);
	file_density.precision(16);
	if (!file_density) {
		cerr<< "error: unable to create density file to write:"<<endl;
		return -1;
	}
	for ( int j=0; j!=used_bandnum; ++j ) {
		for ( int i=0; i!=k_pointnum; ++i ) {
			for ( int k=0; k!=gridnum; ++k ) {
				file_density.precision(6);
				file_density<<setw(10)<<x(k);
				file_density.precision(16);
				file_density<<setw(30)<<band_structure(j,i).density(k);
				file_density<<" "<<i<<" "<<j<<endl;
			}
			file_density<<endl;
		}
	}
	file_density.clear();
	file_density.close();
	return 0;
}