#include <vector>
#include <iomanip>
#include <fstream>
#include "utlt.h"
#include "OneD_exciton_TDDFT_Linear_resp.h"

using std::norm;
using std::ofstream;
using std::string;
using std::cerr;
using std::cout;
using std::endl;
using std::setw;
using std::conj;
using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::VectorXcd;
using Eigen::SelfAdjointEigenSolver;

//exciton base class-----------------------------------
///////////////////////////////////////////////////////////////////////////////
exciton_base::exciton_base ( const OneD_system &sys, int e_from, int e_to ) :
	name(sys.name), excit_from (e_from), excit_to (e_to), k_pointnum ( sys.k_pointnum ),
	dx (sys.region_width.sum()/sys.gridnum), unit_cell_length (sys.region_width.sum()), 
	shift_solved(false), csys(sys)
{
	gridnum = sys.gridnum;
	dk = 2.0*M_PI/unit_cell_length/(k_pointnum/2*2);
	excitation_unshift.resize ( k_pointnum );
	excitation_shift.resize ( k_pointnum );
	excitation_evectors.resize ( k_pointnum,k_pointnum );
	excitation_mix.resize ( k_pointnum,k_pointnum );
	oscillator_strength.resize ( k_pointnum );
	f_kq.resize ( k_pointnum,k_pointnum );
	for ( int i=0; i!=k_pointnum; ++i ) {
		excitation_unshift(i) = sys.band_structure(excit_to,i).energy - sys.band_structure(excit_from,i).energy;
	}
}

void exciton_base::Solve_excitation_shift ()
{
	MatrixXcd Omega_F(f_kq);
	// Omega_F.resize(f_kq.rows(), f_kq.cols());
	// Omega_F = f_kq;
	for ( int i = 0; i != k_pointnum; ++i ) Omega_F(i,i) += excitation_unshift(i);
	SelfAdjointEigenSolver<MatrixXcd> eigen_solver;
	eigen_solver.compute(Omega_F);
	excitation_shift = eigen_solver.eigenvalues();
	excitation_evectors = eigen_solver.eigenvectors();
	for ( int i=0; i!=k_pointnum; ++i ) {
		excitation_mix.col(i) = (excitation_evectors.col(i).conjugate().array()*excitation_evectors.col(i).array()).real();
	}
	shift_solved = true;
}

int exciton_base::Write_excitation_shift ()
{
	if(!shift_solved) Solve_excitation_shift();
	string file_name = name + "_exciton.dat";
	ofstream file_exciton(file_name.c_str(), ofstream::out);
	if (!file_exciton) {
		cerr<< "error: unable to create exciton file to write:"<<endl;
		return -1;
	}
	for ( int j=0; j!=k_pointnum; ++j ) {
		file_exciton.precision(16);
		file_exciton<<setw(25)<<excitation_unshift(j)<<" "<<0.5<<" "<<setw(25)<<excitation_shift(j)<<" "<<0.8<<endl;
	}
	file_exciton.clear();
	file_exciton.close();
}

int exciton_base::Write_excitation_mix ()
{
	if(!shift_solved) Solve_excitation_shift();
	string file_name = name + "_exciton_mix.dat";
	ofstream file_mix(file_name.c_str(), ofstream::out);
	if (!file_mix) {
		cerr<< "error: unable to create excitation_mix file to write:"<<endl;
		return -1;
	}
	for ( int j=0; j!=k_pointnum; ++j ) {
		for ( int i=0; i!=k_pointnum; ++i ) {
			file_mix.precision(16);
			file_mix<<setw(25)<<excitation_mix(i,j)<<endl;
		}
		file_mix<<endl;
	}
	file_mix.clear();
	file_mix.close();
}

int exciton_base::Write_oscillator_strength ()
{
	string file_name = name + "_osci.dat";
	ofstream file_osci(file_name.c_str(), ofstream::out);
	if (!file_osci) {
		cerr<< "error: unable to create oscillator_strength file to write:"<<endl;
		return -1;
	}
	for ( int j=0; j!=k_pointnum; ++j ) {
		file_osci.precision(16);
		file_osci<<setw(25)<<excitation_unshift[j]<<" "<<setw(25)<<oscillator_strength(j)<<endl;
	}
	file_osci.clear();
	file_osci.close();
}

int exciton_base::Write_matrix_fkq()
{
	string file_name = name + "_fkq.dat";
	ofstream file_matrix_fkq(file_name.c_str(), ofstream::out);
	if (!file_matrix_fkq) {
		cerr<< "error: unable to create matrix fkq file to write:"<<endl;
		return -1;
	}
	for ( int j=0; j!=k_pointnum; ++j ) {
		for ( int i=0; i!=k_pointnum; ++i ) {
			file_matrix_fkq.precision(16);
			file_matrix_fkq<<"("<<setw(25)<<f_kq(i,j).real()<<","<<setw(25)<<f_kq(i,j).imag()<<")"<<" | ";
		}
		file_matrix_fkq<<endl;
	}
	file_matrix_fkq.clear();
	file_matrix_fkq.close();
	file_name = name + "_fkq_abs.dat";
	ofstream file_matrix_fkq_abs(file_name.c_str(), ofstream::out);
	if (!file_matrix_fkq_abs) {
		cerr<< "error: unable to create matrix_fkq_abs file to write:"<<endl;
		return -1;
	}
	for ( int j=0; j!=k_pointnum; ++j ) {
		for ( int i=0; i!=k_pointnum; ++i ) {
			file_matrix_fkq_abs.precision(16);
			file_matrix_fkq_abs<<norm(f_kq(i,j))<<endl;
		}
		file_matrix_fkq_abs<<endl;
	}
	file_matrix_fkq_abs.clear();
	file_matrix_fkq_abs.close();
}

int exciton_base::Write_eff_fkq()
{
	string file_name = name + "_eff_fkq.dat";
	ofstream file_eff_fkq(file_name.c_str(), ofstream::out);
	if (!file_eff_fkq) {
		cerr<< "error: unable to create matrix fkq file to write:"<<endl;
		return -1;
	}
	for ( int j=0; j!=k_pointnum; ++j ) {
		for ( int i=0; i!=k_pointnum; ++i ) {
			file_eff_fkq.precision(16);
			file_eff_fkq<<f_kq(i,j).real() + f_kq(k_pointnum-1-i,j).imag()<<endl;
		}
		file_eff_fkq<<endl;
	}
	file_eff_fkq.clear();
	file_eff_fkq.close();
}

int exciton_base::Write_sum_fkq()
{
	// MatrixXd eff_fkq(f_kq);
	string file_name = name + "_sum_fkq_column.dat";
	ofstream file_sum_fkq(file_name.c_str(), ofstream::out);
	if (!file_sum_fkq) {
		cerr<< "error: unable to create matrix fkq file to write:"<<endl;
		return -1;
	}
	// for ( int j=0; j!=k_pointnum; ++j ) {
		// for ( int i=0; i!=k_pointnum; ++i ) {
			// eff_fkq(i,j) = f_kq(i,j).real() + f_kq(k_pointnum-1-i,j).imag()<<endl;
		// }
	// }
	for ( int i=0; i!=k_pointnum; ++i ) {
		file_sum_fkq.precision(16);
		file_sum_fkq<<f_kq.cwiseAbs().col(i).sum()<<"  "<<excitation_unshift(i)<<endl;
	}
	file_sum_fkq.clear();
	file_sum_fkq.close();
}

//exciton non-interacting-------------------------------
///////////////////////////////////////////////////////////////////////////////
exciton_non_interaction::exciton_non_interaction ( const OneD_system &sys, int e_from, int e_to ) :
		exciton_base( sys, e_from, e_to )
{
}

void exciton_non_interaction::Solve_non_interaction_f_kq ()
{
	f_kq = MatrixXcd::Zero(k_pointnum,k_pointnum);
}

void exciton_non_interaction::Solve_f_kq ()
{
	exciton_non_interaction::Solve_non_interaction_f_kq();
}

void exciton_non_interaction::Cal_oscillator_strength ()
{
	for ( int i=0; i!=k_pointnum; ++i ) {
		double x;
		dcmplx temp(0,0);
		const VectorXcd &wfs_c(csys.band_structure(excit_to,i).wf);
		const VectorXcd &wfs_v(csys.band_structure(excit_from,i).wf);
		for ( int j=0; j!=gridnum; ++j ) {
			x=(j+0.5)*dx;
			temp += x*dx*wfs_c(j)*conj(wfs_v(j));
		}
		oscillator_strength(i) = 2.0*excitation_unshift(i)*norm(temp);
	}
}

//Dummy f_kq (unknown interacting)------------------------
///////////////////////////////////////////////////////////////////////////////
exciton_dummy_f_kq::exciton_dummy_f_kq ( const OneD_system &sys, int e_from, int e_to ) :
		exciton_base( sys, e_from, e_to )
{
}

void exciton_dummy_f_kq::Solve_dummy_f_kq ( double amp )
{
	double temp,dk_sq(pow(dk,2)),sigma(150.0),sigma_sq(pow(sigma,2));
	for ( int j=0; j!=k_pointnum; ++j ) {
		for ( int i=0; i!=k_pointnum; ++i ) {
			temp = -dk*amp*exp((-pow(static_cast<double>(i-k_pointnum/2),2)*dk_sq-pow(static_cast<double>(j-k_pointnum/2),2)*dk_sq)/sigma_sq)/unit_cell_length;
			f_kq(i,j) = dcmplx(temp,0);
		}
	}
}

void exciton_dummy_f_kq::Solve_f_kq (double amp)
{
	exciton_dummy_f_kq::Solve_dummy_f_kq (amp);
}

//TDDFT with Tamm Dancoff delta----------------------------
///////////////////////////////////////////////////////////////////////////////
exciton_TDDFT_Tamm_Dancoff_delta::exciton_TDDFT_Tamm_Dancoff_delta ( const OneD_system &sys, int e_from, int e_to ) :
		exciton_base( sys, e_from, e_to )
{
}

void exciton_TDDFT_Tamm_Dancoff_delta::Solve_delta_kernel_f_kq ( double amp )
{
	dcmplx temp;
	double pre_fac(-dk*amp*2/pow(unit_cell_length,2.0));
	for ( int i=0; i!=k_pointnum; ++i ) {
		for ( int j=0; j!=k_pointnum; ++j ) {
			temp=dcmplx(0,0);
			const VectorXcd &wfs_c_i(csys.band_structure(excit_to,i).wf);
			const VectorXcd &wfs_c_j(csys.band_structure(excit_to,j).wf);
			const VectorXcd &wfs_v_i(csys.band_structure(excit_from,i).wf);
			const VectorXcd &wfs_v_j(csys.band_structure(excit_from,j).wf);
			for ( int k=0; k!=gridnum; ++k ) {
				temp += dx*conj(wfs_c_i(k))*wfs_v_i(k)*conj(wfs_v_j(k))*wfs_c_j(k);
			}
			f_kq(j,i) = pre_fac*temp;
		}
	}
}

void exciton_TDDFT_Tamm_Dancoff_delta::Solve_f_kq ( double amp )
{
	exciton_TDDFT_Tamm_Dancoff_delta::Solve_delta_kernel_f_kq( amp );
}


    // k_num=size(wfs_v,2);gridnum=size(wfs_v,1)
    // delta=(a+b)/dble(gridnum)
    // fkq=(0.D0,0.D0)
    // do i=1,k_num
      // do j=1,k_num
        // tempre=0.D0;tempim=0.D0
        // do l=1,gridnum
          // temp=dconjg(wfs_c(l,i))*wfs_v(l,i)*dconjg(wfs_v(l,j))*wfs_c(l,j)
          // tempre=tempre+dble(temp)*delta
          // tempim=tempim+dimag(temp)*delta
        // enddo
        // fkq(i,j)=dcmplx(tempre,tempim)
      // enddo
    // enddo
    // fkq=-dk*AA*2.D0/(a+b)**2*fkq


//TDDFT with Tamm Dancoff soft Coulumb----------------------------
///////////////////////////////////////////////////////////////////////////////
exciton_TDDFT_Tamm_Dancoff_Soft_Coulumb::exciton_TDDFT_Tamm_Dancoff_Soft_Coulumb ( const OneD_system &sys, int e_from, int e_to ) :
		exciton_base( sys, e_from, e_to )
{
}

void exciton_TDDFT_Tamm_Dancoff_Soft_Coulumb::Solve_Soft_Coulumb_kernel_f_kq ( double amp )
{
	dcmplx temp;
	double pre_fac(-dk*amp*2*dx/pow(unit_cell_length,2.0));
	for ( int i=0; i!=k_pointnum; ++i ) {
		for ( int j=0; j!=k_pointnum; ++j ) {
			temp=dcmplx(0,0);
			const VectorXcd &wfs_c_i(csys.band_structure(excit_to,i).wf);
			const VectorXcd &wfs_c_j(csys.band_structure(excit_to,j).wf);
			const VectorXcd &wfs_v_i(csys.band_structure(excit_from,i).wf);
			const VectorXcd &wfs_v_j(csys.band_structure(excit_from,j).wf);
			for ( int k=0; k!=gridnum; ++k ) {
				for ( int l=0; l!=gridnum; ++l ) {
					temp += conj(wfs_c_i(k))*wfs_v_i(k)*(1/sqrt(pow((l-k)*dx,2)+0.01))*conj(wfs_v_j(l))*wfs_c_j(l);
				}
			}
			f_kq(j,i) = pre_fac*temp;
		}
	}
}

void exciton_TDDFT_Tamm_Dancoff_Soft_Coulumb::Solve_f_kq ( double amp )
{
	exciton_TDDFT_Tamm_Dancoff_Soft_Coulumb::Solve_Soft_Coulumb_kernel_f_kq( amp );
}

exciton_TDDFT_Tamm_Dancoff_any_kernel::exciton_TDDFT_Tamm_Dancoff_any_kernel( const OneD_system &sys, int e_from, int e_to ) :
		exciton_base( sys, e_from, e_to )
{
}

void exciton_TDDFT_Tamm_Dancoff_any_kernel::Solve_any_kernel_f_kq ( double amp )
{
	dcmplx temp;
	double pre_fac(-dk*amp*2*dx/pow(unit_cell_length,2.0));
	for ( int i=0; i!=k_pointnum; ++i ) {
		for ( int j=0; j!=k_pointnum; ++j ) {
			temp=dcmplx(0,0);
			const VectorXcd &wfs_c_i(csys.band_structure(excit_to,i).wf);
			const VectorXcd &wfs_c_j(csys.band_structure(excit_to,j).wf);
			const VectorXcd &wfs_v_i(csys.band_structure(excit_from,i).wf);
			const VectorXcd &wfs_v_j(csys.band_structure(excit_from,j).wf);
			for ( int k=0; k!=gridnum; ++k ) {
				for ( int l=0; l!=gridnum; ++l ) {
					temp += conj(wfs_c_i(k))*wfs_v_i(k)*(Kernel((k-0.5)*dx,(l-0.5)*dx))*conj(wfs_v_j(l))*wfs_c_j(l);
				}
			}
			f_kq(j,i) = pre_fac*temp;
		}
	}
	
}

double exciton_TDDFT_Tamm_Dancoff_any_kernel::Kernel ( double r, double rp )
{
	// double val(1/sqrt(pow((r-rp),2)+1));
	// double val(pow(r/2.0,2)+pow(rp/2.0,2));
	// double val(abs(r-rp));
	// double val(abs(3.0-r-rp));
	// double val(1.0);
	// double val(1/sqrt(pow((r+rp),2)+1));
	// double val(-1/sqrt(pow((2.0-r-rp),2)+1));
	// double val(1/sqrt(pow((r-rp),2)+100));
	// double val(1/sqrt(pow((r-rp),2)+0.01));
	// double val(1/sqrt(pow((r-rp),2)+0.001));
	// double val(1/sqrt(pow((r-rp),2)+1)+10.0);
	// double val(1/sqrt(pow((r-rp),2)+0.1)+20.0);
	// double val(1/sqrt(pow((r-rp+0.3),2)+0.1));
	// double val(1/sqrt(pow((r-rp),4)+0.1));
	// double val(1/sqrt(pow((r-rp),4)+1));
	// double val(1/sqrt(pow((r-rp),8)+0.01));
	// double val(1/sqrt(pow((r-rp+0.7),2)+0.1)+1/sqrt(pow((rp-r+0.7),2)+0.1));
	// double val(1/sqrt(pow((r-rp+1.5),4)+0.01)+1/sqrt(pow((rp-r+1.5),4)+0.01));
	// double val(1/sqrt(pow((r-rp+0.4),2)+0.1)+1/sqrt(pow((rp-r+0.4),2)+0.1));
	// double val(1/sqrt(pow((r-rp+0.4),2)+0.01)+1/sqrt(pow((rp-r+0.4),2)+0.01));
	// double val(1/sqrt(pow((r-rp+2.4),2)+0.1)+1/sqrt(pow((rp-r+2.4),2)+0.1));
	// double val(1/sqrt(pow((r-rp+3.4),3)+1)+1/sqrt(pow((rp-r+3.4),3)+1));
	// double val(1/sqrt(pow((r-rp+3.4),3)+20)+1/sqrt(pow((rp-r+3.4),3)+20));
	// double val(-10.0/sqrt(pow((r-rp+3.4),3)+20)-10.0/sqrt(pow((rp-r+3.4),3)+20));
	// double val(-10.0/sqrt(pow((r-rp+3.4),3)+20)-10.0/sqrt(pow((rp-r+3.4),3)+20));
	// double val(1.0/sqrt(pow((3-r-rp+3.4),3)+20)+1.0/sqrt(pow((-3+rp+r+3.4),3)+20));
	// double val(10.0/sqrt(pow((3-r-rp+3.4),3)+20)+10.0/sqrt(pow((-3+rp+r+3.4),3)+20));
	// double val(10.0/sqrt(pow((3-r-rp+0.4),3)+20)+10.0/sqrt(pow((-3+rp+r+0.4),3)+20));
	double val(50.0*pow(r,2)+50.0*pow(rp,2));
	return val;
}

void exciton_TDDFT_Tamm_Dancoff_any_kernel::Solve_f_kq ( double amp )
{
	exciton_TDDFT_Tamm_Dancoff_any_kernel::Solve_any_kernel_f_kq( amp );
}	
//
///////////////////////////////////////////////////////////////////////////////