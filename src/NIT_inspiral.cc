// NIT_inspiral - code to rapidly compute extreme mass-ratio inspirals using self-force results
// Copyright (C) 2017  Niels Warburton
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include <NIT_inspiral.h>
#include <libconfig.h++>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_min.h>
#include <math.h>
#include <iomanip>

// ComplexAmpTest include statements:
#include <fstream>
#include <iostream>
#include <complex>
 // #include <iomanip>
#include <algorithm>
 // #include "Interpolant.h"


#define Sin(x)          (sin((double)(x)))
#define Cos(x)          (cos((double)(x)))
#define Power(x, y)     (pow((double)(x), (double)(y)))
#define Sqrt(x)         (sqrt((double)(x)))
#define Pi				M_PI
#define Conjugate(x)	(conj(x))

#define CHI_MAX 1000000.0

// precision of output
#define OUT_PREC 12

// numerical tolerance for ODE integrations
#define NUM_TOL 1.0e-11

// smallest value for p-2*e-6
#define Y_MIN 0.1

using namespace libconfig;

void compute_waveform(string insp_filename, string out_filename);
Complex waveform_h(double p, double e, double v, double phi, double Z, double Phi);
Complex waveform_T(double Ar,double Ai, double v, double phi,double w_phi, double w_r,double t, double tt, double m, double n);

int mode;				// Selects the behavior of the code. Set to either FUll_INSPIRAL, NIT_INSPIRAL, DECOMPOSE or CONSTRUCT_Fs
double q;				// The mass ratio q = m_1/m_2
ofstream fout;			// The output file
Config cfg;

int main(int argc, char* argv[]){	
	double p0, e0;			// The initial (p,e) values for the evolutions
	
	if(argc  > 1){	
		if( !strcmp(argv[1], "-f") ){
			mode = FULL_INSPIRAL;
			cout << "# Mode: inspiral using full EoM" << endl;
		} else if( !strcmp(argv[1], "-f0") ){
			mode = FULL_INSPIRAL_DEFAULT;
			cout << "# Mode: inspiral using default full EoM" << endl;
		} else if( !strcmp(argv[1], "-n") ){
			mode = NIT_INSPIRAL;
			cout << "# Mode: inspiral using NIT EoM" << endl;
		} else if( !strcmp(argv[1], "-n0") ){
			mode = NIT_INSPIRAL_DEFAULT;
			cout << "# Mode: inspiral using default NIT EoM" << endl;
		} else if( !strcmp(argv[1], "-d") ){
			mode = DECOMPOSE;
			cout << "# Mode: decompose the self-force data to compute the F's and f's in NIT EoM" << endl;
		}else if( !strcmp(argv[1], "-c") ){
			mode = CONSTRUCT_Fs;
			cout << "# Mode: Compute the F's and f's in the NIT EoM" << endl;
		}
		
		
		
		else if( !strcmp(argv[1], "-w") ){
			if(argc == 6){
				p0 = atof(argv[2]);
				e0 = atof(argv[3]);	
				q  = atof(argv[4]);
				if( !strcmp(argv[5], "-n") ){
					mode = WAVEFORM_NIT; 
					cout << "Compute the NIT kludge waveform" << endl;
				}else if( !strcmp(argv[5], "-f") ){
					 mode = WAVEFORM_FULL;
					 cout << "Compute the Full kludge waveform" << endl;
				}else{
					cout << "Unrecognized waveform flag" << endl;
					exit(0);
				}
			}else{
				cout << "For the waveform mode please enter initial p, e and q values and a '-f' or '-n' flag for Full or NIT." << endl;
				exit(0);
			}
		}
		// Teukolsky waveform
				else if( !strcmp(argv[1], "-t") ){
			if(argc == 6){
				p0 = atof(argv[2]);
				e0 = atof(argv[3]);	
				q  = atof(argv[4]);
				if( !strcmp(argv[5], "-n") ){
					mode = T_WAVEFORM_NIT; 
					cout << "Compute the NIT Teukolsky waveform" << endl;
				}else if( !strcmp(argv[5], "-f") ){
					 mode = T_WAVEFORM_FULL;
					 cout << "Compute the Full Teukolsky waveform" << endl;
				}else{
					cout << "Unrecognized waveform flag" << endl;
					exit(0);
				}
			}else{
				cout << "For the waveform mode please enter initial p, e and q values and a '-f' or '-n' flag for Full or NIT." << endl;
				exit(0);
			}
		}
		
		else{
			cout << "Unrecognized flag. Run with no arguments for instructions." << endl;
		}
		if(mode == FULL_INSPIRAL || mode == FULL_INSPIRAL_DEFAULT || mode == NIT_INSPIRAL || mode == NIT_INSPIRAL_DEFAULT){
			if(argc == 5){
				p0 = atof(argv[2]);
				e0 = atof(argv[3]);	
				q  = atof(argv[4]);
			}else{
				cout << "For inspiral modes please enter initial p, e value and a q value." << endl;
				exit(0);
			}
		}
	}else{
		cout << "Necessary parameters:" << endl;
		cout << "\t1. flag      '-f', '-f0', '-n', '-n0', '-w', '-t', '-d' or '-c' " << endl;
		cout << "\t   '-f'      Full inspiral" << endl;
		cout << "\t   '-f0'     Default Full inspiral" << endl;
		cout << "\t   '-n'      NIT inspiral" << endl;
		cout << "\t   '-n0'     Default NIT inspiral" << endl;
		cout << "\t   '-w'      Compute the kludge waveform" << endl;
		cout << "\t   '-t'      Compute the Teukolsky waveform" << endl;
		cout << "\t   '-d'      Decompose the self-force data into Fourier modes" << endl;
		cout << "\t   '-c'      Compute F's and f's in NIT EoM" << endl;
		cout << "\t2. p         Initial semi-latus rectum for '-f' or '-n' inspiral options" << endl;
		cout << "\t3. e         Initial eccentricity for '-f' or '-n' inspiral options" << endl;
		cout << "\t4. q         Mass ratio for '-f' or '-n' inspiral options" << endl;
		cout << "\t5. flag      Waveform only, '-f' for full, '-n' for NIT" << endl;
		exit(0);
	}	
	
	
	// Read in the parameters from the configuration file
	cout << "Reading configuration file" << endl;
	
	// Check the config file exists and is well formatted
    try{
		 cfg.readFile("config/parameters.cfg");
	} catch(const FileIOException &fioex){
	     std::cerr << "Config file could not be found at 'config/parameters.cfg'" << std::endl;
		 exit(0);
    } catch(const ParseException &pex){
    	std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine() << " - " << pex.getError() << std::endl;
  	    exit(0);
    }
	
	
	// Load the model and check that it loaded correctly	
	int model_status = lib_Sch_GSF_load_model();
	if(model_status == MODEL_LOAD_FAIL) exit(0);
	
	// Disable the GSL error handler (only for production version of code)
	// gsl_set_error_handler_off();
	
	// Output all the data at specified digits of precision using scientific notation
	fout.precision(OUT_PREC);
	cout.precision(OUT_PREC);
	fout << scientific;
	cout << scientific;
	
	// Perform the calculation depending upon which mode is selected
	if(mode == FULL_INSPIRAL){
		
		ostringstream filename;
		filename << "output/Inspiral_Full_p" << p0 << "_e" << e0 << "_q" << q <<".dat";
		cout << "Outputting inspiral trajectory to " << filename.str() << endl;	
		
		fout.open(filename.str());
		fout << "# Full Inspiral" << endl;
		integrate_osc_eqs(p0, e0);	
			
	}else if(mode == FULL_INSPIRAL_DEFAULT){
		
		ostringstream filename;
		filename << "output/Inspiral_Full_p" << p0 << "_e" << e0 << "_q" << q <<".dat";
		cout << "Outputting inspiral trajectory to " << filename.str() << endl;	
		
		fout.open(filename.str());
		fout << "# Full Inspiral" << endl;
		integrate_osc_eqs_default(p0, e0);
			
	}else if(mode == NIT_INSPIRAL_DEFAULT){
		
		ostringstream filename;
		filename << "output/Inspiral_NIT_p" << p0 << "_e" << e0 << "_q" << q <<".dat";
		cout << "Outputting inspiral trajectory to " << filename.str() << endl;
	
		fout.open(filename.str());
		fout << "# Default NIT Inspiral" << endl;
		interpolate_Fs_and_integrate_NIT_EoM(mode, p0, e0);
			
	}else if(mode == NIT_INSPIRAL){
		
		ostringstream filename;
		filename << "output/Inspiral_NIT_p" << p0 << "_e" << e0 << "_q" << q <<".dat";
		cout << "Outputting inspiral trajectory to " << filename.str() << endl;
	
		fout.open(filename.str());
		fout << "# NIT Inspiral" << endl;
		interpolate_Fs_and_integrate_NIT_EoM(mode, p0, e0);
		
	}else if(mode == WAVEFORM_NIT){
		
		ostringstream out_filename, insp_filename;
		out_filename << "output/Waveform_NIT_p" << p0 << "_e" << e0 << "_q" << q <<".dat";
		insp_filename << "output/Inspiral_NIT_p" << p0 << "_e" << e0 << "_q" << q <<".dat";
		compute_waveform(insp_filename.str(), out_filename.str());
			
	}else if(mode == WAVEFORM_FULL){
		
		ostringstream out_filename, insp_filename;
		out_filename << "output/Waveform_Full_p" << p0 << "_e" << e0 << "_q" << q <<".dat";
		insp_filename << "output/Inspiral_Full_p" << p0 << "_e" << e0 << "_q" << q <<".dat";
		compute_waveform(insp_filename.str(), out_filename.str());
			
	}else if(mode == T_WAVEFORM_NIT){
		
		ostringstream out_filename, insp_filename;
		out_filename << "output/T_Waveform_NIT_p" << p0 << "_e" << e0 << "_q" << q <<".dat";
		insp_filename << "output/Inspiral_NIT_p" << p0 << "_e" << e0 << "_q" << q <<".dat";
		compute_waveform(insp_filename.str(), out_filename.str()); //should include Teukolsky waveform function
			
	}else if(mode == T_WAVEFORM_FULL){
		
		ostringstream out_filename, insp_filename;
		out_filename << "output/T_Waveform_Full_p" << p0 << "_e" << e0 << "_q" << q <<".dat";
		insp_filename << "output/Inspiral_Full_p" << p0 << "_e" << e0 << "_q" << q <<".dat";
		compute_waveform(insp_filename.str(), out_filename.str()); //should include Teukolsky waveform function
			
	}else if(mode == DECOMPOSE){
		FFT_self_force_over_parameter_space();
	}else if(mode == CONSTRUCT_Fs){
		construct_tilde_Fs();
	}
	
	
	// Close the output file
	fout.close();
}

// Functions to integrate the NIT EoM

// Used to pass the interpolants to the NIT_EoM function
struct interp_params{
	Interpolant *Fp1;
	Interpolant *Fe1;
	Interpolant *fv1;
	Interpolant *Fp2;
	Interpolant *Fe2;
	Interpolant *U1;
	Interpolant *V1;
};

int NIT_EoM (double v, const double y[], double f[], void *params){
	(void)(v); /* avoid unused parameter warning */
	
	struct interp_params *interps = (struct interp_params *)params;

	// p = y[0],  e = y[1]

	f[0] = q*interps->Fp1->eval(y[0]-2*y[1], y[1]) + q*q*interps->Fp2->eval(y[0]-2*y[1], y[1]);
	f[1] = q*interps->Fe1->eval(y[0]-2*y[1], y[1]) + q*q*interps->Fe2->eval(y[0]-2*y[1], y[1]);
	f[2] = 1.0 + q*interps->fv1->eval(y[0]-2*y[1], y[1]);
	f[3] = T_r(y[0], y[1])/(2.0*M_PI) + q*interps->U1->eval(y[0]-2*y[1], y[1]);
	f[4] = Phi(y[0], y[1])/(2.0*M_PI) + q*interps->V1->eval(y[0]-2*y[1], y[1]);
	
	if(y[0]-6-2*y[1] > Y_MIN) return GSL_SUCCESS;
	else return GSL_SUCCESS + 1;
}

void interpolate_Fs_and_integrate_NIT_EoM(int mode, double p0, double e0){
	
	char Ftildes_file[32];
	
	if(mode == NIT_INSPIRAL_DEFAULT) sprintf(Ftildes_file,"data/Ftildes_0.dat");
	else sprintf(Ftildes_file,"data/Ftildes.dat");
	
	ifstream Ftilde_file(Ftildes_file);
	
	// Load the data for the F/f/U/V on the RHS of the NIT EoM
	string Ftilde_string;
	vector<double> ys, es, Fp1s, Fe1s, fv1s, Fp2s, Fe2s, U1s, V1s, Xv1s, Yp1s, Ye1s;
	double y, e, Fp1, Fe1, fv1, Fp2, Fe2, U1, V1, Xv1, Yp1, Ye1;
	while(getline(Ftilde_file, Ftilde_string)){
		
		stringstream Ftilde_ss(Ftilde_string);
		
		Ftilde_ss >> y >> e >> Fp1 >> Fe1 >> fv1 >> Fp2 >> Fe2 >> U1 >> V1 >> Xv1 >> Yp1 >> Ye1;
		
		ys.push_back(y);
		es.push_back(e);
		Fp1s.push_back(Fp1);
		Fe1s.push_back(Fe1);
		fv1s.push_back(fv1);
		Fp2s.push_back(Fp2);
		Fe2s.push_back(Fe2);
		U1s.push_back(U1);
		V1s.push_back(V1);
		Xv1s.push_back(Xv1);
		Yp1s.push_back(Yp1);
		Ye1s.push_back(Ye1);
	}
	
	// Interpolate the data
	Interpolant Fp1_interp(ys, es, Fp1s);
	Interpolant Fe1_interp(ys, es, Fe1s);
	Interpolant fv1_interp(ys, es, fv1s);
	Interpolant Fp2_interp(ys, es, Fp2s);
	Interpolant Fe2_interp(ys, es, Fe2s);
	Interpolant U1_interp(ys, es, U1s);
	Interpolant V1_interp(ys, es, V1s);
	Interpolant Xv1_interp(ys, es, Xv1s);
	Interpolant Yp1_interp(ys, es, Yp1s);
	Interpolant Ye1_interp(ys, es, Ye1s);
	
	// Numerically integrate the NIT inspiral below
	
	struct interp_params interps = {&Fp1_interp, &Fe1_interp, &fv1_interp, &Fp2_interp, &Fe2_interp, &U1_interp, &V1_interp};
		
	double t0 = 0;
	double phi0 = 0;
	
	double chi = 0;

	// Initial conditions
	double ptilde0 = p0 + q*Yp1_interp.eval(p0-2*e0, e0);
	double etilde0 = e0 + q*Ye1_interp.eval(p0-2*e0, e0);
	double vtilde0 = 0  + q*Xv1_interp.eval(p0-2*e0, e0);
	double y1[5] = {ptilde0, etilde0, vtilde0, t0, phi0};
	
	// Output the initial parameters
	fout << "# Format: chi    p    e    xi    t    phi" << endl;
	fout << chi << " " << y1[0] << " " << y1[1] << " " << y1[2] << " " << y1[3] << " " << y1[4] << " " << endl;
	
	int dense_output;
	try{
		dense_output 		= cfg.lookup("Dense_output");		
 	}catch(const SettingNotFoundException &nfex){
   		cerr << "'Dense_output' setting missing from configuration file." << endl; exit(0);
 	}
	
	if(dense_output == 1){				// Dense output
		
		double n_per_orbit;
		try{
			n_per_orbit 		= cfg.lookup("n_per_orbit");		
	 	}catch(const SettingNotFoundException &nfex){
	   		cerr << "'n_per_orbit' setting missing from configuration file." << endl; exit(0);
	 	}

		double chi_i = 0;
		double delta_chi = 2.0*M_PI/n_per_orbit;
	
	    gsl_odeiv2_system sys 	= {NIT_EoM, NULL, 5, &interps};	
	    gsl_odeiv2_driver *d 	= gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, NUM_TOL, NUM_TOL, 0.0);
	
	    high_resolution_clock::time_point t1 = high_resolution_clock::now();
	
		while(1){
			chi_i += delta_chi;
			int status = gsl_odeiv2_driver_apply (d, &chi, chi_i, y1);

			if (status != GSL_SUCCESS){
				//printf ("error, return value=%d\n", status);
				break;
			}
			fout << chi << " " << y1[0] << " " << y1[1] << " " << y1[2] << " " << y1[3] << " " << y1[4] << " " << endl;
		}
		
	  high_resolution_clock::time_point t2 = high_resolution_clock::now();
	  
	  duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	  fout << "# Computing the inspiral took: " << time_span.count() << " seconds." << endl;

		gsl_odeiv2_driver_free (d);
		
	}else{		// Sparse output
	    const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk8pd;

	    gsl_odeiv2_step *s 		= gsl_odeiv2_step_alloc (T, 5);
	    gsl_odeiv2_control *c 	= gsl_odeiv2_control_y_new (NUM_TOL, NUM_TOL);
	    gsl_odeiv2_evolve *e 	= gsl_odeiv2_evolve_alloc (5);

	    gsl_odeiv2_system sys = {NIT_EoM, NULL, 5, &interps};

	    double chi_max = CHI_MAX;	//FIXME make this something like 100/q
	    double h = 1e-1;

	    high_resolution_clock::time_point t1 = high_resolution_clock::now();
		
		double chi_prev = 1.;
	    while (chi < chi_max) {
	        int status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &chi, chi_max, &h, y1);
			if (status != GSL_SUCCESS){
				//printf ("error, return value=%d\n", status);
				break;
			}

			fout << chi << " " << y1[0] << " " << y1[1] << " " << y1[2] << " " << y1[3] << " " << y1[4] << " " << endl;
			
			// Stop output of lots of data near the separatrix when the time step gets very small
			if(fabs(chi_prev/chi - 1.0) < NUM_TOL) break;
			
		    chi_prev = chi;	
	      }
		  high_resolution_clock::time_point t2 = high_resolution_clock::now();
		  
		  duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
		  fout << "# Computing the inspiral took: " << time_span.count() << " seconds." << endl;
		  
		  
		  
		  gsl_odeiv2_evolve_free (e);
		  gsl_odeiv2_control_free (c);
		  gsl_odeiv2_step_free (s);
	}
	
	
}


// Functions for computing the RHS of the NIT EoM
void FFT_self_force_over_parameter_space(){
	double p, e;
	int N = 50;
    fftw_complex *in, *out;
    fftw_plan plan;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	
	// Open file to output the Fourier data
    ofstream fv_file, Fp_file, Fe_file;
	ofstream dFp_dp_file, dFp_de_file, dFe_dp_file, dFe_de_file;
	ofstream dV0_dp_file, dV0_de_file, dV0_dv_file;
	ofstream dU0_dp_file, dU0_de_file, dU0_dv_file;
	
	fv_file.open("data/Fourier_fv.dat");
	Fp_file.open("data/Fourier_Fp.dat");
	Fe_file.open("data/Fourier_Fe.dat");
	
	dFp_dp_file.open("data/Fourier_dFp_dp.dat");	
	dFp_de_file.open("data/Fourier_dFp_de.dat");
	dFe_dp_file.open("data/Fourier_dFe_dp.dat");
	dFe_de_file.open("data/Fourier_dFe_de.dat");
	
	dV0_dp_file.open("data/Fourier_dV0_dp.dat");
	dV0_de_file.open("data/Fourier_dV0_de.dat");
	dV0_dv_file.open("data/Fourier_dV0_dv.dat");
	
	dU0_dp_file.open("data/Fourier_dU0_dp.dat");
	dU0_de_file.open("data/Fourier_dU0_de.dat");
	dU0_dv_file.open("data/Fourier_dU0_dv.dat");

	double de = 0.002;
	double dp = 0.05;
	e = 0;
	for(int i=1; i <= 100; i++){
		e += de;
		for(int j=0; j < 100; j++){
			p = 6 + 2*e + 0.05  + j*dp;
			FFT_self_force(p, e, &plan, N, in, out, &fv_file, &Fp_file, &Fe_file, &dFp_dp_file, &dFp_de_file, &dFe_dp_file, &dFe_de_file, &dV0_dp_file, &dV0_de_file, &dV0_dv_file, &dU0_dp_file, &dU0_de_file, &dU0_dv_file);
		}
	}
	
	// Close all the files after writing the data out    
	fv_file.close();
	Fp_file.close();
	Fe_file.close();	
	
	dFp_dp_file.close();	
	dFp_de_file.close();
	dFe_dp_file.close();
	dFe_de_file.close();
	
	dV0_dp_file.close();
	dV0_de_file.close();
	dV0_dv_file.close();
	
	dU0_dp_file.close();
	dU0_de_file.close();
	dU0_dv_file.close();
		
    fftw_destroy_plan(plan);
    fftw_free(in); fftw_free(out);
	
}

void FFT_self_force(double p, double e, fftw_plan *plan, int N, fftw_complex *in, fftw_complex *out, ofstream *fv_file, ofstream *Fp_file, ofstream *Fe_file, ofstream *dFp_dp_file, ofstream *dFp_de_file, ofstream *dFe_dp_file, ofstream *dFe_de_file,ofstream *dV0_dp_file, ofstream *dV0_de_file, ofstream *dV0_dv_file, ofstream *dU0_dp_file, ofstream *dU0_de_file, ofstream *dU0_dv_file){
		
	double *Fr = new double[N];
	double *Fphi = new double[N];
	
	double *dFr_dp = new double[N];
	double *dFr_de = new double[N];
	double *dFphi_dp = new double[N];
	double *dFphi_de = new double[N];
	
	double *v = new double[N];
	
	// The number of Fourier modes to output
	int N_out = 10;
	
	// Output the coordinates in phase space, y = p-2e, e 
	*Fp_file << p-2*e << " " << e << " ";
	*Fe_file << p-2*e << " " << e << " ";
	*fv_file << p-2*e << " " << e << " ";
	
	*dFp_dp_file << p-2*e << " " << e << " ";
	*dFp_de_file << p-2*e << " " << e << " ";
	*dFe_dp_file << p-2*e << " " << e << " ";
	*dFe_de_file << p-2*e << " " << e << " ";
	
	*dV0_dp_file << p-2*e << " " << e << " ";
	*dV0_de_file << p-2*e << " " << e << " ";
	*dV0_dv_file << p-2*e << " " << e << " ";
	
	*dU0_dp_file << p-2*e << " " << e << " ";
	*dU0_de_file << p-2*e << " " << e << " ";
	*dU0_dv_file << p-2*e << " " << e << " ";
	
	
	// Cache SF data and compute fv1
	for(int i = 0; i < N; i++){
		v[i] = i*2.0*M_PI/N;
		Fr[i] 	= lib_Sch_GSF_Fr_diss(e, p, v[i]) + lib_Sch_GSF_Fr_cons(e, p, v[i]);
		Fphi[i] = lib_Sch_GSF_Fphi_diss(e, p, v[i]) + lib_Sch_GSF_Fphi_cons(e, p, v[i]);
		
		in[i][0] = dw_dchi(p, e, v[i], Fphi[i], Fr[i]);
		in[i][1] = 0.0;
	}

    fftw_execute(*plan); /* repeat as needed */
	for(int i = 0; i < N_out; i++){
   	 	*fv_file << out[i][0]/N << " " << out[i][1]/N << " ";
	}
	*fv_file << endl;
	
	// Compute F1p		
	for(int i = 0; i < N; i++){
		in[i][0] = dp_dchi(p, e, v[i], Fphi[i], Fr[i]);
	}
    fftw_execute(*plan); /* repeat as needed */
	for(int i = 0; i < N_out; i++){
   	 	*Fp_file << out[i][0]/N << " " << out[i][1]/N << " ";
	}
	*Fp_file << endl;
	
	// Compute F1e
	for(int i = 0; i < N; i++){
		in[i][0] = de_dchi(p, e, v[i], Fphi[i], Fr[i]);
	}
    fftw_execute(*plan); /* repeat as needed */
	for(int i = 0; i < N_out; i++){
   	 	*Fe_file << out[i][0]/N << " " << out[i][1]/N << " ";
	}
	*Fe_file << endl;
	
	// The derivatives of F and U w.r.t. p, e, and v are below
	
	// Compute dF1p_dp
	for(int i = 0; i < N; i++){
		dFr_dp[i] = lib_Sch_GSF_dFr_cons_dp(e, p, v[i]) + lib_Sch_GSF_dFr_diss_dp(e, p, v[i]);
		dFr_de[i] = lib_Sch_GSF_dFr_cons_de(e, p, v[i]) + lib_Sch_GSF_dFr_diss_de(e, p, v[i]);
		
		dFphi_dp[i] = lib_Sch_GSF_dFphi_cons_dp(e, p, v[i]) + lib_Sch_GSF_dFphi_diss_dp(e, p, v[i]);
		dFphi_de[i] = lib_Sch_GSF_dFphi_cons_de(e, p, v[i]) + lib_Sch_GSF_dFphi_diss_de(e, p, v[i]);
		
		in[i][0] = dF1p_dp(p, e, v[i], Fphi[i], Fr[i], dFphi_dp[i], dFr_dp[i]);
	}
	fftw_execute(*plan); /* repeat as needed */
	for(int i = 0; i < N_out; i++){
   	 	*dFp_dp_file << out[i][0]/N << " " << out[i][1]/N << " ";
	}
	*dFp_dp_file << endl;
	
	//Compute dF1p_de
	for(int i = 0; i < N; i++){
		in[i][0] = dF1p_de(p, e, v[i], Fphi[i], Fr[i], dFphi_de[i], dFr_de[i]);
	}
	fftw_execute(*plan); /* repeat as needed */
	for(int i = 0; i < N_out; i++){
   	 	*dFp_de_file << out[i][0]/N << " " << out[i][1]/N << " ";
	}
	*dFp_de_file << endl;
	
	//Compute dF1e_dp
	for(int i = 0; i < N; i++){
		in[i][0] = dF1e_dp(p, e, v[i], Fphi[i], Fr[i], dFphi_dp[i], dFr_dp[i]);
	}
	fftw_execute(*plan); /* repeat as needed */
	for(int i = 0; i < N_out; i++){
   	 	*dFe_dp_file << out[i][0]/N << " " << out[i][1]/N << " ";
	}
	*dFe_dp_file << endl;
	
	//Compute dF1e_de
	for(int i = 0; i < N; i++){
		in[i][0] = dF1e_de(p, e, v[i], Fphi[i], Fr[i], dFphi_de[i], dFr_de[i]);
	}
	fftw_execute(*plan); /* repeat as needed */
	for(int i = 0; i < N_out; i++){
   	 	*dFe_de_file << out[i][0]/N << " " << out[i][1]/N << " ";
	}
	*dFe_de_file << endl;
	
	
	
	
	//Compute dV0_dp
	for(int i = 0; i < N; i++){
		in[i][0] = dV0_dp(p, e, v[i]);
	}
	fftw_execute(*plan); /* repeat as needed */
	for(int i = 0; i < N_out; i++){
   	 	*dV0_dp_file << out[i][0]/N << " " << out[i][1]/N << " ";
	}
	*dV0_dp_file << endl;
	
	//Compute dV0_de
	for(int i = 0; i < N; i++){
		in[i][0] = dV0_de(p, e, v[i]);
	}
	fftw_execute(*plan); /* repeat as needed */
	for(int i = 0; i < N_out; i++){
   	 	*dV0_de_file << out[i][0]/N << " " << out[i][1]/N << " ";
	}
	*dV0_de_file << endl;
	
	//Compute dV0_dv
	for(int i = 0; i < N; i++){
		in[i][0] = dV0_dv(p, e, v[i]);
	}
	fftw_execute(*plan); /* repeat as needed */
	for(int i = 0; i < N_out; i++){
   	 	*dV0_dv_file << out[i][0]/N << " " << out[i][1]/N << " ";
	}
	*dV0_dv_file << endl;
	
	
	
	
	
	//Compute dU0_dp
	for(int i = 0; i < N; i++){
		in[i][0] = dU0_dp(p, e, v[i]);
	}
	fftw_execute(*plan); /* repeat as needed */
	for(int i = 0; i < N_out; i++){
   	 	*dU0_dp_file << out[i][0]/N << " " << out[i][1]/N << " ";
	}
	*dU0_dp_file << endl;
	
	//Compute dU0_de
	for(int i = 0; i < N; i++){
		in[i][0] = dU0_de(p, e, v[i]);
	}
	fftw_execute(*plan); /* repeat as needed */
	for(int i = 0; i < N_out; i++){
   	 	*dU0_de_file << out[i][0]/N << " " << out[i][1]/N << " ";
	}
	*dU0_de_file << endl;
	
	//Compute dU0_dv
	for(int i = 0; i < N; i++){
		in[i][0] = dU0_dv(p, e, v[i]);
	}
	fftw_execute(*plan); /* repeat as needed */
	for(int i = 0; i < N_out; i++){
   	 	*dU0_dv_file << out[i][0]/N << " " << out[i][1]/N << " ";
	}
	*dU0_dv_file << endl;
	
}

void construct_tilde_Fs(){
	// Open the Fourier data files
	ifstream fv_file("data/Fourier_fv.dat");
	ifstream Fp_file("data/Fourier_Fp.dat");
	ifstream Fe_file("data/Fourier_Fe.dat");
	
	ifstream dFpdp_file("data/Fourier_dFp_dp.dat");
	ifstream dFpde_file("data/Fourier_dFp_de.dat");
	ifstream dFedp_file("data/Fourier_dFe_dp.dat");
	ifstream dFede_file("data/Fourier_dFe_de.dat");
	
	ifstream dU0dp_file("data/Fourier_dU0_dp.dat");
	ifstream dU0de_file("data/Fourier_dU0_de.dat");
	ifstream dU0dv_file("data/Fourier_dU0_dv.dat");
	ifstream dV0dp_file("data/Fourier_dV0_dp.dat");
	ifstream dV0de_file("data/Fourier_dV0_de.dat");
	ifstream dV0dv_file("data/Fourier_dV0_dv.dat");
	
	// The output file for the Ftilde data
	ofstream Ftilde_file("data/Ftildes.dat");
	Ftilde_file.precision(OUT_PREC);
	
	string fv_string, Fp_string, Fe_string, dFpdp_string, dFpde_string, dFedp_string, dFede_string;
	string dU0dp_string, dU0de_string, dU0dv_string, dV0dp_string, dV0de_string, dV0dv_string;
	double p, e, fv1, Fp1, Fe1, null;
	while(getline(fv_file, fv_string)){			// We assume all the Fourier files are on identical grids
		getline(Fp_file, Fp_string);
		getline(Fe_file, Fe_string);
		getline(dFpdp_file, dFpdp_string);
		getline(dFpde_file, dFpde_string);
		getline(dFedp_file, dFedp_string);
		getline(dFede_file, dFede_string);	
		getline(dU0dp_file, dU0dp_string);
		getline(dU0de_file, dU0de_string);
		getline(dU0dv_file, dU0dv_string);
		getline(dV0dp_file, dV0dp_string);
		getline(dV0de_file, dV0de_string);
		getline(dV0dv_file, dV0dv_string);	
		
		stringstream fv_ss(fv_string);
		stringstream Fp_ss(Fp_string);
		stringstream Fe_ss(Fe_string);
		stringstream dFpdp_ss(dFpdp_string);
		stringstream dFpde_ss(dFpde_string);
		stringstream dFedp_ss(dFedp_string);
		stringstream dFede_ss(dFede_string);
		stringstream dU0dp_ss(dU0dp_string);
		stringstream dU0de_ss(dU0de_string);
		stringstream dU0dv_ss(dU0dv_string);
		stringstream dV0dp_ss(dV0dp_string);
		stringstream dV0de_ss(dV0de_string);
		stringstream dV0dv_ss(dV0dv_string);
		
		Fp_ss >> p >> e >> Fp1 >> null;
		Fe_ss >> p >> e >> Fe1 >> null;
		fv_ss >> p >> e >> fv1 >> null;
		fv1 = -fv1;			// Note the sign change as xi = chi - chi0
		
		dFpdp_ss >> p >> e >> null >> null;
		dFpde_ss >> p >> e >> null >> null;
		dFedp_ss >> p >> e >> null >> null;
		dFede_ss >> p >> e >> null >> null;
		
		dU0dp_ss >> p >> e >> null >> null;
		dU0de_ss >> p >> e >> null >> null;
		dU0dv_ss >> p >> e >> null >> null;
		dV0dp_ss >> p >> e >> null >> null;
		dV0de_ss >> p >> e >> null >> null;
		dV0dv_ss >> p >> e >> null >> null;
		
		double re, im;
		double Fp2 = 0, Fe2 = 0, U1 = 0, V1 = 0, Xv1 = 0, Yp1 = 0, Ye1 = 0;
		Complex Fpk, Fek, fvk, dFpdpk, dFpdek, dFedpk, dFedek;
		Complex dU0dpk, dU0dek, dU0dvk, dV0dpk, dV0dek, dV0dvk;
		for(int k = 1; k < 10; k++){
			Fp_ss >> re >> im;
			Fpk = re + 1i*im;
			
			Fe_ss >> re >> im;
			Fek = re + 1i*im;
			
			fv_ss >> re >> im;
			fvk = -re - 1i*im;			// Note the sign change as xi = chi - chi0
			
			dFpdp_ss >> re >> im;
			dFpdpk = re + 1i*im;
			
			dFpde_ss >> re >> im;
			dFpdek = re + 1i*im;
			
			dFedp_ss >> re >> im;
			dFedpk = re + 1i*im;
			
			dFede_ss >> re >> im;
			dFedek = re + 1i*im;
			
			dU0dp_ss >> re >> im;
			dU0dpk = re + 1i*im;
			
			dU0de_ss >> re >> im;
			dU0dek = re + 1i*im;
			
			dU0dv_ss >> re >> im;
			dU0dvk = re + 1i*im;
			
			dV0dp_ss >> re >> im;
			dV0dpk = re + 1i*im;
			
			dV0de_ss >> re >> im;
			dV0dek = re + 1i*im;
			
			dV0dv_ss >> re >> im;
			dV0dvk = re + 1i*im;
			
			Fp2 += 2.*( 1.i/((double)k)*(dFpdpk*conj(Fpk) + dFpdek*conj(Fek)) - Fpk*conj(fvk) ).real();
			Fe2 += 2.*( 1.i/((double)k)*(dFedpk*conj(Fpk) + dFedek*conj(Fek)) - Fek*conj(fvk) ).real();
			
			U1 += 2.*(dU0dpk*conj(Fpk) + dU0dek*conj(Fek) + dU0dvk*conj(fvk)).real();
			V1 += 2.*(dV0dpk*conj(Fpk) + dV0dek*conj(Fek) + dV0dvk*conj(fvk)).real();
			
			// The below is for the initial condition matching and assumes that v0 = 0
			Xv1 += 2.*(1.i/((double)k)*fvk).real();
			Yp1 += 2.*(1.i/((double)k)*Fpk).real();
			Ye1 += 2.*(1.i/((double)k)*Fek).real();			
		}
		Ftilde_file  << p << " " << e << " " << Fp1 << " " << Fe1 << " " << fv1 << " " << Fp2 << " " << Fe2 << " " << U1 << " " << V1 << " " << Xv1 << " " << Yp1 << " " << Ye1 << endl;
	};
	
	
	fv_file.close();
	Fp_file.close();
	Fe_file.close();
	Ftilde_file.close();
}


// Code for explicit integration of the {p, e, v, t, phi} equations using the full self-force

int osc_eqs (double chi, const double y[], double f[], void *params){
	
	double p = y[0];
	double e = y[1];
	double chi0 = y[2];
	double v = chi - chi0;
	
	double Fr 	= q*(lib_Sch_GSF_Fr_diss(e, p, v) + lib_Sch_GSF_Fr_cons(e, p, v));
	double Fphi = q*(lib_Sch_GSF_Fphi_diss(e, p, v) + lib_Sch_GSF_Fphi_cons(e, p, v));

	f[0] = dp_dchi(p, e, v, Fphi, Fr);
	f[1] = de_dchi(p, e, v, Fphi, Fr);
	f[2] = dw_dchi(p, e, v, Fphi, Fr);
	f[3] = dt_dchi(p, e, v);
	f[4] = dphi_dchi(p, e, v);	
	
	if(p-6-2*e > Y_MIN) return GSL_SUCCESS;
	else return GSL_SUCCESS + 1;
}

void integrate_osc_eqs(double p0, double e0){
	
	double chi = 0;
	
	// y[0] = p, y[1] = e, y[2] = chi0, y[3] = t, y[4] = phi
	double y[5] = {p0, e0, 0, 0, 0};
	
	// Output the initial parameters
	fout << "# Format: chi    p    e    chi0    t    phi" << endl;
	fout << chi << " " << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << " " << y[4] << " " << endl;


	int dense_output;
	try{
		dense_output 		= cfg.lookup("Dense_output");		
 	}catch(const SettingNotFoundException &nfex){
   		cerr << "'Dense_output' setting missing from configuration file." << endl; exit(0);
 	}


	if(dense_output == 1){
		
		double n_per_orbit;
		try{
			n_per_orbit 		= cfg.lookup("n_per_orbit");		
	 	}catch(const SettingNotFoundException &nfex){
	   		cerr << "'n_per_orbit' setting missing from configuration file." << endl; exit(0);
	 	}
		
	    gsl_odeiv2_system sys 	= {osc_eqs, NULL, 5, NULL};	
	    gsl_odeiv2_driver *d 	= gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, NUM_TOL, NUM_TOL, 0.0);
		
		double chi_i = 0;
		double delta_chi = 2.0*M_PI/n_per_orbit;
		
		while(1){
			chi_i += delta_chi;
			int status = gsl_odeiv2_driver_apply (d, &chi, chi_i, y);

			if (status != GSL_SUCCESS){
				printf ("error, return value=%d\n", status);
				break;
			}

			fout << chi << " " << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << " " << y[4] << " " << endl;
			
		}

		gsl_odeiv2_driver_free (d);
	}else{		// Sparse output
	    const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk8pd;

	    gsl_odeiv2_step *s 		= gsl_odeiv2_step_alloc (T, 5);
	    gsl_odeiv2_control *c 	= gsl_odeiv2_control_y_new (NUM_TOL, NUM_TOL);
	    gsl_odeiv2_evolve *e 	= gsl_odeiv2_evolve_alloc (5);

	    gsl_odeiv2_system sys = {osc_eqs, NULL, 5, NULL};

	    double chi_max = CHI_MAX;
		double h = 1e-6;

	    high_resolution_clock::time_point t1 = high_resolution_clock::now();

	    while (chi < chi_max) {
	        int status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &chi, chi_max, &h, y);
			if (status != GSL_SUCCESS){
				//printf ("error, return value=%d\n", status);
				break;
			}

			fout << chi << " " << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << " " << y[4] << " " << endl;
			
	      }
		  high_resolution_clock::time_point t2 = high_resolution_clock::now();
		  
		  duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
		  fout << "# Computing the inspiral took: " << time_span.count() << " seconds." << endl;
		  
		  
		  
		  gsl_odeiv2_evolve_free (e);
		  gsl_odeiv2_control_free (c);
		  gsl_odeiv2_step_free (s);
	}
	
}
//-------------------------------------------------------------
//  -f0 additions for full self-force with high eccentricity
//-------------------------------------------------------------

// interpolants for each self-force component (global variables)
vector<Interpolant*> F_interp_cos[3], F_interp_sin[3];

double insp_at(double chi, double e, double p){
return -(sqrt(-4*e*e+(-2+p)*(-2+p))*(3+e*e-p)*(6+2*e*e-p)*p*pow(1+e*cos(chi),2)*(2-p+2*e*cos(chi)));}

double insp_aphi(double chi, double e, double p){
return (1-e*e)*(4*e*e+(-6+p)*(-2+p))*(3+e*e-p)*p*p*sqrt(p);}

double insp_bt(double chi, double e, double p){
return -2*e*sqrt(-4*e*e+(-2+p)*(-2+p))*(-3-e*e+p)*(-2+p-2*e*cos(chi))*p*p*pow(1+e*cos(chi),2);}

double insp_bphi(double chi, double e, double p){
return 2*e*(-4+p)*(-4+p)*p*p*p*sqrt(p)*(-3-e*e+p);}

double insp_q(double chi, double e, double p){
return e*(-6-2*e+p)*(-6+2*e+p)*sqrt(-6+p-2*e*cos(chi))*pow(1+e*cos(chi),4);}

double getdpdchi(double e, double p, double v, double Ft, double Fphi){
return (insp_bt(v,e,p)*Ft+insp_bphi(v,e,p)*Fphi)/insp_q(v,e,p);}

double getdedchi(double e, double p, double v, double Ft, double Fphi){
return (insp_at(v,e,p)*Ft+insp_aphi(v,e,p)*Fphi)/insp_q(v,e,p);}

double getF(double e, double p, double v, vector<Interpolant*> F_interp_cos_comp, vector<Interpolant*> F_interp_sin_comp){
double y = p-2*e;
int n;
int n_max = F_interp_cos_comp.size()-1;
double F = F_interp_cos_comp[0]->eval(y,e);
for(n=1; n<=n_max; n++)	F += F_interp_cos_comp[n]->eval(y,e)*cos(n*v) + F_interp_sin_comp[n]->eval(y,e)*sin(n*v);
return F;}

int osc_eqs_default(double chi, const double y[], double f[], void *params){
	double p = y[0];
	double e = y[1];
	double chi0 = y[2];
	double v = chi - chi0;
	
	//double Fr = q*(lib_Sch_GSF_Fr_diss(e, p, v) + lib_Sch_GSF_Fr_cons(e, p, v));
	//double Fphi = q*(lib_Sch_GSF_Fphi_diss(e, p, v) + lib_Sch_GSF_Fphi_cons(e, p, v));
	double Ft = q*getF(e, p, v, F_interp_cos[0], F_interp_sin[0]);
    double Fr 	= q*getF(e, p, v, F_interp_cos[1], F_interp_sin[1]);
	double Fphi = q*getF(e, p, v, F_interp_cos[2], F_interp_sin[2]);

	//f[0] = dp_dchi(p, e, v, Fphi, Fr);
	//f[1] = de_dchi(p, e, v, Fphi, Fr);
	f[0] = getdpdchi(e, p, v, Ft, Fphi);
    f[1] = getdedchi(e, p, v, Ft, Fphi);
	f[2] = dw_dchi(p, e, v, Fphi, Fr);
	f[3] = dt_dchi(p, e, v);
	f[4] = dphi_dchi(p, e, v);	
	
	if(p-6-2*e > Y_MIN) return GSL_SUCCESS;
	else return GSL_SUCCESS + 1;
}

void integrate_osc_eqs_default(double p0, double e0){
	// import self-force data
	ifstream inFile;
	double y_col, e_col;
	vector<double> ys, es;
	string input_names[3], line, temp;
	input_names[0] = "GSF_model_data/Ft.dat";
	input_names[1] = "GSF_model_data/Fr.dat";
	input_names[2] = "GSF_model_data/Fphi.dat";
	complex<double> cos_sin; 
	int trp, n_max, test, n;
	int j = 0;
	for(trp=0;trp<3;trp++)
		{
		inFile.open(input_names[trp]);
		if (!inFile)
			{
			cout << "Unable to open file \n";
			exit(1); // terminate with error
			}
		// count the number of columns and advance past 1st row
		getline(inFile,line);
		istringstream ss(line);
		ss >> temp;
		ss >> temp;
		while(ss >> test) j++;
		if(trp==0) n_max = j-1;
		vector<double> F_data_cos[n_max+1], F_data_sin[n_max+1];
				
		// loop through each orbit (rows in data file)
		while(getline(inFile,line))
			{
			istringstream ss(line);

			ss >> y_col;    // extracts 1st col, which is y = p-2e
			ss >> e_col;    // extracts 2nd col, which is e
			// assumes the 3 files have the same orbits
			if(trp==0)
				{
				ys.push_back(y_col);
				es.push_back(e_col);
				}
			
			for(n=0; n<=n_max; n++) // loop through each n (columns in data file)
				{
				ss >> cos_sin;
				// no need for another vector (n_re & n_im)
				// can access the components of N_re and N_imag directly
				F_data_cos[n].push_back(cos_sin.real());
				F_data_sin[n].push_back(cos_sin.imag());
				}
			}
		inFile.close();
		// interpolate self-force data
		for(n=0; n<=n_max; n++)
			{
			F_interp_cos[trp].push_back(new Interpolant(ys, es, F_data_cos[n]));
			F_interp_sin[trp].push_back(new Interpolant(ys, es, F_data_sin[n]));
			}
		}
	
	double chi = 0;
	
	// y[0] = p, y[1] = e, y[2] = chi0, y[3] = t, y[4] = phi
	double y[5] = {p0, e0, 0, 0, 0};
	
	// Output the initial parameters
	fout << "# Format: chi    p    e    chi0    t    phi" << endl;
	fout << chi << " " << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << " " << y[4] << " " << endl;


	int dense_output;
	try{
		dense_output 		= cfg.lookup("Dense_output");		
 	}catch(const SettingNotFoundException &nfex){
   		cerr << "'Dense_output' setting missing from configuration file." << endl; exit(0);
 	}

	if(dense_output == 1){
		
		double n_per_orbit;
		try{
			n_per_orbit 		= cfg.lookup("n_per_orbit");		
	 	}catch(const SettingNotFoundException &nfex){
	   		cerr << "'n_per_orbit' setting missing from configuration file." << endl; exit(0);
		}
	
	    gsl_odeiv2_system sys 	= {osc_eqs_default, NULL, 5, NULL};	
	    gsl_odeiv2_driver *d 	= gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, NUM_TOL, NUM_TOL, 0.0);
		
		double chi_i = 0;
		double delta_chi = 2.0*M_PI/n_per_orbit;
		
		while(1){
			chi_i += delta_chi;
			int status = gsl_odeiv2_driver_apply (d, &chi, chi_i, y);

			if (status != GSL_SUCCESS){
				printf ("error, return value=%d\n", status);
				break;
			}

			fout << chi << " " << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << " " << y[4] << " " << endl;
			
		}

		gsl_odeiv2_driver_free (d);
	}else{		// Sparse output
	    const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk8pd;

	    gsl_odeiv2_step *s 		= gsl_odeiv2_step_alloc (T, 5);
	    gsl_odeiv2_control *c 	= gsl_odeiv2_control_y_new (NUM_TOL, NUM_TOL);
	    gsl_odeiv2_evolve *e 	= gsl_odeiv2_evolve_alloc (5);

	    gsl_odeiv2_system sys = {osc_eqs_default, NULL, 5, NULL};

	    double chi_max = CHI_MAX;
		double h = 1e-6;

	    high_resolution_clock::time_point t1 = high_resolution_clock::now();

	    while (chi < chi_max) {
	        int status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &chi, chi_max, &h, y);
			if (status != GSL_SUCCESS){
				//printf ("error, return value=%d\n", status);
				break;
			}

			fout << chi << " " << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << " " << y[4] << " " << endl;
			
	      }
		  high_resolution_clock::time_point t2 = high_resolution_clock::now();
		  
		  duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
		  fout << "# Computing the inspiral took: " << time_span.count() << " seconds." << endl;
		  
		  
		  
		  gsl_odeiv2_evolve_free (e);
		  gsl_odeiv2_control_free (c);
		  gsl_odeiv2_step_free (s);
	}
	
	// free memory from interpolants
	for(trp=0; trp<3; trp++)
		{
		for(n=0; n<=n_max; n++)
			{
			delete F_interp_cos[trp][n];
			delete F_interp_sin[trp][n];
			}
		}
	
}


// Waveform generation functions below

void compute_waveform(string insp_filename, string out_filename){

	// Check if the associated inspiral trajectory file exists
	ifstream insp(insp_filename);
	if(!insp){
		cout << "Inspiral file: " << insp_filename << " does not exist." << endl;
		exit(0);
	}

	// Load and interpolate the inspiral trajectory
	cout << "Loading and interpolating inspiral trajectory data"  << endl;
	
	// FIXME only load data up to t_max (rather than the entire phase space trajectory as is currently done)	
	string insp_string;
	vector<double> ps, es, vs, ts, phis;
	double chi, p, e, v, t, phi;
	int test = 0;
	while(getline(insp, insp_string)){
				
		if(insp_string.at(0) == '#') continue;
		
		stringstream insp_ss(insp_string);
		
		insp_ss >> chi >> p >> e >> v >> t >> phi;
						
		if(mode == WAVEFORM_FULL || mode == T_WAVEFORM_FULL)
			v = chi-v;
		
		if(test == 0) test = 1;
		else if(v <= vs.back()) break;
		
		vs.push_back(v);
		
		ps.push_back(p);
		es.push_back(e);
		ts.push_back(t);
		phis.push_back(phi);
	}
	
	Interpolant p_interp(vs, ps);
	Interpolant e_interp(vs, es);
	Interpolant t_interp(vs, ts);
	Interpolant phi_interp(vs, phis);
		
	// Read in the parameters and give an error if they do not exist
	double M_solar, Deltat_sec;
	int i_max;
	try{
		M_solar 		= cfg.lookup("M_solar");		// Mass of the primary in solar masses}
 	}catch(const SettingNotFoundException &nfex){
   		cerr << "'M_solar' setting missing from configuration file." << endl; exit(0);
 	}
	
	try{
		Deltat_sec 	= cfg.lookup("Deltat_sec");		// Time step in seconds
	}catch(const SettingNotFoundException &nfex){
   		cerr << "'Delta_sec' setting missing from configuration file." << endl; exit(0);
 	}
	
	try{
		i_max 	= cfg.lookup("i_max");				// The number of time steps to take
	}catch(const SettingNotFoundException &nfex){
   		cerr << "'i_max' setting missing from configuration file." << endl; exit(0);
 	}
			 	
	double Msolar_sec	= 4.9253908380897e-6;	// Solar mass in seconds				
	double M 			= M_solar * Msolar_sec;	// Mass of the primary in seconds
	double Deltat		= Deltat_sec/M;			// Time step in units of M	
	
	cout << "Computing the waveform" << endl;
	fout.open(out_filename);
	
	fout << "# M     = " << M_solar << " [Solar Masses]" << endl;
	fout << "# q     = " << q << endl;
	fout << "# Δt    = " << Deltat_sec << " [Seconds]" << endl;
	fout << "# t_end = " << Deltat_sec*i_max/60./60./24./30. << " [30 day months]" << endl;
	
	
	// Resample t more densely
	// FIXME only resample up to t_max
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
	double t_max = i_max*Deltat;
	vector<double> t_dense, v_dense;
	
	v = vs[0];
		
	while(v < vs[vs.size() - 1]){
	
		p = p_interp.eval(v);
		e = e_interp.eval(v);
			
		v_dense.push_back(v);
		if(mode == WAVEFORM_NIT || mode == T_WAVEFORM_NIT) t = t_interp.eval(v) - U0(p,e,v);
		else t = t_interp.eval(v);
		
		t_dense.push_back(t);
		
		if(t > t_max) break;	
		
		v += 2.0*M_PI/10.;
	
	}
	 // fprintf(stderr, "test 1\n");
	Interpolant v_interp(t_dense, v_dense);
	high_resolution_clock::time_point t2 = high_resolution_clock::now();

	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	fout << "# Resampling t took: " << time_span.count() << " seconds." << endl;

	if(mode == WAVEFORM_NIT){
		fout << "# Format: t h+ h×" << endl;
	} else if (mode == T_WAVEFORM_NIT){
		fout << "# Format: t Re(HTeuk) Im(HTeuk)" << endl;
	}

	// --------------------- prep for Teukolsky waveform calculations ---------------------------- //

	//if (mode  == T_WAVEFORM_NIT) // This runs the importing and creation of interpolants in the T_WAVEFORM case
	//{
		// First pre-calculate the n_modes and then execute the complex amp code for creating the interpolants:

		vector<int> n_modes;

		for (int i = -40; i <= 40; i++)
    	{
			n_modes.push_back(i);
    	}

		// ------------------------- Interpolant of Complex Amplitudes ------------------------- //

        //----------------------------
        //  read Almn data from file
        //----------------------------
        ifstream inFile;

        // initialize variables outside loops (speeds up by avoiding extra initialization cost)
        string line, temp;
        double y_col, e_col, stest;
        int iter, input, num_cols;
        vector<double> ys, Es; //y = p-2*e
        complex<double> Almn; 

        int col_count = 2; // This will be the total number of columns
        inFile.open("GSF_model_data/A_l2_m2.dat"); 
        if (!inFile)
            {
            cout << "Unable to open file \n";
            exit(1); // terminate with error
            }
        else
            { // count the number of columns and advance past 1st row
            // (don't need to erase 1st row later)
            getline(inFile,line);
            istringstream ss(line);
            ss >> temp;
            ss >> temp;
            while(ss >> stest) col_count++;
            }
            
        // initialize 2 vectors (real and imaginary) for each n
        // used an array of vectors (instead of a vector of vectors) because we know the number of columns at this point
        vector<double> N_re[col_count-2], N_im[col_count-2]; // 1st index is each n, 2nd index is each orbit

        // loop through each orbit (rows in data file)
        while(getline(inFile,line))
            {
            istringstream ss(line);


            // I am using y instead of x to match the variable names in the existing code
            ss >> y_col;    // extracts 1st col, which is y = p-2e
            ss >> e_col;    // extracts 2nd col, which is e

            ys.push_back(y_col);
            Es.push_back(e_col);
            
            for(iter=0; iter<col_count-2; iter++) // loop through each n (columns in data file)
                {
                ss >> Almn;
                // no need for another vector (n_re & n_im)
                // can access the components of N_re and N_imag directly
                N_re[iter].push_back(Almn.real());
                N_im[iter].push_back(Almn.imag());
                }
            }
        inFile.close();  
		
		// Create the interpolants:
		vector<Interpolant*> N_Interp_re, N_Interp_im;
        vector<vector<double>> A_re, A_im;

            for (int i = 0; i < col_count - 2; i++)
            {
                // Interpolant part:
                N_Interp_re.push_back(new Interpolant(ys, Es, N_re[i]));
                N_Interp_im.push_back(new Interpolant(ys, Es, N_im[i]));
			}
	//}
	 // As of now we have a vector of interpolant objects which are callable to find the Almn's
	
	t1 = high_resolution_clock::now();

	// ------------------------- waveform computation begins here --------------------------------- //
	Complex h, HTeuk;
	vector<double> hplus, hcross, HTeuk_re, HTeuk_im;

	// The new parameter initializations:
	double Ar, Ai, tt,t_b, w_phi, w_r;
	double p_before,e_before,v_before,phi_before,t_before; // needed for delta quantities in omegas

	double m = 2;

	int i_max_test = floor(t_dense.back()/Deltat)-1;
	if(i_max > i_max_test) i_max = i_max_test;

	// For loop which goes for a very large number of iterations
	for(int i = 0; i <= i_max; i ++){
		t 		= i*Deltat; // t is normal t in this loop, so define tt as t tilde
		v 	= v_interp.eval(t);

	// We calculate the individual quantities usedat each iteration rather than pulling from them all stored somewhere:	
		p = p_interp.eval(v);
		e = e_interp.eval(v);
		phi = phi_interp.eval(v);

	// Need phi for -w -n, but keep phi-tilde for -t -n
		if(mode == WAVEFORM_NIT)
			phi -= V0(p,e,v);						
		
	// We calculate the complex amplitudes of the waveform using the corresponding equation:
		if(mode == WAVEFORM_NIT || mode == WAVEFORM_FULL){

			h = waveform_h(p, e, v, phi, 0., 0.);
		
			hplus.push_back(h.real());
			hcross.push_back(h.imag());

		} else if (mode == T_WAVEFORM_NIT){
			if (i == 0){
				continue; // if we are on the initial iteration, pass to the next (for delta quantities)
			}
			else{
				t_b = (i-1)*Deltat;
				v_before = v_interp.eval(t_b);
				t_before = t_interp.eval(v_before); // t-tilde

				p_before = p_interp.eval(v_before);
				e_before = e_interp.eval(v_before);

				phi_before = phi_interp.eval(v_before); // again keep phi-tilde for -t -n

				//Calculate the omega_phi and omega_r for this iteration:

				tt = t_interp.eval(v); // this is t tilde since we are not applying the inverse transform

				w_r = (v - v_before)/(tt - t_before);
				w_phi = (phi-phi_before)/(tt-t_before);

				// Calculate the un-tilded t:

				

				// loop to sum all values over the n_modes to find each HTeuk value:

				// Calculate the corresponding Almn for this iteration: (loop these in such that it will 
				// generate each appropriate a during the summation process)

				Complex sum;

				for (int it = 0; it < N_Interp_re.size(); it++){
					Ar = N_Interp_re[it]->eval(p-2*e, e);
					Ai = N_Interp_im[it]->eval(p-2*e, e);

					double power = (m * phi) + (n_modes[it] * v) + (m * w_phi + n_modes[it] * w_r)*(t-tt); // t is normal t and tt is t tilde 
    				Complex coeff = Complex(Ar , Ai);
					sum += coeff * (Cos(power)-Complex(0,Sin(power))); // use Euler's theorem (since power is negative cos-isin)
				}

				// calculate the waveform stuff for Teukolsky:
				HTeuk = sum;

				HTeuk_re.push_back(HTeuk.real());
				HTeuk_im.push_back(HTeuk.imag());
			}
		} else if (mode == T_WAVEFORM_FULL){
			
			// ------------------------ j(v) and j related calculations --------------------------- //
			if (i == 0){
				v_before = v;
				continue;
			} else if(i != 0){
			
			double j = floor(v/(2*M_PI)); // j as a function of v

			if(j == 0){
				j = 1;
			}
			
			double vj = 2*M_PI*j;
			double j_before = floor(v_before/(2*M_PI)); // j as a function of v
			double vj_before = 2*M_PI*j_before;
			double tj_before = t_interp.eval(vj_before); 
			double phij_before = phi_interp.eval(vj_before);

			double w_r = 0;
			double w_phi = 0;

			//cout << j << endl;

			// keep vj of previous j delta with the last value of vj not for each iteration of vj(maybe calculate frequencies here?)
			// when j turns, keep with the omega at the turniing point fromj to j+1
			//Conditional statement to update j and vj if j_before changes
			if (j != j_before){
			
				vj = 2*M_PI*j;
				j_before = j; 
				double tj = t_interp.eval(vj);
				double phij = phi_interp.eval(vj);
				w_r = (vj - vj_before)/(tj - tj_before);
				w_phi = (phij - phij_before)/(tj - tj_before);
				}
			
			double tj = t_interp.eval(vj);
			double phij = phi_interp.eval(vj);
			
			// ---------------------- Calculaint HTeuk using 3.4,3.5,3.6 and 3.13 ----------------- // 
			if (i < 40){
			cout << w_r << " " << w_phi << endl;
			cout << j << endl;
			}
				Complex sum;

				for (int it = 0; it < N_Interp_re.size(); it++){
					Ar = N_Interp_re[it]->eval(p-2*e, e);
					Ai = N_Interp_im[it]->eval(p-2*e, e);

					double power = (m * phij) + (n_modes[it] * vj) + (m * w_phi + n_modes[it] * w_r)*(t-tj); // t is normal t and tt is t tilde 
    				Complex coeff = Complex(Ar , Ai);
					sum += coeff * (Cos(power)-Complex(0,Sin(power))); // use Euler's theorem (since power is negative cos-isin)
				}

				// calculate the waveform stuff for Teukolsky:
				HTeuk = sum;

				HTeuk_re.push_back(HTeuk.real());
				HTeuk_im.push_back(HTeuk.imag());

			}


		}
	} // This is the end of the loop that is used to create the complex waveform values
	
	t2 = high_resolution_clock::now();

	time_span = duration_cast<duration<double>>(t2 - t1);
	fout << "# Computing the waveform took: " << time_span.count() << " seconds." << endl;
	
	// ----------------- Output the waveform data --------------------------
	cout << "Outputting waveform to " << out_filename << endl;
	double t_sec;
	for(int i = 0; i <= i_max; i ++){
		
		t_sec 	= i*Deltat_sec;
		
		if (mode  == WAVEFORM_NIT || mode == WAVEFORM_FULL){
			fout << t_sec << " " << hplus[i] << " " << hcross[i] << endl;\

		} else if (mode == T_WAVEFORM_NIT || mode == T_WAVEFORM_FULL){
			fout << t_sec << " " << HTeuk_re[i] << " " << HTeuk_im[i] << endl;
		}
	}
	
}

// Definitions for using CForm'ed output from Mathematica
#define Sin(x)          (sin((double)(x)))
#define Cos(x)          (cos((double)(x)))
#define Power(x, y)     (pow((double)(x), (double)(y)))
#define Sqrt(x)         (sqrt((double)(x)))
#define Conjugate(x)	(conj(x))
#define Pi				M_PI

// Returns the quadrupolar waveform, Z = Cos Theta and Theta and Phi are the angles between the detector and the source
Complex waveform_h(double p, double e, double v, double phi, double Z, double Phi){
	
	Complex XX = (5*Power(e,4) + Power(e,2)*(52 - 10*p) - 8*(-2 + p)*p)*Sqrt(-6 + p - 2*e*Cos(v)) + 
   e*(48 - Power(e,2)*(-56 + p) + 4*(4 - 3*p)*p)*Cos(v)*Sqrt(-6 + p - 2*e*Cos(v)) + 
   2*Power(e,2)*(30 + 4*Power(e,2) + (5 - 2*p)*p)*Sqrt(-6 + p - 2*e*Cos(v))*Cos(2*v) + Power(e,3)*(24 + p)*Sqrt(-6 + p - 2*e*Cos(v))*Cos(3*v) + 
   3*Power(e,4)*Sqrt(-6 + p - 2*e*Cos(v))*Cos(4*v) - Complex(0,8)*e*Sqrt(p)*(-6 + p - 2*e*Cos(v))*(-1 + p - e*Cos(v))*(1 + e*Cos(v))*Sin(v);
	
	return ((exp(2.*Complex(0,1.)*(-phi + Phi))*XX*Power(1 - Z,2) + exp(2.*Complex(0,1.)*(phi - Phi))*Power(1 + Z,2)*Conjugate(XX))*(-2 + p - 2*e*Cos(v)))/
    (4.*(-4*Power(e,2) + Power(-2 + p,2))*Power(p,1.5)*Sqrt(-6 + p - 2*e*Cos(v))) - 
   (e*(-1 + Power(Z,2))*(2 - p + 2*e*Cos(v))*((Power(e,2)*(56 - 13*p) + 4*(12 - 8*p + Power(p,2)))*Cos(v) + 
        e*(52 + 5*Power(e,2) - 34*p + 4*Power(p,2) + 2*(30 + 4*Power(e,2) - 7*p)*Cos(2*v) - 3*e*(-8 + p)*Cos(3*v) + 3*Power(e,2)*Cos(4*v))))/
    (2.*(-4*Power(e,2) + Power(-2 + p,2))*Power(p,2));
}