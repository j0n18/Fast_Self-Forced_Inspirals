// ---------------------- Get Error Data from NIT and Full inspirals ------------------------------
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

#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include <vector>
#include <algorithm>
#include <chrono>
#include <math.h>
#include <cstring>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_matrix.h>

#include <fftw3.h>

using namespace std;
using namespace std::literals;
using namespace std::chrono;

typedef complex<double> Complex;

// Used to set the mode the code runs in
#define FULL_INSPIRAL  			0
#define FULL_INSPIRAL_DEFAULT	1
#define NIT_INSPIRAL   			2
#define NIT_INSPIRAL_DEFAULT	3
#define DECOMPOSE      			4
#define CONSTRUCT_Fs   			5
#define WAVEFORM_FULL  			6
#define WAVEFORM_NIT   			7
#define T_WAVEFORM_FULL  		8
#define T_WAVEFORM_NIT   		9


// -------------------------------- Include Statements --------------------------------------------
#include <Interpolant.h>
 // #include "../include/NIT_inspiral.h"
#include <libconfig.h++>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_sf_ellint.h> // This gives the elliptic functions
#include <math.h>
#include <iomanip>

// ComplexAmpTest include statements:
#include <fstream>
#include <string>
#include <iostream>
#include <complex>
 // #include <iomanip>
#include <algorithm>

#define OUT_PREC 12
#define Sin(x)          (sin((double)(x)))
#define Cos(x)          (cos((double)(x)))
#define Power(x, y)     (pow((double)(x), (double)(y)))
#define Sqrt(x)         (sqrt((double)(x)))
#define Pi				M_PI
#define Conjugate(x)	(conj(x))

using namespace libconfig;
Config cfg;



// precision of output
void GetErrors(vector<string> insp_filenames, string out_filename);
double Tr(double p, double e, double v);
double Z0t(double p, double e, double v);
double Z0phi(double p, double e, double v);
double Phir(double p, double e);
double EllipticK(double k);
double EllipticF(double phi, double k);
double EllipticE(double k);
double EllipticEIncomp(double phi, double k);
double EllipticPi(double n, double k);
double EllipticPiIncomp(double n, double phi, double k);

// ------------------------------------------------------------------------- //
/*
 * Erase First Occurrence of given substring from main string.
 */
void eraseSubStr(std::string & mainStr, const std::string & toErase)
{
    // Search for the substring in string
    size_t pos = mainStr.find(toErase);
    if (pos != std::string::npos)
    {
        // If found then erase it from string
        mainStr.erase(pos, toErase.length());
    }
}
// ------------------------------------------------------------------------- //


// ----------------------------- Generate the Inspiral Data ---------------------------------------
int main(int argc, char* argv[])
{

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


string fileloc, NITfile0, Fullfile0, NITfile, Fullfile, outNITfile0, outFullfile0, outNITfile, outFullfile,NITfile1, Fullfile1;
vector<string> NITStringsList, FullStringsList, outNITStringsList, outFullStringsList;

fileloc = "/mnt/c/Users/JonnyM/Desktop/GravityWaves/Fast_Self-Forced_Inspirals/";

 // Initial strings
NITfile0 = "NIT_inspiral -n0 10 0.1 0.001";
Fullfile0 = "NIT_inspiral -f0 10 0.1 0.001"; // The first part should be NIT_inspiral for directory reasons. -f0 signifies Full.

NITfile1 = "NIT_inspiral -n 10 0.1 0.001";
Fullfile1 = "NIT_inspiral -f 10 0.1 0.001";

double dp = 1; // make 1
double de = 0.08;
int pmax = 11; // make 14 // p > 11 does not work for -n and (reliably) for -f
int pmin = 10; // make 8
double emax = 0.20; // make 0.75; This will be rounded to 0.7
double emin = 0.1; // make 0.1
double floatValue, in_val, in_val1,in_val2;

string half = ".5"; // not going to need
string in_str, sub, sub1;
string newNIT, newFull,newNIT_e, newFull_e, newoutNIT, newoutFull,newoutNIT_e, newoutFull_e;
Fullfile = Fullfile0;
NITfile = NITfile0;

vector<string> OutputFiles;
string outputFile, newoutputFile,newoutputFile_e;

string outputFile0 = "output/ErrorData_p10_e0.1_q0.001.dat";

outputFile = outputFile0;

outNITfile0 = "output/Inspiral_NIT_p12_e0.7_q0.001.dat";
outFullfile0 = "output/Inspiral_Full_p12_e0.7_q0.001.dat";

outFullfile = outFullfile0;
outNITfile = outNITfile0;
int reg_counter = 0;

// -----------------------------
double reg = 0;
// -----------------------------

if (reg == 0){

	Fullfile = Fullfile1;
	NITfile = NITfile1;

} else
{
	Fullfile = Fullfile0;
	NITfile = NITfile0;
}

for ( int i = 1; i <= floor((pmax - pmin)*(1/dp)); i++) // eventually include evething upt to and including GetErrors(inputStrings,outputFile); in this loop
{
	in_val = pmin + i*dp;
	in_val1 = in_val;
	in_str = to_string(in_val);

	sub = in_str.substr(0,4);

	if (in_val < 10){
		sub1 = in_str.substr(0,1);	
	}
	else{
		sub1 = in_str.substr(0,2);
	}

	if (reg == 0){

		newNIT = NITfile.replace(16,2,sub);
		newFull = Fullfile.replace(16,2,sub);

	}
	
	else{

		newNIT = NITfile.replace(17,2,sub);
		newFull = Fullfile.replace(17,2,sub);

	}

		newoutNIT = outNITfile.replace(21,2,sub1);
		newoutFull = outFullfile.replace(22,2,sub1);
		newoutputFile = outputFile.replace(18,2,sub1);	

	for ( int j = 0; j <= floor((emax - emin)*(1/de)); j++) // eventually include evething upt to and including GetErrors(inputStrings,outputFile); in this loop
	{	
		in_val = emin + j*de;
		in_val2 = in_val;
		//cout << in_val << endl;
		in_str = to_string(in_val);

		if (reg == 0){
			sub = in_str.substr(0,4);
		}
		else{
			sub = in_str.substr(0,3);
		}

		if (reg == 0){

			newNIT_e = newNIT.replace(21,4,sub);
			newFull_e = newFull.replace(21,4,sub);

		if (reg_counter == 0){
			reg_counter += 1;
			newNIT_e = newNIT.insert(25," ");
			newFull_e = newFull.insert(25," ");	
		}

		}
		else{

			newNIT_e = newNIT.replace(22,3,sub);
			newFull_e = newFull.replace(22,3,sub);

		}

			newoutNIT_e = newoutNIT.replace(25,3,sub);
			newoutFull_e = newoutFull.replace(26,3,sub);
			newoutputFile_e = newoutputFile.replace(22,3,sub);

		if (reg == 0){
			newoutNIT_e.erase(28,1);
			newoutFull_e.erase(29,1);
			newoutputFile_e.erase(25,1);
		}

		if (reg == 0 && fmod(in_val2,0.1) != 0)
		{
			newoutNIT_e.replace(27,2,sub.substr(2,3));
			newoutFull_e.replace(28,2,sub.substr(2,3));
			newoutputFile_e.replace(24,2,sub.substr(2,3));

		}


		if (in_val1 < 10){
		newoutNIT_e.erase(24,1);
		newoutFull_e.erase(25,1);
		newoutputFile_e.erase(21,1);

		newoutNIT_e.insert(27,"_");
		newoutFull_e.insert(28,"_");
		newoutputFile_e.insert(24,"_");

		}
		
		NITStringsList.push_back(newNIT_e);
		FullStringsList.push_back(newFull_e);
		outNITStringsList.push_back(newoutNIT_e);
		outFullStringsList.push_back(newoutFull_e);
		OutputFiles.push_back(newoutputFile_e);
	}

	string in_str, sub, sub1;
	string newNIT, newFull, newNIT_e, newFull_e, newoutNIT, newoutFull, newoutNIT_e, newoutFull_e, newoutputFile, newoutputFile_e;

	if (reg == 0){

		Fullfile = Fullfile1;
		NITfile = NITfile1;

	}
	else{

		Fullfile = Fullfile0;
		NITfile = NITfile0;

	}

	outFullfile = outFullfile0;
	outNITfile = outNITfile0;
	outputFile = outputFile0;
}

cout << "Running GetErrorData.cc"<< endl;
 // cout << OutputFiles[1] << endl;

for ( int i = 0; i < outFullStringsList.size(); i++){
	cout << outNITStringsList[i] << endl;
	cout << outFullStringsList[i] << endl;
	cout << OutputFiles[i] << endl;
}


for (int i = 0; i < FullStringsList.size(); i++){


	string run_pre = fileloc + NITStringsList[i];
	cout << run_pre << endl;
	const char* run = run_pre.c_str();
	system(run);
	//system("/mnt/c/Users/JonnyM/Desktop/GravityWaves/Fast_Self-Forced_Inspirals/NIT_inspiral -n0 12 0.7 0.001");
	run_pre = fileloc + FullStringsList[i];
	cout << run_pre << endl;
	run = run_pre.c_str();
	system(run);
	//system("/mnt/c/Users/JonnyM/Desktop/GravityWaves/Fast_Self-Forced_Inspirals/NIT_inspiral -f0 12 0.7 0.001");
	vector<string> inputStrings;
	inputStrings.push_back(fileloc + outNITStringsList[i]);
	inputStrings.push_back(fileloc + outFullStringsList[i]);

	//inputStrings.push_back("/mnt/c/Users/JonnyM/Desktop/GravityWaves/Fast_Self-Forced_Inspirals/output/Inspiral_NIT_p12_e0.7_q0.001.dat");
	//inputStrings.push_back("/mnt/c/Users/JonnyM/Desktop/GravityWaves/Fast_Self-Forced_Inspirals/output/Inspiral_Full_p12_e0.7_q0.001.dat");
	//outputFile = "testOutputA.dat";

	GetErrors(inputStrings,OutputFiles[i]); // PROBLEM - Only saving the first file to the output folder??

	//GetErrors(inputStrings,outputFile);
	// fileloc = "/mnt/c/Users/JonnyM/Desktop/GravityWaves/Fast_Self-Forced_Inspirals/";
	cout << " " << endl;
	}

} 

// ----------------------------- Import the Inspiral Data -----------------------------------------

void GetErrors(vector<string> insp_filenames, string out_filename){
	ofstream fout;

	fout.precision(OUT_PREC);
	cout.precision(OUT_PREC);
	fout << scientific;
	cout << scientific;

	// Input: the name of eht NIT file and the name of the Full file
	vector<double> chis, ps, es, vs, ts, phis, chits, pts, ets, vts, tts, phits;
	// FIXME only load data up to t_max (rather than the entire phase space trajectory as is currently done)	
	string insp_string;
	
	double chi, p, e, v, t, phi;
	double mode = NIT_INSPIRAL;

	vector<vector<double>> NIT, Full;
	vector<Interpolant*> NIT_interp, Full_interp;

for ( int j = 0; j <= 1; j++ )
{
	vector<double> chis, ps, es, vs, ts, phis, chits, pts, ets, vts, tts, phits;
	// First pass for NIT, second pass for Full:
	//Store in the same order read in from file: chi, p, e, v, t, phi

	string insp_filename = insp_filenames[j];

	if (j == 0){
		mode = NIT_INSPIRAL;
	}else if (j == 1){
		mode = FULL_INSPIRAL;
	}
	// cout << mode << endl;
	// Check if the associated inspiral trajectory file exists
	ifstream insp(insp_filename);
	if(!insp){
		cout << "Inspiral file: " << insp_filename << " does not exist." << endl;
		exit(0);
	}
	
	// Load the NIT and Full inspiral trajectory
    if(mode == NIT_INSPIRAL){

	    cout << "Loading NIT inspiral trajectory data"  << endl;

    }else if (mode == FULL_INSPIRAL){
        
        cout << "Loading Full inspiral trajectory data"  << endl;
    }

	int test = 0;
	while(getline(insp, insp_string)){
				
		if(insp_string.at(0) == '#') continue;
		
		stringstream insp_ss(insp_string);
		
		insp_ss >> chi >> p >> e >> v >> t >> phi;
						
		if(mode == FULL_INSPIRAL)
			v = chi-v; // This makes the v used for Full actually v, meaning vError = abs(vt-v)
		
		if(test == 0) test = 1;
		else if(v <= vs.back()) break;
		else if(chi <= chis.back()) break;
		
		vs.push_back(v);
		chis.push_back(chi);
		ps.push_back(p);
		es.push_back(e);
		ts.push_back(t);
		phis.push_back(phi);


	}

vector<double> w_rs, w_phis;

for (int i = 1; i < vs.size(); i++){
	w_rs.push_back((vs[i]-vs[i-1])/(ts[i]-ts[i-1]));
	w_phis.push_back((phis[i]-phis[i-1])/(ts[i]-ts[i-1]));
}
double ccount = 0;

if(j == 0){
	NIT.push_back(chis);
	NIT.push_back(ps);
	NIT.push_back(es);
	NIT.push_back(vs);
	NIT.push_back(ts);
	NIT.push_back(phis);
	NIT.push_back(w_rs);
	NIT.push_back(w_phis);

	NIT_interp.push_back(new Interpolant(vs,chis));
	NIT_interp.push_back(new Interpolant(chis,ps));
	NIT_interp.push_back(new Interpolant(chis,es));
	NIT_interp.push_back(new Interpolant(chis,ts));
	NIT_interp.push_back(new Interpolant(chis,phis));
	NIT_interp.push_back(new Interpolant(chis,w_rs));
	NIT_interp.push_back(new Interpolant(chis,w_phis));
	NIT_interp.push_back(new Interpolant(chis,vs));

} else if (j == 1){

	Full.push_back(chis);
	Full.push_back(ps);
	Full.push_back(es);
	Full.push_back(vs);
	Full.push_back(ts);
	Full.push_back(phis);
	Full.push_back(w_rs);
	Full.push_back(w_phis);

	Full_interp.push_back(new Interpolant(vs,chis));
	Full_interp.push_back(new Interpolant(chis,ps));
	Full_interp.push_back(new Interpolant(chis,es));
	Full_interp.push_back(new Interpolant(chis,ts));
	Full_interp.push_back(new Interpolant(chis,phis));
	Full_interp.push_back(new Interpolant(chis,vs));
}

/* Remember:

chit_interp = NIT_interp[0];
pt_interp = NIT_interp[1];
et_interp = NIT_interp[2];
tt_interp = NIT_interp[3];
phit_interp = NIT_interp[4];
w_rt_interp = NIT_interp[5];
w_phit_interp = NIT_interp[6];
vt_interp = NIT_interp[7];

chi_interp = Full_interp[0];
p_interp = Full_interp[1];
e_interp = Full_interp[2];
t_interp = Full_interp[3];
phi_interp = Full_interp[4];
*/

}

chis = Full[0];
//ps = Full[1];
//es = Full[2];
vs = Full[3];		// Should only need the vs and vts
//ts = Full[4];
//phis = Full[5];

chits = NIT[0];
//pts = NIT[1];
//ets = NIT[2];
vts = NIT[3];
//tts = NIT[4];
//phits = NIT[5];

// omega*t,phi and v error functions:

// obtain the dense vs to use to plug into all the interpolants:
int i_max;
double M_solar, Deltat_sec;
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

double t_max = i_max*Deltat;
double chit, pt, et, vt, tt, phit;
vector<double> t_dense, tt_dense, v_dense, vt_dense, chit_dense, chi_vec;
	
	//v = vs[0];
	if (chits.back() > chis.back()){
		chi_vec = chis;
	} else if (chits.back() < chis.back()){
		chi_vec = chits;
	}

	chit = chi_vec[0];

	//cout << "hi" << endl;
	while(chit < chi_vec[chi_vec.size() -1]){
		v = Full_interp[5]->eval(chit);
		vt = NIT_interp[7]->eval(chit);
		p = Full_interp[1]->eval(chit);
		e = Full_interp[2]->eval(chit);

		pt = NIT_interp[1]->eval(chit);
		et = NIT_interp[2]->eval(chit);
			
		v_dense.push_back(v);
		vt_dense.push_back(vt);
		chit_dense.push_back(chit);

		t = Full_interp[3]->eval(chit);
		tt = NIT_interp[3]->eval(chit);
		
		t_dense.push_back(t);
		tt_dense.push_back(tt);
		
		if(t > t_max) break;	
		// v += 2.0*M_PI/10.;
		chit += 2.0*M_PI/10.;// NIT_interp[0]->eval(v);	
	}

	 // cout << "hey" << endl;

/*
for ( int i = 0; i < v_dense.size();i++){
	cout << v_dense[i] - vt_dense[i] << endl; // vt is essentially the same as v 
}
*/

	int i_max_test = floor(t_dense.back()/Deltat)-1;
	if(i_max > i_max_test) i_max = i_max_test;

	vector<double> wtR_Error, wtPhi_Error, phi_Error, v_Error;
	double dwtR,dwtPhi, dphi, dv;
	double w_r, w_phi;
	double DeltaChi = 2.0*M_PI/10.;

	for(int i = 0; i <= chit_dense.size()-1; i ++)
	{
	//v = v_dense[i];
	//vt = vt_dense[i];
	chit = i*DeltaChi;

	//chi = Full_interp[0]->eval(v);
	v = Full_interp[5]->eval(chit);
	p = Full_interp[1]->eval(chit);
	e = Full_interp[2]->eval(chit);
	t = Full_interp[3]->eval(chit);
	phi = Full_interp[4]->eval(chit);

	 // chit = NIT_interp[0]->eval(vt);
	vt = NIT_interp[7]->eval(chit);
	pt = NIT_interp[1]->eval(chit);
	et = NIT_interp[2]->eval(chit);
	tt = NIT_interp[3]->eval(chit);
	phit = NIT_interp[4]->eval(chit);
	w_r = NIT_interp[5]->eval(chit);
	w_phi = NIT_interp[6]->eval(chit);

// Big for loop for calculating all the error values
	dwtR = w_r*abs(((1/M)*((tt-Z0t(pt,et,vt)) - t)));
	dwtPhi = w_phi*abs(((1/M)*((tt-Z0t(pt,et,vt)) - t)));
	dphi = abs((phit-Z0phi(pt,et,vt)) - phi);
	dv = abs(vt-v);

	wtR_Error.push_back(dwtR);
	wtPhi_Error.push_back(dwtPhi);
	phi_Error.push_back(dphi);
	v_Error.push_back(dv);

	}// end of big for loop for calculating all the error values


// Output error data to out_filename:
fout.open(out_filename);
cout << "Outputting error data to " << out_filename << endl;
	double t_sec;
	for(int i = 0; i <= chit_dense.size() -1; i ++){
		
		fout << chit_dense[i] << " " << wtR_Error[i] << " " << wtPhi_Error[i] << " " << phi_Error[i] << " " << v_Error[i] << endl;

	}

} // End of GetErrors

// ------------------------------- Error Functions ----------------------------------- //

//Tr: --- For this one we need elliptical functions: EllipticK, EllipticPi,EllipticF
double Tr(double p, double e, double v)
{
	return (2*Sqrt((-4*Power(e,2) + Power(-2 + p,2))/(-6 + 2*e + p))*
     ((-1 - e)*(2 + 2*e - p)*(-4 + p)*p*(-6 + 2*e + p)*EllipticE((4*e)/(-6 + 2*e + p)) + 
       (1 + e)*(2 + 2*e - p)*p*(36 + 2*Power(e,2)*(-2 + p) + (-14 + p)*p)*EllipticK((4*e)/(-6 + 2*e + p)) + 
       2*(6 + 2*e - p)*((-2 - 2*e + p)*(8 + p - Power(p,2) + Power(e,2)*(-8 + 3*p))*
           EllipticPi((2*e*(-4 + p))/((1 + e)*(-6 + 2*e + p)),(4*e)/(-6 + 2*e + p)) - 
          (-1 + e)*Power(1 + e,2)*Power(-4 + p,2)*EllipticPi((16*e)/(12 + 8*e - 4*Power(e,2) - 8*p + Power(p,2)),
            (4*e)/(-6 + 2*e + p)))))/((-1 + e)*Power(1 + e,2)*(2 + 2*e - p)*Power(-4 + p,2));
}

//Phir:
double Phir(double p, double e)
{
	return 4*Sqrt(p/(-6 + 2*e + p))*EllipticK((4*e)/(-6 + 2*e + p));
}

//Z0t: takes tilded values as arguments
double Z0t(double p, double e, double v)
{
return -((p*Sqrt((-4*Power(e,2) + Power(-2 + p,2))*(-6 + 2*e + p))*EllipticEIncomp((Pi - v)/2.,(4*e)/(-6 + 2*e + p)))/((-1 + Power(e,2))*(-4 + p))) +
   (p*Sqrt((-4*Power(e,2) + Power(-2 + p,2))/(-6 + 2*e + p))*EllipticF((Pi - v)/2.,(4*e)/(-6 + 2*e + p)))/(-1 + Power(e,2)) +
   ((-Pi + v)*((-2*p*Sqrt((-4*Power(e,2) + Power(-2 + p,2))*(-6 + 2*e + p))*EllipticE((4*e)/(-6 + 2*e + p)))/((-1 + Power(e,2))*(-4 + p)) +
        (2*p*Sqrt((-4*Power(e,2) + Power(-2 + p,2))/(-6 + 2*e + p))*EllipticK((4*e)/(-6 + 2*e + p)))/(-1 + Power(e,2)) -
        (4*Sqrt((-4*Power(e,2) + Power(-2 + p,2))/(-6 + 2*e + p))*(8 + p - Power(p,2) + Power(e,2)*(-8 + 3*p))*EllipticPi((2*e)/(-1 + e),(4*e)/(-6 + 2*e + p)))/
         (Power(-1 + e,2)*(1 + e)*(-4 + p)) + (16*Sqrt((-4*Power(e,2) + Power(-2 + p,2))/(-6 + 2*e + p))*EllipticPi((4*e)/(-2 + 2*e + p),(4*e)/(-6 + 2*e + p)))/
         (-2 + 2*e + p)))/(2.*Pi) + (2*Sqrt((-4*Power(e,2) + Power(-2 + p,2))/(-6 + 2*e + p))*(8 + p - Power(p,2) + Power(e,2)*(-8 + 3*p))*
      EllipticPiIncomp((2*e)/(-1 + e),(-Pi + v)/2.,(4*e)/(-6 + 2*e + p)))/(Power(-1 + e,2)*(1 + e)*(-4 + p)) -
   (8*Sqrt((-4*Power(e,2) + Power(-2 + p,2))/(-6 + 2*e + p))*EllipticPiIncomp((4*e)/(-2 + 2*e + p),(-Pi + v)/2.,(4*e)/(-6 + 2*e + p)))/(-2 + 2*e + p) -
   (e*p*Sqrt((-4*Power(e,2) + Power(-2 + p,2))*(-6 + p - 2*e*Cos(v)))*Sin(v))/((-1 + Power(e,2))*(-4 + p)*(1 + e*Cos(v)));
}

//Z0phi: takes tildeed values as arguments
double Z0phi(double p, double e, double v)
{
	return (2*Sqrt(p/(-6 + 2*e + p))*(Pi*EllipticF((Pi - v)/2.,(4*e)/(-6 + 2*e + p)) + (-Pi + v)*EllipticK((4*e)/(-6 + 2*e + p))))/Pi;
}



// Elliptice Integrals and functions:

double EllipticK(double k){
	return gsl_sf_ellint_Kcomp(sqrt(k), GSL_PREC_DOUBLE);
}

double EllipticF(double phi, double k){
	return gsl_sf_ellint_F(phi, sqrt(k), GSL_PREC_DOUBLE) ;
}

double EllipticE(double k){
	return gsl_sf_ellint_Ecomp(sqrt(k), GSL_PREC_DOUBLE);
}

double EllipticEIncomp(double phi, double k){
	return gsl_sf_ellint_E(phi, sqrt(k), GSL_PREC_DOUBLE) ;
}

double EllipticPi(double n, double k){
	return gsl_sf_ellint_Pcomp(sqrt(k), -n, GSL_PREC_DOUBLE);
}

double EllipticPiIncomp(double n, double phi, double k){
	return gsl_sf_ellint_P(phi, sqrt(k), -n, GSL_PREC_DOUBLE);
}

