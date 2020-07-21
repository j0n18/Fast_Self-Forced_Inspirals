// ---------------------- Get Error Data from NIT and Full inspirals ------------------------------

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

//-------------------------------- Instantiate functions ----------------------------------------- //
vector<string> InputStrings(int which, int reg); //, int pmin, int pmax, int emin, int emax
void GetErrors(vector<string> NIT_output_filenames, vector<string> Full_output_filenames, vector<string> Error_out_filenames);
vector<double> ImportData(string insp_filename,	double mode, int which_param);
double Z0t(double p, double e, double v);
double Z0phi(double p, double e, double v);
double EllipticK(double k);
double EllipticF(double phi, double k);
double EllipticE(double k);
double EllipticEIncomp(double phi, double k);
double EllipticPi(double n, double k);
double EllipticPiIncomp(double n, double phi, double k);


int main(int argc, char* argv[])
{
 cout << "Running GetErrorData.cc" << endl;
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

// ----------------------------- Generate the Inspiral Data --------------------------------------- //

// Acquire command strings and output file strings:
cout << "Building command strings for NIT and Full inspirals" << endl;

string fileloc = "/mnt/c/Users/JonnyM/Desktop/GravityWaves/Fast_Self-Forced_Inspirals/";
int reg = 1;


vector<string> NITcommands = InputStrings(0,reg);
vector<string> Fullcommands = InputStrings(1,reg);
vector<string> NIToutputFiles = InputStrings(2,reg);
vector<string> FulloutputFiles = InputStrings(3,reg);
vector<string> OutputFiles = InputStrings(4,reg);

// Check: 
//cout << NITcommands[0]<< endl;
//cout << Fullcommands[0]<< endl;
//cout << NIToutputFiles[0]<< endl;
//cout << FulloutputFiles[0]<< endl;
//cout << OutputFiles[0]<< endl;

//cout << "testing ImportData" << endl;
//string yoyo = fileloc + FulloutputFiles[0];
//vector<double> test_vec = ImportData(yoyo,FULL_INSPIRAL,0); 
//cout << test_vec[0] << endl;


for ( int i = 0; i < OutputFiles.size();i++){
	cout << NIToutputFiles[i] << endl;
}
cout << " " << endl;

// Execute commands:
cout << "Executing commands for NIT anf Full inspirals" << endl;
for (int i = 0; i < Fullcommands.size(); i++){

	string run_pre = fileloc + NITcommands[i];
	cout << run_pre << endl;
	const char* run = run_pre.c_str();
	system(run);

	run_pre = fileloc + Fullcommands[i];
	cout << run_pre << endl;
	run = run_pre.c_str();
	system(run);

}

cout << "Calling GetErrors() ..." << endl;
cout << " " << endl;
	GetErrors(NIToutputFiles, FulloutputFiles, OutputFiles); 

}


// ---------------------------------------- Functions -------------------------------------------------- //

//vector<string> ImportStrings(string NITfile0, string Fullfile0, string outNITfile0, string outFullfile0, string outputFile0, int reg, double dp, double de, int pmax, int pmin, double emax, double emin){
vector<string> InputStrings(int which, int reg){
string fileloc, NITfile0, Fullfile0, NITfile, Fullfile, outNITfile0, outFullfile0, outNITfile, outFullfile,NITfile1, Fullfile1;
vector<string> NITStringsList, FullStringsList, outNITStringsList, outFullStringsList;

fileloc = "/mnt/c/Users/JonnyM/Desktop/GravityWaves/Fast_Self-Forced_Inspirals/";

 // Initial strings
NITfile0 = "NIT_inspiral -n0 10 0.1 0.001";
Fullfile0 = "NIT_inspiral -f0 10 0.1 0.001"; // The first part should be NIT_inspiral for directory reasons. -f0 signifies Full.

NITfile1 = "NIT_inspiral -n 10 0.1 0.001";
Fullfile1 = "NIT_inspiral -f 10 0.1 0.001";

double dp = 1; 
double de = 0.1;
int pmax = 12; // make 14 // p > 11 does not work for -n and (reliably) for -f
int pmin = 10; // make 8
double emax = 0.75; // make 0.75; This will be rounded to 0.7
double emin = 0.1; // make 0.1
double floatValue, in_val, in_val1,in_val2;

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

vector<string> out;

if (which == 0){

out = NITStringsList;

} else if (which == 1){

out = FullStringsList;

}else if (which == 2){

out = outNITStringsList;

}else if (which == 3){

out = outFullStringsList;

}else if (which == 4){

out = OutputFiles;

}

return out;

}

//GetErrors:

void GetErrors(vector<string> NIT_output_filenames, vector<string> Full_output_filenames, vector<string> Error_out_filenames){
string fileloc = "/mnt/c/Users/JonnyM/Desktop/GravityWaves/Fast_Self-Forced_Inspirals/";
vector<double> chis, ps, es, vs, ts, phis, chits, pts, ets, vts, tts, phits, w_rs, w_phis,w_rts,w_phits;
vector<Interpolant*> NIT_interp, Full_interp;
string inNIT, inFull, Error_out_filename;
ofstream fout; // this inside the for loop or here?

fout.precision(OUT_PREC);
cout.precision(OUT_PREC);
fout << scientific;
cout << scientific;

for ( int i = 0; i < NIT_output_filenames.size(); i++){
	
	inFull = fileloc + Full_output_filenames[i];
	inNIT = fileloc + NIT_output_filenames[i];

	chis = ImportData(inFull,FULL_INSPIRAL,0);
	ps = ImportData(inFull,FULL_INSPIRAL,1);
	es = ImportData(inFull,FULL_INSPIRAL,2);
	vs = ImportData(inFull,FULL_INSPIRAL,3);
	ts = ImportData(inFull,FULL_INSPIRAL,4);
	phis = ImportData(inFull,FULL_INSPIRAL,5);
	w_rs = ImportData(inFull,FULL_INSPIRAL,6);
	w_phis = ImportData(inFull,FULL_INSPIRAL,7);

	chits = ImportData(inNIT,NIT_INSPIRAL,0);
	pts = ImportData(inNIT,NIT_INSPIRAL,1);
	ets = ImportData(inNIT,NIT_INSPIRAL,2);
	vts = ImportData(inNIT,NIT_INSPIRAL,3);
	tts = ImportData(inNIT,NIT_INSPIRAL,4);
	phits = ImportData(inNIT,NIT_INSPIRAL,5);
	w_rts = ImportData(inNIT,NIT_INSPIRAL,6);
	w_phits = ImportData(inNIT,NIT_INSPIRAL,7);

	NIT_interp.push_back(new Interpolant(vts,chits));
	NIT_interp.push_back(new Interpolant(chits,pts));
	NIT_interp.push_back(new Interpolant(chits,ets));
	NIT_interp.push_back(new Interpolant(chits,tts));
	NIT_interp.push_back(new Interpolant(chits,phits));
	NIT_interp.push_back(new Interpolant(chits,w_rts));
	NIT_interp.push_back(new Interpolant(chits,w_phits));
	NIT_interp.push_back(new Interpolant(chits,vts));

	Full_interp.push_back(new Interpolant(vs,chis));
	Full_interp.push_back(new Interpolant(chis,ps));
	Full_interp.push_back(new Interpolant(chis,es));
	Full_interp.push_back(new Interpolant(chis,ts));
	Full_interp.push_back(new Interpolant(chis,phis));
	Full_interp.push_back(new Interpolant(chis,vs));

	double DeltaChi = 2.0*M_PI/10.;
	double chi_max;
	if (chits.back() < chis.back()){
		chi_max = chits.back();
	} else if (chits.back() > chis.back()){
		chi_max = chis.back();
	}

	double j_max = chi_max/DeltaChi;
	double chit,pt,et,vt,tt,phit;
	double chi,p,e,v,t,phi,w_r,w_phi;
	double dwtR,dwtPhi,dphi,dv;

	// for now:
	double M  = 1; // This likely must change 
	Error_out_filename = Error_out_filenames[i];

	fout.open(Error_out_filename); // fout opens here, is populated in the sub loop and closed at the beginning of the next iteration

	for(int j = 0; j < j_max; j++)
	{

	chit = j*DeltaChi;

	v = Full_interp[5]->eval(chit);
	p = Full_interp[1]->eval(chit);
	e = Full_interp[2]->eval(chit);
	t = Full_interp[3]->eval(chit);
	phi = Full_interp[4]->eval(chit);

	vt = NIT_interp[7]->eval(chit);
	pt = NIT_interp[1]->eval(chit);
	et = NIT_interp[2]->eval(chit);
	tt = NIT_interp[3]->eval(chit);
	phit = NIT_interp[4]->eval(chit);
	w_r = NIT_interp[5]->eval(chit);
	w_phi = NIT_interp[6]->eval(chit);

	dwtR = w_r*abs(((1/M)*((tt-Z0t(pt,et,vt)) - t)));
	dwtPhi = w_phi*abs(((1/M)*((tt-Z0t(pt,et,vt)) - t)));
	dphi = abs((phit-Z0phi(pt,et,vt)) - phi);
	dv = abs(vt-v);

	fout << chit << " " << dwtR << " " << dwtPhi << " " << dphi << " " << dv << endl;

	}// end of sub loop

} // end of outer loop

} // end of GetErrorData

// Import data:
vector<double> ImportData(string insp_filename,	double mode, int which_param){
	//ofstream fout; \\ this doesnt seem like it should be here
	
	string insp_string;
	
	double chi, p, e, v, t, phi;

	vector<double> chis, ps, es, vs, ts, phis;

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

		chis.push_back(chi);
		ps.push_back(p);
		es.push_back(e);
		vs.push_back(v);
		ts.push_back(t);
		phis.push_back(phi);

	}

vector<double> w_rs, w_phis;

for (int i = 1; i < vs.size(); i++){

	if (which_param == 6)
	{
		w_rs.push_back((vs[i]-vs[i-1])/(ts[i]-ts[i-1]));
	} 
	else if (which_param == 7)
	{
		w_phis.push_back((phis[i]-phis[i-1])/(ts[i]-ts[i-1]));
	}
}

vector<double> out;

if (which_param == 0){

	out = chis;

}else if (which_param == 1){

		out = ps;

}else if (which_param == 2){

		out = es;

}else if (which_param == 3){

		out = vs;

}else if (which_param == 4){

		out = ts;
}
else if (which_param == 5){
		out = phis;
}
else if (which_param == 6){
		out = w_rs;
}
else if (which_param == 7){
		out = w_phis;
}

return out;

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