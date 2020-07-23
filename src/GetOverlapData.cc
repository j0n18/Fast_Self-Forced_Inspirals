// ---------------------- Get Overlap Data from NIT and Full inspirals ------------------------------

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


// -------------------------------- Include Statements -------------------------------------------- //
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
#define Tanh(x)			(tanh((double) (x)))

using namespace libconfig;
Config cfg;

// Instantiate functions:
void GetOverlap(int Comp, vector<string> insp_filenames, string out_filename, string NITwaveform_filename, string Fullwaveform_filename);
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
vector<complex<double>> waveform_DFT(vector<complex<double>> Amps);
double Eqn32(vector<complex<double>> x,vector<complex<double>> h, vector<double> f);
double Eqn34(double xx, double hh, double xh);
complex<double> Sh(double f1);
vector<complex<double>> DataLoad(string filename);
vector<complex<double>> waveform_inverse_DFT(vector<complex<double>> Amps);
double RunOverlap(vector<complex<double>> FullAmps, vector<complex<double>> NITAmps, vector<double> fs);
vector<string> StringBuilder(int Comp, int returntype);
void fft(int N, fftw_complex *in, fftw_complex *out);
void ifft(int N, fftw_complex *in, fftw_complex *out);
vector<double> TimesLoad(string filename);
int maximum_index(vector<double> ts);
double complex_modulus(complex<double> x);
complex<double> Sc(double f);
double Poms(double f);
double Pacc(double f);
// ----------------------------- Generate the Inspiral Data --------------------------------------- //
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

string fileloc = "/mnt/c/Users/JonnyM/Desktop/GravityWaves/Fast_Self-Forced_Inspirals/";

int Comp;

int wntn = 0; // comparing -w -n with -t -n
int wftf = 1; // comparing -w -f with -t -f
int tntf = 0; // comparing -t -n with -t -f
int wnwf = 0; // comparing -w -n with -w -f

if(wntn == 1){

	Comp = 0;

} else if (wftf == 1){

	Comp = 1;

} else if (tntf == 1){

	Comp = 2;

} else if (wnwf == 1){

	Comp = 3;

}

vector<string> NITStringsList = StringBuilder(Comp,0);
vector<string> FullStringsList = StringBuilder(Comp,1);
vector<string> outNITStringsList = StringBuilder(Comp,2);
vector<string> outFullStringsList = StringBuilder(Comp,3);
vector<string> OutputFiles = StringBuilder(Comp,4);
vector<string> WaveformCommands = StringBuilder(Comp,5);
vector<string> WaveformOutputFilesNIT = StringBuilder(Comp,6);
vector<string> WaveformOutputFilesFull = StringBuilder(Comp,7);

string run_pre;
const char* run;
cout << "Running GetOverlapData.cc"<< endl;
for (int i = 0; i < NITStringsList.size(); i++){
	if (wftf != 1){
	run_pre = fileloc + NITStringsList[i];
	cout << run_pre << endl;
	run = run_pre.c_str();
	system(run);
	}
	//system("/mnt/c/Users/JonnyM/Desktop/GravityWaves/Fast_Self-Forced_Inspirals/NIT_inspiral -n0 12 0.7 0.001");
	if (wntn != 1){
	run_pre = fileloc + FullStringsList[i];
	cout << run_pre << endl;
	run = run_pre.c_str();
	system(run);
	}
	//system("/mnt/c/Users/JonnyM/Desktop/GravityWaves/Fast_Self-Forced_Inspirals/NIT_inspiral -f0 12 0.7 0.001");
	//---------------------------------------------------------------------------------------------------------- //
	//run the waveform commands here:
	run_pre = fileloc + WaveformCommands[2*i];
	cout << run_pre << endl;
	run = run_pre.c_str();
	system(run);

	//run the waveform commands here:
	run_pre = fileloc + WaveformCommands[2*i + 1];
	cout << run_pre << endl;
	run = run_pre.c_str();
	system(run);
	
	vector<string> inputStrings;
	if (wftf != 1){
	inputStrings.push_back(fileloc + outNITStringsList[i]);
	}

	if(wntn != 1){
	inputStrings.push_back(fileloc + outFullStringsList[i]);
	}
	
	//cout << inputStrings[0] << endl;
	
	GetOverlap(Comp,inputStrings,OutputFiles[i],WaveformOutputFilesNIT[i],WaveformOutputFilesFull[i]); // WaveformCommands is the command string for -t or -w

	cout << " " << endl;
	}

} 

// ----------------------------- Import the Inspiral Data -----------------------------------------

void GetOverlap(int Comp, vector<string> insp_filenames, string out_filename, string NITwaveform_filename, string Fullwaveform_filename){
	ofstream fout;

	fout.precision(OUT_PREC);
	cout.precision(OUT_PREC);
	fout << scientific;
	cout << scientific;

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
	
// -------------------------------------------- Waveform Fourier Analysis ------------------------------------ //
// 1. Calculate waveforms for all different waveform strings:


// 2. load in the waveform output (ampitude) data:

vector<complex<double>> FullAmps, NITAmps, FullAmpsIn, NITAmpsIn;
vector<double> ts;

cout << NITwaveform_filename << endl;
cout << Fullwaveform_filename << endl;

FullAmps = DataLoad(Fullwaveform_filename);
NITAmps = DataLoad(NITwaveform_filename);

ts = TimesLoad(Fullwaveform_filename);

// make sure xand h are the same length and of an odd length:
	int n;

	if ( FullAmps.size() < NITAmps.size())
	{
		n = FullAmps.size();
	} 

	else
	{
		n = NITAmps.size();
	}
	

	if (n % 2 == 0){
		n -= 1;
	}

	for ( int i = 1; i <= n; i++){ // including N gets the last 0 just like mathematica
		FullAmpsIn.push_back(FullAmps[i]);
		NITAmpsIn.push_back(NITAmps[i]);
	}



// Staying within the double for loop:
cout << "Computing the Fourier Amplitudes" << endl;

// 3. find the Fourier Amplitudes:
	vector<complex<double>> x, h, xP, hP, xN, hN;
	vector<double> fP, fN;

	x = waveform_DFT(NITAmpsIn); // NIT is x 

	h = waveform_DFT(FullAmpsIn); // Full is h

	cout << n << endl;

	double tsmaxind = maximum_index(ts); // max_element does not work
	double df = 1/ts[tsmaxind];//(fmax-fmin)/n; //This needs to come from the times of the TimesLoad

	vector<double> fs;

	int iter = (n-1)/2;

	for ( int i = 0; i <= iter; i++)
	{
		fs.push_back(i*df); // starting at 1 to remove the zero frequency value
	}

	for ( int i = -iter; i <= -1; i++)
	{
		fs.push_back(i*df); // starting at 1 to remove the zero frequency value
	}
	
	fs[0] = 1e-8; //matches mathematica
	
	cout << "beginning of fs" << endl;
	for ( int i = 0; i < 20; i++){
		cout << fs[i] << endl;
	}
	cout <<"..." << endl;
	cout << "end of fs" << endl;
	for ( int i = fs.size()-20; i < fs.size(); i++){
		cout << fs[i] << endl;
	}
	
	
	// the negative frequencies seem to be correct...

	for ( int i = 0; i < iter; i++) // i = 0 should correspond in fs to a frequency of f = df since we removed the zero frequency element
	{
		xP.push_back(x[i]);
		hP.push_back(h[i]);
		fP.push_back(fs[i]);
	}

	for ( int i = iter; i < x.size(); i++){ // x.size() should correspond to n-1 since we erased the zero frequency element
		xN.push_back(x[i]);
		hN.push_back(h[i]);
		fN.push_back(fs[i]);
	}

cout << fN.size() << endl;
cout << fP.size() << endl;
// 4. Calculate the Overlap using Eqn 32 and 34:

	double Exx, Ehh, Exh,ExxP, EhhP, ExhP, ExxN, EhhN, ExhN;

	Exx = Eqn32(x,x,fs);
	Ehh = Eqn32(h,h,fs);
	Exh = Eqn32(x,h,fs);

	// Splitting into Positive and Negative seems to be giving more reliable answers... 
	ExxP = Eqn32(xP,xP,fP);
	EhhP = Eqn32(hP,hP,fP);
	ExhP = Eqn32(xP,hP,fP);

	ExxN = Eqn32(xN,xN,fN);
	EhhN = Eqn32(hN,hN,fN);
	ExhN = Eqn32(xN,hN,fN);

	
	
	/*
	cout << "HELLO" << endl;
	cout << Exx << endl;
	cout << Ehh << endl;
	cout << Exh << endl;
	cout << "HELLO" << endl;
	*/

	double similarity, similarity1, similarity2;

	similarity = Eqn34(Exx, Ehh, Exh);
	similarity1 = Eqn34(ExxP + ExxN, EhhP + EhhN, ExhP + ExhN); // values are not consistent with more e's
	similarity2 = RunOverlap(FullAmpsIn,NITAmpsIn,fs); // This line causes the values of similarity 1 and 2 to change??
														// similarity values not consistent 

 // end of loop through waveform filenames
	string fileloc = "/mnt/c/Users/JonnyM/Desktop/GravityWaves/Fast_Self-Forced_Inspirals/";
	// make a general file to output overlap scores to:
	
	 // string fileloc = "/mnt/c/Users/JonnyM/Desktop/GravityWaves/Fast_Self-Forced_Inspirals/";
	 string name1, name2;

	name1 = insp_filenames[0].substr(fileloc.size(),insp_filenames[0].size());
	
	if (Comp == 2 || Comp == 3){
	name2 = insp_filenames[1].substr(fileloc.size(),insp_filenames[1].size());
	}
	
	fout.open(out_filename);


	cout << "Outputting error data to " << out_filename << endl;

	cout << "Overlap Score using Exx... for " << name1 << " and " << name2 << "\n" << similarity << "\n" << endl;
	cout << "Overlap Score using ExxP and ExxN ... for " << name1 << " and " << name2 << "\n" << similarity1 << "\n" << endl;
	cout << "Overlap Score using SimulationTools Method for " << name1 << " and " << name2 << "\n" << similarity2 << "\n" << endl;
	fout << "Overlap Score for " << name1 << " and " << name2 << "\n" << similarity << "\n" << endl;
} // End of GetOverlap

// ------------------------------- Functions ----------------------------------- //

vector<complex<double>> waveform_DFT(vector<complex<double>> Amps){ // SOMETHING MAY BE HERE

	int N = Amps.size();
  //  fftw_complex *in, *out;
  //  fftw_plan plan;

  fftw_complex in[N];
  fftw_complex out[N];

/*
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
*/
	// Build input:
	//The input to this should just be a complex number whose real part is N_re and whose imaginary part is N_im
	for ( int i = 0; i <Amps.size(); i++)
	{
		in[i][0] = Amps[i].real();
		in[i][1] = Amps[i].imag();
	}

	//fftw_execute(plan); // this populates the out vector with the fourier amplitudes

	fft(N,in,out);

	vector<complex<double>> Out1, Out;

	for ( int i = 0; i < Amps.size(); i++){
		Out.push_back(Complex (out[i][0]/Sqrt(N),out[i][1]/Sqrt(N)));// the /Sqrt(N) makes this yield the same output as mathematica
	}


	for ( int i = 0; i < Amps.size(); i++){
		if( i == 0){
			Out1.push_back(Out[i]);
		}else{
			Out1.push_back(Out[Out.size() - i]);
		}
	}

	return Out1; // This is reordered to return in the same format as mathematica. Gives same overlap as before.

} // end of waveform_DFT


// -------------------------------------------------------------------------------------

vector<complex<double>> waveform_inverse_DFT(vector<complex<double>> Amps){

	int N = Amps.size();
   // fftw_complex *in, *out;
  //  fftw_plan plan;
	fftw_complex in[N];
	fftw_complex out[N];
   // in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
   // out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
   // plan = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

	// Build input:
	//The input to this should just be a complex number whose real part is N_re and whose imaginary part is N_im
	for ( int i = 0; i <Amps.size(); i++)
	{
		in[i][0] = Amps[i].real();
		in[i][1] = Amps[i].imag();
	}

	//fftw_execute(plan); // this populates the out vector with the inverse fourier amplitudes

	ifft(N,in,out);

	vector<complex<double>> Out1, Out;

	for ( int i = 0; i < Amps.size(); i++){

		Out.push_back(Complex (out[i][0],out[i][1])); //Not sure about this /N
	}

		for ( int i = 0; i < Amps.size(); i++){
		if( i == 0){
			Out1.push_back(Out[i]);
		}else{
			Out1.push_back(Out[Out.size() - i]);
		}
	}

	return Out1;

} // end of waveform_inverse_DFT


void fft(int N, fftw_complex *in, fftw_complex *out){

fftw_plan plan  = fftw_plan_dft_1d(N,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
fftw_execute(plan);
fftw_destroy_plan(plan);
fftw_cleanup();

}

void ifft(int N, fftw_complex *in, fftw_complex *out){

fftw_plan plan  = fftw_plan_dft_1d(N,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);
fftw_execute(plan);
fftw_destroy_plan(plan);
fftw_cleanup();
/*
for (int i = 0; i < N; i++){
		out[i][0] /= N; // properly scales output
		out[i][1] /= N;
	}
*/
}


// ---------------------------------------------------------------------------------------

double RunOverlap(vector<complex<double>> FullAmps, vector<complex<double>> NITAmps, vector<double> fs){

vector<complex<double>> x, h;

 // 1. Fourier[]
	x = waveform_DFT(NITAmps); // NIT is x ; "NITAmps" are values of the waveform itself; this one has agreement between mathematica and the DFT c++ functions
	h = waveform_DFT(FullAmps); // Full is h

	 //cout << x.size() << endl;
	 cout << FullAmps.size() << endl;
	 //cout << x[0] << endl;
	 //cout << fs.size() << endl;

// 2. Compute Normalizations:
	vector<complex<double>> xIn, hIn, invx_pre, invh_pre;
	vector<double> invx, invh;
	complex<double> inputx, inputh; // double or complex<double> ?

	for ( int i = 0; i < x.size(); i++){

		inputx = (Power(abs(x[i]),2))/Sh(fs[i]);
		inputh = (Power(abs(h[i]),2))/Sh(fs[i]); 

		//somewhere in this process, a complex part is lost that is present in mathematica.				
		// This is likely accounting for some or all of the discrepancy in the values when applying the inverse transformation
		// the real part is always shared and correct but in this code the complex part always remains zero. 

		xIn.push_back(inputx);
		hIn.push_back(inputh);
	}
	
	cout << "beginning of h abs^2/sn" << endl;
	for ( int i = 0; i < 20; i++){
		cout << hIn[i] << endl;
	}
	cout << " " << endl;
	cout << "end of h abs^2/sn" << endl;
	for ( int i = hIn.size()-20; i < hIn.size(); i++){
		cout << hIn[i] << endl;
	}
	cout << " " << endl;
	cout << "Sc[-1.265110159467e-6]:" << endl;
	cout << Sc(-1.265110159467e-6) << endl;
	
	invx_pre = waveform_inverse_DFT(xIn); // off from what mathematica returns but x and h are relatively close in the transform.
	invh_pre = waveform_inverse_DFT(hIn);

	for ( int i = 0; i < invx_pre.size(); i++){
		invx.push_back(abs(invx_pre[i]));
		invh.push_back(abs(invh_pre[i]));
	}	

	int norm1_maxind = maximum_index(invx);
	int norm2_maxind = maximum_index(invh);

	double norm1 = 4 * invx[norm1_maxind]; // using Max instead of First
	double norm2 = 4 * invh[norm2_maxind];

	cout << "norm magnitudes" << endl;
	cout << norm1 << endl;
	cout << norm2 << endl;
	cout << "end of norm magnitudes"<< endl;
// 3. Compute ther inner product integrand:

	vector<complex<double>> integrand;
	complex<double> integrand_pre;
	for ( int i = 0; i < x.size(); i++){
		integrand_pre = (conj(x[i]) * h[i])/Sh(fs[i]); // The integrand values are correct
		integrand.push_back(integrand_pre);
	}

// 4. Compute the match:

 vector<complex<double>> ret;
 vector<double> retout;

 ret = waveform_inverse_DFT(integrand); // The inverse transform returns values that are slightly off but are ~ correct in magnitude

	for ( int i = 0; i < ret.size(); i++){
		retout.push_back(abs(ret[i])); 
	}	

	double retmaxind = maximum_index(retout); // This is the index of the highest value in retout
	return 4 * (retout[retmaxind])/Sqrt((norm1 * norm2));

}


// ---------------------------------------------------------------------------------------

double Eqn32(vector<complex<double>> x,vector<complex<double>> h, vector<double> f){
	double sum = 0;

	int K = (1/2) * (x.size() - 1);  // by the time this function is called x should be of length n
	for ( int i = 0; i <= K; i++){
		sum += abs(( x[i]*conj(h[i]) + conj(x[i])*h[i] )/(Sh(f[i]))); // abs() of a complex number gives the modulus which we want
	}
	
	return 2*sum;
}

double Eqn34(double xx, double hh, double xh){
	return xh/(Sqrt(abs(xx) * abs(hh))); // xx, xh and hh shold be the results of plugging into Eqn32 for (x,x) (x,h) and (h,h)
}

//Functions within Sh:

double Poms(double f){
	return 2.25e-22*(1 + 1.6e-11/Power(f,4));
}

double Pacc(double f){
	return (9*(1 + 1.6e-7/Power(f,2))*(1 + 2.44140625e8*Power(f,4)))/1.e30;
}

complex<double> Sc(double f){
	f = abs(f);
	return (9*Power(exp(1),-Power(f,0.165) + 299*f*Sin(611*f))*(1 + Tanh(1340*(0.00173 - f))))/(1.e45*Power(f,2.3333333333333335));
	
}


complex<double> Sh(double f){ // works but is dropping any imaginary part
	f = abs(f);
	double fStar = 0.01909;
	double L = 2.5e9;
	double A = 9e-45;
	double alpha = 0.165;
	int beta = 299;
	int kappa = 611;
	int gamma = 1340;
	double fk = 0.00173;
	return 10/(3*L*L) * (Poms(f) + 4 * Pacc(f)/Power(2*Pi*f,4)) * (1 + (1/10)*Power(f/fStar,2)) + Sc(f);
}


// basic maximum function:
int maximum_index(vector<double> ts){
	int maxind = 0;
	for (int i = 0; i < ts.size(); i++){

		if(ts[maxind] < ts[i]){

			maxind = i;

		}else if (ts[maxind] >= ts[i]){

			continue;
		}

	}
	return maxind;
}//end of maximum_index

double complex_modulus(complex<double> x){
	double a = x.real() * x.real();
	double b = x.imag() * x.imag();
	return Sqrt(a + b);
}

vector<complex<double>> DataLoad(string filename){

	ifstream inFile;
	string line;
	double t_col, Re, Im;
	int iter, input, num_cols;
	vector<double> t_vals;

	inFile.open(filename); 
	if (!inFile)
		{
		cout << "Unable to open Waveform data file \n";
		exit(1); // terminate with error
		}

	vector<complex<double>> HTeukAmps; //HTeuk values 
	int k = 0;
	// loop through each orbit (rows in data file)
	while(getline(inFile,line))
		{
		istringstream ss(line);

		if ( k < 6){
			k+=1;
			continue;
		}
	
		ss >> t_col;    // extracts 1st col, which are the times
		ss >> Re;    // extracts 2nd col, which are the real parts of the amplitude
		ss >> Im; 

		t_vals.push_back(t_col);
		HTeukAmps.push_back(Complex(Re,Im));

		
		}
	inFile.close();

	return HTeukAmps;
}//end of DataLoad

vector<double> TimesLoad(string filename){

	ifstream inFile;
	string line;
	double t_col, Re, Im;
	int iter, input, num_cols;
	vector<double> t_vals;

	inFile.open(filename); 
	if (!inFile)
		{
		cout << "Unable to open Waveform data file \n";
		exit(1); // terminate with error
		}

	vector<complex<double>> HTeukAmps; //HTeuk values 
	int k = 0;
	// loop through each orbit (rows in data file)
	while(getline(inFile,line))
		{
		istringstream ss(line);

		if ( k < 6){
			k+=1;
			continue;
		}
	
		ss >> t_col;    // extracts 1st col, which are the times
		ss >> Re;    // extracts 2nd col, which are the real parts of the amplitude
		ss >> Im; 

		t_vals.push_back(t_col);
		HTeukAmps.push_back(Complex(Re,Im));

		
		}
	inFile.close();

	return t_vals;
}//end of TimesLoad

vector<string> StringBuilder(int Comp, int returntype){
// 0 means compare -w -n with -t -n, 1 means compare -n -t with -f -t

double dp = 1; // make 1
double de = 0.1;
int pmax = 12; // make 14 
int pmin = 11; // make 8 
double emax = 0.7; // This will be rounded to 0.7
double emin = 0.2; // make 0.1
double in_val, in_val1;
int str_offset0 = 1;
int str_offset1 = 1;

string in_str, sub, sub1;
string NITfile,Fullfile, outNITfile, outFullfile, newwaveformFile;
string outwaveformNITFile,outwaveformFullFile,outTwaveformNITFile,outTwaveformFullFile;
string outputFile, outputFile1;
string append_strNIT, append_strFull;
string NITfileSave, FullfileSave;

 // Initial strings
string NITfile0 = "NIT_inspiral -n0 10 0.1 0.001";
string Fullfile0 = "NIT_inspiral -f0 10 0.1 0.001";

string outNITfile0 = "output/Inspiral_NIT_p12_e0.7_q0.001.dat";
string outFullfile0 = "output/Inspiral_Full_p12_e0.7_q0.001.dat";

string outwaveformNITFile0 = "output/Waveform_NIT_p12_e0.7_q0.001.dat";
string outwaveformFullFile0 = "output/Waveform_Full_p12_e0.7_q0.001.dat";

string outTwaveformNITFile0 = "output/T_Waveform_NIT_p12_e0.7_q0.001.dat";
string outTwaveformFullFile0 = "output/T_Waveform_Full_p12_e0.7_q0.001.dat";

string outputFile0 = "output/OverlapData_p10_e0.1_q0.001.dat";
string outputFile01 = "output/OverlapData_p8_e0.1_q0.001.dat";


NITfile = NITfile0;
outNITfile = outNITfile0;
outwaveformNITFile = outwaveformNITFile0;
outTwaveformNITFile = outTwaveformNITFile0;

Fullfile = Fullfile0;
outFullfile = outFullfile0;

outwaveformFullFile = outwaveformFullFile0;
outTwaveformFullFile = outTwaveformFullFile0;

outputFile = outputFile0;
outputFile1 = outputFile01;



vector<string> NITStringsList, FullStringsList, outNITStringsList, outFullStringsList;
vector<string> OutputFiles, WaveformOutputFilesNIT, WaveformOutputFilesFull, WaveformCommands;

for ( int i = 0; i <= floor((pmax - pmin)*(1/dp)); i++) 
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


	NITfile.replace(17,2,sub);
	outNITfile.replace(21,2,sub1);
	
	Fullfile.replace(17,2,sub);
	outFullfile.replace(22,2,sub1);

	outwaveformNITFile.replace(21,2,sub1);
	outwaveformFullFile.replace(22,2,sub1);
	outTwaveformNITFile.replace(23,2,sub1);
	outTwaveformFullFile.replace(24,2,sub1);
	
	outputFile.replace(20,2,sub1);	
	

	for ( int j = 0; j <= floor((emax - emin)*(1/de)); j++) 
	{	
		in_val = emin + j*de;
		in_str = to_string(in_val);
		sub = in_str.substr(0,3);

		NITfile.replace(22,3,sub);
		if(in_val1 < 10){
			outNITfile.replace(25-str_offset0,3,sub);
		}else{
			outNITfile.replace(25,3,sub);	
		}

		Fullfile.replace(22,3,sub);

		if (in_val1 < 10){
			outFullfile.replace(26-str_offset1,3,sub);
		}
		else{
			outFullfile.replace(26,3,sub);
		}
	
		if ( in_val1 < 10){
			outwaveformNITFile.replace(25-str_offset0,3,sub);
			outwaveformFullFile.replace(26-str_offset0,3,sub);
			outTwaveformNITFile.replace(27-str_offset0,3,sub);
			outTwaveformFullFile.replace(28-str_offset0,3,sub);
		}
		else {
			outwaveformNITFile.replace(25,3,sub);
			outwaveformFullFile.replace(26,3,sub);
			outTwaveformNITFile.replace(27,3,sub);
			outTwaveformFullFile.replace(28,3,sub);
		}

		if ( in_val1 < 10){

			outputFile.replace(24-str_offset1,3,sub);

		} 
		else{

			outputFile.replace(24,3,sub);

		}

		NITfileSave = NITfile;
		FullfileSave = Fullfile;

	NITStringsList.push_back(NITfile);
	outNITStringsList.push_back(outNITfile);
	OutputFiles.push_back(outputFile);
	
	FullStringsList.push_back(Fullfile);
	outFullStringsList.push_back(outFullfile);

	if (Comp == 0){
		WaveformOutputFilesNIT.push_back(outwaveformNITFile);
		WaveformOutputFilesFull.push_back(outTwaveformNITFile); 
	} else if (Comp == 1){
		WaveformOutputFilesNIT.push_back(outwaveformFullFile);
		WaveformOutputFilesFull.push_back(outTwaveformFullFile); 
	} else if (Comp == 2){
		WaveformOutputFilesNIT.push_back(outTwaveformNITFile);
		WaveformOutputFilesFull.push_back(outTwaveformFullFile);
	} else if (Comp == 3){
		WaveformOutputFilesNIT.push_back(outwaveformNITFile);
		WaveformOutputFilesFull.push_back(outwaveformFullFile);
	}

		//waveform command string:
		if (Comp == 0){
			newwaveformFile = NITfile.replace(14,2,"w ");
			newwaveformFile.insert(newwaveformFile.size()," -n");
			WaveformCommands.push_back(newwaveformFile);
			NITfile = NITfileSave; // This line refreshes the changes to NITfile to the version before they were modified to make newwaveformFile

			newwaveformFile = NITfile.replace(14,2,"t ");
			newwaveformFile.insert(newwaveformFile.size()," -n");
			WaveformCommands.push_back(newwaveformFile);
			NITfile = NITfileSave;

		} else if (Comp == 1){
			newwaveformFile = Fullfile.replace(14,2,"w ");
			newwaveformFile.insert(newwaveformFile.size()," -f");
			WaveformCommands.push_back(newwaveformFile);
			Fullfile = FullfileSave; // This line refreshes the changes to NITfile to the version before they were modified to make newwaveformFile

			newwaveformFile = Fullfile.replace(14,2,"t ");
			newwaveformFile.insert(newwaveformFile.size()," -f");
			WaveformCommands.push_back(newwaveformFile);
			Fullfile = FullfileSave;
		} else if (Comp == 2){
			newwaveformFile = NITfile.replace(14,2,"t ");
			newwaveformFile.insert(newwaveformFile.size()," -n");
			WaveformCommands.push_back(newwaveformFile);
			NITfile = NITfileSave; // This line refreshes the changes to NITfile to the version before they were modified to make newwaveformFile

			newwaveformFile = Fullfile.replace(14,2,"t ");
			newwaveformFile.insert(newwaveformFile.size()," -f");
			WaveformCommands.push_back(newwaveformFile);
			Fullfile = FullfileSave;
		} else if (Comp == 3){
			newwaveformFile = NITfile.replace(14,2,"w ");
			newwaveformFile.insert(newwaveformFile.size()," -n");
			WaveformCommands.push_back(newwaveformFile);
			NITfile = NITfileSave; // This line refreshes the changes to NITfile to the version before they were modified to make newwaveformFile

			newwaveformFile = Fullfile.replace(14,2,"w ");
			newwaveformFile.insert(newwaveformFile.size()," -f");
			WaveformCommands.push_back(newwaveformFile);
			Fullfile = FullfileSave;	


			}

		} 

	

	string in_str, sub, sub1, append_strNIT, append_strFull, newwaveformFile;
	Fullfile = Fullfile0;
	NITfile = NITfile0;
	outFullfile = outFullfile0;
	outNITfile = outNITfile0;
	outputFile = outputFile0;

	outwaveformNITFile = outwaveformNITFile0;
	outwaveformFullFile = outwaveformFullFile0;
	outTwaveformNITFile = outTwaveformNITFile0;
	outTwaveformFullFile = outTwaveformFullFile0;

	}


vector<string> out;
if (returntype == 0){
	out = NITStringsList;
} else if (returntype == 1){
	out = FullStringsList;
} else if (returntype == 2){
	out = outNITStringsList;
} else if (returntype == 3){
	out = outFullStringsList;
} else if (returntype == 4){
	out = OutputFiles;
} else if (returntype == 5){
	out = WaveformCommands;
}else if (returntype == 6){
	out = WaveformOutputFilesNIT;
}else if (returntype == 7){
	out = WaveformOutputFilesFull;
}

return out;

 
}