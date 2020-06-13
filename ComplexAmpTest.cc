// -----------------Interpolant of Complex Amplitudes ------------------------- //

#include <fstream>
#include <iostream>
#include <complex>
#include <iomanip>
#include <algorithm>
#include "Interpolant.h"

int main()
{

//----------------------------
//  read Almn data from file
//----------------------------
ifstream inFile;

// initialize variables outside loops (speeds up by avoiding extra initialization cost)
string line, temp;
double y_col, e_col, test;
int iter, input, num_cols;
vector<double> ys, es; //y = p-2*e
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
	while(ss >> test) col_count++;
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
	es.push_back(e_col);
	
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
// We have now read all the data, one row at a time, into the program, only leaving the interpolation 


    // --------------------- Interpolation ------------------------------------- //

    // Generate the data used for interpolation:

    vector<double> unique_y, unique_e;

    unique_y.push_back(ys[0]);
    unique_e.push_back(es[0]);

    // Find the unique x values:
    for (int i = 1; i < ys.size(); i++)
    {
        if (ys[i] == ys[i - 1])
        {
            continue;
        }
        else
        {
            unique_y.push_back(ys[i]);
        }
    }

    // Find the unique e values:
    for (int i = 1; i < col_count - 2; i++) // col_count -2 instead of -1 fixes the 0.402 extra line error.
    {
        unique_e.push_back(es[i]);
    }

    vector<double> y_test_pre, e_test_pre;
    double test_iter_y, test_iter_e;

    // Find midpoints of given data:

    vector<double> test_vec;

    double tot = (26.05-6.05)/0.25; 
    double tot_e = (0.802-0.002)/0.01;
    int divisions = 4;
    double inc_y = 0.25/divisions;// the 0.25 increment would be the same as the original x data
    double inc_e = 0.01/divisions;// the 0.01 increment would be the same as the original e data

    for(int i = 0; i < (4*tot)+1;i++)
    {
        if ((i & divisions-1) != 0) // This avoids rewriting the original values, which were already used in training
        {
        y_test_pre.push_back(i*inc_y + unique_y[0]);
        e_test_pre.push_back(i*inc_e + unique_e[0]);
        } 
        
        else
        {
            continue;
        }
    }


// unique_x and unique_e contain the original unique values, so thay can be used with find to remove them from x_test_pre,e_test_pre

// Test vector for finer sampling:
    for (int i = 0; i < y_test_pre.size(); i++)
    {
        cout << "x test pre values: " << y_test_pre[i] << endl; // Index goes from 0 to 80 (captures the data for all 81 n's)
    }

    // From here we can use a double for loop to generate the values of x and e as new columns: x_test and e_test

    vector<double> y_test, e_test;

    for (int i = 0; i < y_test_pre.size() - 1; i++)
    {
        for (int j = 0; j < e_test_pre.size() - 1; j++)
        {
            y_test.push_back(y_test_pre[i]);
            e_test.push_back(e_test_pre[j]);
        }
    }

    // ---------------------- Interpolation ---------------------------------- //

    // Do a 2d interpolation of the data using x,e and n for each n.

    // Construct a 2D interpolant of f(x,y)
    // Interpolant::Interpolant(Vector x, Vector y, Vector f)

    vector<Interpolant*> N_Interp_re, N_Interp_im;
    vector<vector<double>> A_re, A_im;

    for (int i = 0; i < col_count - 2; i++)
    {
        // Interpolant part:
        N_Interp_re.push_back(new Interpolant(ys, es, N_re[i]));
        N_Interp_im.push_back(new Interpolant(ys, es, N_im[i]));

        // Interpolated data part:
        vector<double> a_re, a_im;

        for (int k = 0; k < y_test.size(); k++)
        {
            a_re.push_back(N_Interp_re[i]->eval(y_test[k], e_test[k]));
            a_im.push_back(N_Interp_im[i]->eval(y_test[k], e_test[k]));
        }

        A_re.push_back(a_re);
        A_im.push_back(a_im);
    }

    // --------------- Output Interpolation results to file -------------------- //

    // This is a test that the data has been stored in A_re and A_im (Still needs to be written to an output file)
    // Write A_re and A_im to output file as members of a complex double (same way as given data):

    // Output controls stream to file 'A_l2_m2_interp_test.dat'.
    ofstream output;
    output.open("A_l2_m2_interp_test.dat");

    // Output the column headers
    output << "p-2*e e";

    for (int i = -40; i <= 40; i++)
    {
        output << " " << i;
    }

    output << endl;

    // Loop over test rows
    for (int i = 0; i < y_test.size(); i++)
    {
        // Add the test values to the beginning of the row.
        output << y_test[i] << " " << e_test[i];

        // Loop through the columns and add each harmonic mode.
        for (int j = 0; j < A_re.size(); j++)
        {
            output << " (" << A_re[j][i] << "," << A_im[j][i] << ")";
        }

        output << endl;

    }
	
// close file
output.close();

// free memory to prevent memory leaks in future applications
for(iter = 0; iter < col_count-2; iter++)
	{
	delete N_Interp_re[iter];
	delete N_Interp_im[iter];
	}

    // --------------------------- End ----------------------------------------- //

} // This is the brace for int main; end of the program