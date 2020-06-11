// -----------------Interpolant of Complex Amplitudes ------------------------- //

#include <ComplexAmpTest.h>

int main()
{

    ifstream inFile;

    inFile.open("A_l2_m2.dat");
    if (!inFile)
    {
        cout << "Unable to open file \n";
        exit(1); // terminate with error
    }

    // ------------------ gather important parameters --------------------------------- //

    int col_count = 83; // This is the total number of columns

    //----------------- read file and create e and x column  -------------------------- //

    vector<double> x, e; // x will be the value of p-2e
    string line;

    while (getline(inFile, line)) // going to eventually have too throw away first line of output for each (title line)
    {
        istringstream ss(line);

        double x_col, e_col;

        ss >> x_col; // extracts 1st col, which is the p-2e column
        ss >> e_col; // extracts 2nd col, which is the e column

        // output the x-values and e-values to check:
        // cout << x_col << " " << e_col << endl; Print e's and x's to verify

        x.push_back(x_col);
        e.push_back(e_col);
    }
    // end of use of the file stream
    inFile.close();

    vector<vector<double>> N_re, N_im; // This needs to stay outside the loop as to not erase it everytime

    for (int num_cols = 1; num_cols < col_count; num_cols++)
    {
        cin.clear(); // This line is SUPER IMPORTANT: This refreshes the flags on input stream allowing the reopening of A_l2_m2.dat

        //----------------- n column loop ( needs to loop over all n)  -------------------------- //

        inFile.open("A_l2_m2.dat"); // This goes back to the directory which contains the file
        if (!inFile)
        {
            cout << "Unable to open file \n";
            exit(1); // terminate with error
        }

        vector<double> n_re, n_im;

        int input;
        while (getline(inFile, line))
        {
            istringstream ss(line);

            complex<double> n_col;

            ss >> n_col; // extracts 1st col, which is the p-2e column
            ss >> n_col; // extracts 2nd col, which is the e column

            int iter = 0;

            // This while loop allows the loading of the columns from the first to the col_count-th column
            while (iter < num_cols)
            {
                ss >> n_col;
                iter++;
            }

            // output the n-values to check:
            // cout << n_col << endl; // verify n's

            n_re.push_back(n_col.real());
            n_im.push_back(n_col.imag());
        }

        N_re.push_back(n_re);
        N_im.push_back(n_im);

        inFile.close();
    }

    // --------------------------------------------------------

    // Get rid of the zeros in the beginning:

    x.erase(x.begin());
    e.erase(e.begin());

    for (int i = 0; i < N_re.size(); i++)
    {
        N_re[i].erase(N_re[i].begin());
        N_im[i].erase(N_im[i].begin()); // N_re and N_im are the same size
    }

    // --------------------------------------------------------

    // Up to now, the last displayed value is the final real value in the final column of the .dat file. This means it worked!
    // We have now read all the data, one column at a time, into the program, only leaving the interpolation of x,e and each column.


    // --------------------- Interpolation ------------------------------------- //

    // Generate the data used for interpolation:

    vector<double> unique_x, unique_e;

    unique_x.push_back(x[0]);
    unique_e.push_back(e[0]);

    // Find the unique x values:
    for (int i = 1; i < x.size(); i++)
    {
        if (x[i] == x[i - 1])
        {
            continue;
        }
        else
        {
            unique_x.push_back(x[i]);
        }
    }

    // Find the unique e values:
    for (int i = 1; i < col_count - 2; i++) // col_count -2 instead of -1 fixes the 0.402 extra line error.
    {
        unique_e.push_back(e[i]);
    }

    vector<double> x_test_pre, e_test_pre;
    double test_iter_x, test_iter_e;

    // Find midpoints of given data:

    vector<double> test_vec;

    double tot = (26.05-6.05)/0.25; 
    double tot_e = (0.802-0.002)/0.01;
    int divisions = 4;
    double inc_x = 0.25/divisions;// the 0.25 increment would be the same as the original x data
    double inc_e = 0.01/divisions;// the 0.01 increment would be the same as the original e data

    for(int i = 0; i < (4*tot)+1;i++)
    {
        if ((i & divisions-1) != 0) // This avoids rewriting the original values, which were already used in training
        {
        x_test_pre.push_back(i*inc_x + unique_x[0]);
        e_test_pre.push_back(i*inc_e + unique_e[0]);
        } 
        
        else
        {
            continue;
        }
    }


// unique_x and unique_e contain the original unique values, so thay can be used with find to remove them from x_test_pre,e_test_pre

// Test vector for finer sampling:
    for (int i = 0; i < x_test_pre.size(); i++)
    {
        cout << "x test pre values: " << x_test_pre[i] << endl; // Index goes from 0 to 80 (captures the data for all 81 n's)
    }

    // From here we can use a double for loop to generate the values of x and e as new columns: x_test and e_test

    vector<double> x_test, e_test;

    for (int i = 0; i < x_test_pre.size() - 1; i++)
    {
        for (int j = 0; j < e_test_pre.size() - 1; j++)
        {
            x_test.push_back(x_test_pre[i]);
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
        N_Interp_re.push_back(new Interpolant(x, e, N_re[i]));
        N_Interp_im.push_back(new Interpolant(x, e, N_im[i]));

        // Interpolated data part:
        vector<double> a_re, a_im;

        for (int k = 0; k < x_test.size(); k++)
        {
            a_re.push_back(N_Interp_re[i]->eval(x_test[k], e_test[k]));
            a_im.push_back(N_Interp_im[i]->eval(x_test[k], e_test[k]));
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
    for (int i = 0; i < x_test.size(); i++)
    {
        // Add the test values to the beginning of the row.
        output << x_test[i] << " " << e_test[i];

        // Loop through the columns and add each harmonic mode.
        for (int j = 0; j < A_re.size(); j++)
        {
            output << " (" << A_re[j][i] << "," << A_im[j][i] << ")";
        }

        output << endl;

    }

    // --------------------------- End ----------------------------------------- //

} // This is the brace for int main; end of the program