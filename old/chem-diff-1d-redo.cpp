/*
This function solves chemical and diffusion together.  The diffusion is handled first
using the tridiagonal matrix solver for every chemcial species (easily parallelizable).
Then, the chemistry is solved for every layer (also parallelizable).
*/

// @sec3{Include files}
// Solution class describes a phase consists of a mixture of chemical species
#include <cantera/base/Solution.h>

// ThermoPhase object stores the thermodynamic state
#include <cantera/thermo.h>

// Kinetics object stores the chemical kinetics information
#include <cantera/kinetics.h>

// output stream
#include <iostream>
// output file stream:
#include <fstream>
#include <cantera/numerics/eigen_dense.h>
#include <cantera/numerics/eigen_sparse.h>
#include <vector>
#include <numeric>
#include  <string>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// defines Real
#include <configure.hpp>

#define SIZE 300 // defining the T,P arrays to be 500 layers initially.

using namespace std;
using namespace Eigen;
using namespace Cantera;
// using std::vector;
// using Eigen::MatrixBase;


// If using CO (full network)
enum {
  CO = 0, H2, H2O, O2, H2O2, CH4, H2CO, CH3OH, CO2, CH3OOH, C2H2,
    C2H4, C2H6, CH2CO, CH3CHO, C2H5OH, C2H5OOH, CH3COOOH, cC2H4O, C, CH,
    _1CH2, _3CH2, O3P, H, OH, OOH, CH3, HCO, CH2OH, CH3O, CH3OO, C2H, C2H3,
    C2H5, CHCO, CH2CHO, CH3CO, C2H5O, C2H4OOH, C2H5OO, CH3COOO, CH3OCO,
    CO2H, _1C2H4OH, _2C2H4OH, C3H8, C4H8Y, C4H10, C2H5CHO, C3H7OH, C3H7O,
    C4H9O, C2H6CO, C3H8CO, C2H3CHO, _1C3H7, _2C3H7, _1C4H9, _2C4H9, O1D, N2D,
    NO3, HONO2, CH3ONO, CH3NO2, HNO2, CH3NO, NO2, HONO, HCNN, HCNO, N2O,
    NCO, HNO, HOCN, NNH, H2CN, N4S, CN, HNCO, NO, NH, NH2, HCN, NH3, N2,
    N2O4, N2O3, N2H2, N2H3, N2H4, HNNO, HNOH, HNO3, NH2OH, H2NO, CNN, H2CNO,
    C2N2, HCNH, HON, NCN, He, CH3NH2, CH3NH, CH2NH2, CH2NH
};

// enum {
//     CH3OH, H2, C, CO2H, H2O, He, CH3, C2H6, C2H, C2H5OOH, CO, C2H5, C2H5OO,
//     _1CH2, HCO, C2H2, CO2, CH3COOOH, H2O2, CH3CHO, H, CH, C2H5O, OH, CH4,
//     H2CO, C2H4, CH3O, C2H3, CH3CO, CH2OH, _3CH2
// };

// 35 species and 110 reactions (CHO_Wang_red-35-sp-110-rxn.yaml)
// enum {
//     CH3CO, C2H6, CH3O, _3CH2, C2H4, C2H5, CO2H, H2O2, C2H5OH, CH2CHO, OOH,
//     OH, He, CH2CO, H2CO, O3P, O2, H, CO2, CH3OH, CH4, _1C2H4OH, C2H5O, H2O,
//     _1CH2, HCO, O1D, CH3OO, CH2OH, CH3, CH3OOH, CH3CHO, H2, _2C2H4OH, CO
// };


MatrixXd openData(string fileToOpen)
{
    // the inspiration for creating this function was drawn from here
    // https://stackoverflow.com/questions/34247057/how-to-read-csv-file-and-assign-to-eigen-matrix
    vector<double> matrixEntries;
 
    // in this object we store the data from the matrix
    ifstream matrixDataFile(fileToOpen);
 
    // this variable is used to store the row of the matrix that contains commas 
    string matrixRowString;
 
    // this variable is used to store the matrix entry;
    string matrixEntry;
 
    // this variable is used to track the number of rows
    int matrixRowNumber = 0;
 
    while (getline(matrixDataFile, matrixRowString)) // here we read a row by row of matrixDataFile and store every line into the string variable matrixRowString
    {
        stringstream matrixRowStringStream(matrixRowString); //convert matrixRowString that is a string to a stream variable.
 
        while (getline(matrixRowStringStream, matrixEntry, ',')) // here we read pieces of the stream matrixRowStringStream until every comma, and store the resulting character into the matrixEntry
        {
            matrixEntries.push_back(stod(matrixEntry));   //here we convert the string to double and fill in the row vector storing all the matrix entries
        }
        matrixRowNumber++; //update the column numbers
    }
 
    // here we convet the vector variable into the matrix and return the resulting object, 
    // note that matrixEntries.data() is the pointer to the first memory location at which the entries of the vector matrixEntries are stored;
    return Map<Matrix<double, Dynamic, Dynamic, RowMajor>>(matrixEntries.data(), matrixRowNumber, matrixEntries.size() / matrixRowNumber);
}

void saveData(string fileName, MatrixXd  matrix)
{
    //https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
    const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
 
    std::ofstream file(fileName);
    if (file.is_open())
    {
        file << matrix.format(CSVFormat);
        file.close();
    }
}

// Block tridiagonal solver:
// https://github.com/barkm/tridiagonal/blob/master/tridiagonal/tridiagonal_solver.cpp
vector<Eigen::MatrixXd> solve_tridiagonal(vector<Eigen::MatrixXd> &lower_diagonal,
                            vector<Eigen::MatrixXd>  & diagonal,
                            vector<Eigen::MatrixXd>  & upper_diagonal,
                            vector<Eigen::MatrixXd>  & rhs) {

    unsigned long n = diagonal.size();

    vector<Eigen::MatrixXd> upper_diagonal_prime;
    vector<Eigen::MatrixXd> rhs_prime;
    vector<Eigen::MatrixXd> solution;

    upper_diagonal_prime.resize(n - 1);
    rhs_prime.resize(n);
    solution.resize(n);

    for(int i = 0; i < n - 1; i++){
        if(i == 0){
            // C_0' = B_0.inv() * C_0
            upper_diagonal_prime[0] = diagonal[0].inverse() * upper_diagonal[0];
        }
        else{
            // C_i' = (B_i - A_{i-1} * C_{i-1}').inv() * C[i]
            upper_diagonal_prime[i] = (diagonal[i] - lower_diagonal[i-1] * upper_diagonal_prime[i-1]).inverse()
                                        * upper_diagonal[i];
        }
    }
    // D_0' = B_0.inv() * D_0
    rhs_prime[0] = diagonal[0].inverse() * rhs[0];
    for(int i = 1; i < n; i++){
        // D_i' = (B_i.inv() - A_{i-1} * C_{i-1}').inv() * (D_i - A_{i-1} * D_{i-1}')
        rhs_prime[i] = (diagonal[i] - lower_diagonal[i-1] * upper_diagonal_prime[i-1]).inverse()
                        * (rhs[i] - lower_diagonal[i-1] * rhs_prime[i-1]);
    }
    // X_{n-1} = D_{n-1}'
    solution[solution.size() - 1] = rhs_prime[diagonal.size()-1];
    for(int i = n - 2; i >= 0; i--){
        // X_i = D_{n-1}' - C_i' * X_{i+1}
        solution[i] = rhs_prime[i] - upper_diagonal_prime[i] * solution[i+1];
    }
    return solution;
}

// @sec3{Main program}
int main(int argc, char **argv)
{   
    //Read T,P profile where P is in BARS.
    std::ifstream infile("TP_prof_FP_Saturn_100.txt");
    double Temperature[SIZE], Pressure[SIZE], Molardensity[SIZE], dZ[SIZE];
    std::string line;
    int i=0;
    if(infile.is_open()) 
    {
        while(!infile.eof())
        {
            std::getline(infile, line);
            std::istringstream iss(line);
            double Temp, Press, mol_den, delta_z;
            if (!(iss >> Temp >> Press >> mol_den >> delta_z)) { break; }
            Temperature[i] = Temp; // K
            Pressure[i]    = Press; // bars
            Molardensity[i]= mol_den; // kmol/m3
            dZ[i]          = delta_z; // meters
            i++;
        }
        infile.close();
    }

    // matrix to be loaded from a file
    MatrixXd Start_Matrix;
  std::cout << "Loading the reaction network" << std::endl;
    auto sol     = newSolution("CNOH_reaction_network.yaml", "reaction network");
    // auto sol     = newSolution("CHO_no_N_O2_3P_s32_r73_mod.yaml", "CHO_no_N_O2_O3P");
    // auto sol     = newSolution("CHO_Wang_red-35-sp-110-rxn.yaml", "CHO_Wang_red_35");
    auto gas     = sol->thermo();
    auto gas_kin = sol->kinetics();
    
 
    // load the matrix from the file - initial state of system
    // Can be Equilibrium_Matrix, but can also just be the final state of a previous run:
    // Start_Matrix = openData("../../Disequilibria_along_adiabat/HPO_redux_Saturn_30xSol_O_1D_1000steps_1e5.csv");
    // Start_Matrix = openData("../../../unit_test_diff_gauss.csv");
    
    // Start_Matrix = openData("../../Equilibria_along_adiabat/CHO_no_N_O2_O3P_Equilibrium_Saturn_21xSol_O.csv");
    
    // Start_Matrix = openData("../../Equilibria_along_adiabat/CHO_Wang_red_Equilibrium_Saturn_21xSol_O.csv");
    // Start_Matrix = openData("../../Equilibria_along_adiabat/CNHO_He_Equilibrium_Saturn_21xSol_O_normed.csv");
    Start_Matrix = openData("CNHO_He_Equilibrium_Saturn_21xSol_O_normed.csv");

    int layers = Start_Matrix.cols();
    int nsp      = Start_Matrix.rows();
    std::cout << "Number of atmospheric layers: " << layers << std::endl;
    std::cout << "Number of species layers: " << nsp << std::endl;
    std::cout << "Minimum of Equilibrium_Matrix: " << Start_Matrix.minCoeff() << std::endl;
    int CO_index = gas->speciesIndex("CO");
    std::cout << "CO index: " << CO_index << std::endl;

    double dz_mean = 6000;//delta_z[layers - 1]; // last dz
    double K_eddy = 1E5; // m2/s
    double diff_const = K_eddy / pow(dz_mean, 2);
    // double dt = 1/diff_const; // seconds
    double dt = 1E-8; // seconds - for log stepping
    std::cout << "Mean diffusion constant (no dt) = " << diff_const << std::endl;
    std::cout << "dt = " << dt << std::endl;

    double t_tot = 0.0;
    // double max_Time = 2000*dt; // seconds - testing with 100 seconds
    double max_Time = 1E9; // seconds - for log stepping
    double r_log_incr = 1.2; // used for "log"-timestepping.

    MatrixXd Input_Mole_Matrix = Start_Matrix;
    MatrixXd auxilary_mat = Start_Matrix;
    //Initiating Matrices for Time Evolution
    Eigen::SparseMatrix<double>  m_wjac;
    m_wjac.resize(nsp, nsp);
    VectorXd m_wdot(nsp);


    // Initialize the system:
    double inv_dt = 1/dt;
    
    // Define RHS without top or bottom layers
    vector<Eigen::MatrixXd> RHS(layers - 2);

    // Define BC-vector without top or bottom layers (only index 0 and -1 will be non-zero)
    vector<Eigen::MatrixXd> BC(layers - 2);
    
    // Band of I/dt - J (without top and bottom)
    vector<Eigen::MatrixXd> ImJ_diagonal(layers - 2), Main_diagonal(layers - 2);
    
    // Bands for off-diagonals (these will be 1 less in size than the main diagonal):
    vector<Eigen::MatrixXd> Upper_diagonal(layers - 3), Lower_diagonal(layers - 3);
    
    // Generally used identity matrices and main diffusion diagonal:
    MatrixXd I(nsp, nsp);
    I = MatrixXd::Identity(nsp, nsp);

    vector<Eigen::MatrixXd> Diff_diagonal_block(layers - 2);

    // Steady-state check for time loop:
    MatrixXd stead_state_check(layers, 1);
    MatrixXd prev_step(layers, 1);
    double CO_Eq_Min = Start_Matrix.row(CO_index).minCoeff();
    double delta_X = 1.0; // initialize greater than threshold
    int step_count = 0;
    int sub_step_count = 0;
    int max_Steps = 426;

    // Uniform spacing - setting diagonals:
    for (int p = 0; p < layers - 3; p++){
        Lower_diagonal[p] = -K_eddy/pow(dz_mean,2) * I;
    }
    // Both diagonals are the same and have no time evolution
    Upper_diagonal = Lower_diagonal;
    for (int q = 0; q < layers - 2; q++){
        Diff_diagonal_block[q] = -2*K_eddy/pow(dz_mean,2) * I;
    }

    Input_Mole_Matrix.col(0) = Input_Mole_Matrix.col(1);

    // Bottom zero-flux as well:
    // Input_Mole_Matrix.col(layers - 1) = Input_Mole_Matrix.col(layers - 2);
    
    // Initialize BC vector:
    for (int j = 1; j < layers - 1; j++){
        int q = j - 1;
        if (q == 0){ // TOP ZERO-FLUX CONDITION
            // set top boundary (q=0) to mole fractions of the (q=1) layer:
            BC[q] = Input_Mole_Matrix.col(j);
        }
        else if (q == layers-3){ // BOTTOM CONSTANT-VALUE CONDITION
            BC[q] = Input_Mole_Matrix.col(j+1);
            // set bottom boundary (q = layers-3) to X's of the j layer (index: layer-2):
            // BC[q] = Input_Mole_Matrix.col(j);
        }
        else{
            // Otherwise, BC is 0:
            BC[q] = VectorXd::Zero(nsp);
        }
    }

    std::cout << "Initializations completed... Entering main loop:" << std::endl;

    while (step_count < max_Steps){
    // while (t_tot < max_Time){
        // Loop defined over layers excluding the boundaries:
        for (int j = 1; j < layers - 1; j++){
            // index for main diagonal
            int p = j - 1;
            gas->setState_TPX(Temperature[j], Pressure[j]*1E5, &Input_Mole_Matrix.col(j)[0]);

            // net production rates in [kmol/m^3/s or kmol/m^2/s]:
            gas_kin->getNetProductionRates(&m_wdot[0]);

            // Initialize the tridiagonal matrix bands (excludes top and bottom)
            ImJ_diagonal[p] = I*inv_dt - (gas_kin->netProductionRates_ddX() / gas->molarDensity());
            Main_diagonal[p] = ImJ_diagonal[p] - Diff_diagonal_block[p];
            // Main_diagonal[p] = I - Diff_diagonal_block;

            // BC must be added to RHS. BC contributes 0 except at the top and bottom of the RHS-vector.
            RHS[p] = ImJ_diagonal[p]*Input_Mole_Matrix.col(j) + (m_wdot / gas->molarDensity()) + diff_const*BC[p];
            // RHS[p] = I*Input_Mole_Matrix.col(j) + alpha_const*BC[p];
        }
        vector<Eigen::MatrixXd> solution = solve_tridiagonal(Lower_diagonal,
                                                        Main_diagonal,
                                                        Upper_diagonal,
                                                        RHS);

        for (int l = 0; l < layers - 2; l++){
            auxilary_mat.col(l+1) = solution[l];
        }
        // Check if CO has gone below minimum of equilibrium value anywhere.  If so, reduce dt by 1/2.
        ////////if (auxilary_mat.row(CO_index).minCoeff() < CO_Eq_Min){
        //    dt = 0.5 * dt;
        //    sub_step_count++;
        //    std::cout << "Halved dt..." << std::endl;
        //}
        //else{
            prev_step = Input_Mole_Matrix.row(CO_index);

            for (int l = 0; l < layers - 2; l++){
                gas->setMoleFractions(&solution[l](0));
                gas->getMoleFractions(&Input_Mole_Matrix.col(l+1)[0]);
            }
            // Otherwise Input_Mole_Matrix.col(0) is still at the initial condition:
            Input_Mole_Matrix.col(0) = Input_Mole_Matrix.col(1);
            
            // BC update:
            // Zero-flux condition for top:
            // Input_Mole_Matrix.col(0) = Input_Mole_Matrix.col(1);
            BC[0] = Input_Mole_Matrix.col(1);
            
            t_tot = t_tot + dt;
            // dt grows till threshold is reached and is then linear with 10x(1/diff_const):
            if (r_log_incr*dt > 10/diff_const){
                dt = dt; //10/diff_const;
            }
            else{
                dt    = r_log_incr*dt;
            }
            
            inv_dt = 1/dt;

            stead_state_check = (prev_step - Input_Mole_Matrix.row(CO_index)) * prev_step.transpose().cwiseInverse();
            delta_X = stead_state_check.cwiseAbs().maxCoeff();
            std::cout << delta_X << ", " << step_count << ", " << dt << ", " << sub_step_count << std::endl;

            step_count++;
            sub_step_count = 0;
        //}
    }
    
    // This is how you set vector elements (even when each element is a matrix!)
    // for (int l=0; l<3; l++){
    //     diagonal.at(l) = b2; // just a concrete example
    // }
    // And this is how you print those elements once the vector is set:
    // for (auto i: diagonal){
    //     std::cout << i << std::endl;
    // }

    std::cout << "Computations completed at (s): " << t_tot << std::endl;
    std::cout << "Total Steps: " << step_count << std::endl;
    std::cout << "Final dt: " << dt << std::endl;
    // saveData("../../../unit_test_diff_gauss_step_1000_BCconst.csv", Input_Mole_Matrix);
    // saveData("../../Disequilibria_along_adiabat/CHO_no_N_O2_O3P_Saturn_redo-test_loglint125_uniform.csv", Input_Mole_Matrix);
    // saveData("../../Disequilibria_along_adiabat/CHO_no_N_O2_O3P_Saturn_redo-test_loglint115.csv", Input_Mole_Matrix);
    // saveData("../../Disequilibria_along_adiabat/CHO_Wang_red_35_Saturn_redo-test_loglint115.csv", Input_Mole_Matrix);

    saveData("CNHO_Wang_Saturn_21xSol_O_1D_uniform_526_mwdot.csv", Input_Mole_Matrix);

    return 0;
}
