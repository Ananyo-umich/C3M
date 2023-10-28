// @sec3{Include files}
// Author: Ananyo Bhattacharya
// Affiliation: University of Michigan
// Email: ananyo@umich.edu
// Solution class describes a phase consists of a mixture of chemical species
#include <cantera/base/Solution.h>

// ThermoPhase object stores the thermodynamic state
#include <cantera/thermo.h>

// Kinetics object stores the chemical kinetics information
#include <cantera/kinetics.h>

// output stream
#include <iostream>
#include <fstream>
#include <cantera/numerics/eigen_dense.h>
#include <cantera/numerics/eigen_sparse.h>
#include <vector>
#include <string>
#include <sstream>
#include <regex>

// Athena++ header
#include <parameter_input.hpp>

// C3M header
#include <PhotoChemistry.hpp>
#include <RadTran.hpp>
#include <interpolation.hpp>


// NetCDF Output
#if NETCDFOUTPUT
  #include <netcdf.h>
#endif

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace Cantera;
using namespace std;


int main(int argc, char **argv) {

// Define the size of the matrix
    int size = 5; // Change this to your desired matrix size

    // Define vectors for the diagonals
    Eigen::VectorXd mainDiagonal(size);
    Eigen::VectorXd upperDiagonal(size - 1);
    Eigen::VectorXd lowerDiagonal(size - 1);

    // Fill the diagonals with values (example values)
    for (int i = 0; i < size; i++) {
        mainDiagonal(i) = 2.0; // Main diagonal
        if (i < size - 1) {
            upperDiagonal(i) = -1.0; // Upper diagonal
            lowerDiagonal(i) = -1.0; // Lower diagonal
        }
    }
    

    // Create the tridiagonal matrix
    Eigen::SparseMatrix<double> tridiagonal(size, size);
    tridiagonal.reserve(3 * size - 2); // Reserve memory for non-zero entries

    // Fill the matrix with the diagonals
    for (int i = 0; i < size; i++) {
        tridiagonal.insert(i, i) = mainDiagonal(i); // Main diagonal
        if (i < size - 1) {
            tridiagonal.insert(i, i + 1) = upperDiagonal(i); // Upper diagonal
            tridiagonal.insert(i + 1, i) = lowerDiagonal(i); // Lower diagonal
        }
    }

    // Finalize the matrix
    tridiagonal.makeCompressed();

    // Print the tridiagonal matrix (for demonstration)
    std::cout << "Tridiagonal Matrix:" << std::endl;
    std::cout << tridiagonal << std::endl;

}