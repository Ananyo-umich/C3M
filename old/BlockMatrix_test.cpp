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
    int species = 5; // Change this to your desired matrix size
    int alt = 10;
    // Define vectors for the diagonals
    Eigen::VectorXd mainDiagonal(species);
    Eigen::VectorXd upperDiagonal(species - 1);
    Eigen::VectorXd lowerDiagonal(species - 1);


   for (int i = 0; i < species; i++) {
        mainDiagonal(i) = 2.0; // Main diagonal
        if (i < species - 1) {
            upperDiagonal(i) = -1.0; // Upper diagonal
            lowerDiagonal(i) = -1.0; // Lower diagonal
        }
    }


   // Create the block matrix of species x altitude
   Eigen::SparseMatrix<double> m(species*alt,species*alt);
   
   //Fill in the values of the matrix
   for(int i = 0; i < alt; i++)
   {
     for(int j = 0; j < species; j++){
     m.insert(i*species + j, i*species + j) = mainDiagonal(j);
     std::cout << "main diagonal working" << std::endl;
     if(j < species-1){
     m.insert(i*species + j, i*species + j+1) = upperDiagonal(j);
     m.insert(i*species + j+1, i*species + j) = lowerDiagonal(j);
     std::cout << "storage working" << std::endl;
     }


     }

   }

   m.makeCompressed();
   //Printing out the matrix
   std::cout << "Block matrix: " << std::endl;
   std::cout << m << std::endl;
}
