//The program solves the kinetics for thermochemistry network in 1-D domain

#include <cantera/base/Solution.h>


// ThermoPhase object stores the thermodynamic state
#include <cantera/thermo.h>

// Kinetics object stores the chemical kinetics information
#include <cantera/kinetics.h>

// output stream
#include <iostream>
#include <cantera/numerics/eigen_dense.h>
#include <cantera/numerics/eigen_sparse.h>
#include <vector>
#include "cantera/oneD/Sim1D.h"
#include "cantera/oneD/Boundary1D.h"
#include "cantera/oneD/StFlow.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/transport.h"
#include "cantera/base/Solution.h"
#include "cantera/base/Array.h"
#include "cantera/base/stringUtils.h"
#include <fstream>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace Cantera;

enum {
  PH3 = 0, H3PO4, HOPO, H, PH2, PH, H2O, H2, O2, O, OH, HO2, PO, PO2,
    PO3, HPO, HOPO2, P2O3, P, P2, P4, P2O, P2O2, HPOH, H2POH, He
};

int main()
{
//Kinetics Network
auto sol = newSolution("P_Reaction_Set.yaml", "p_reaction_set");
auto gas = sol->thermo();
auto gas_kin = sol->kinetics();

//Gas Properties
double *mole_fractions = new double [gas->nSpecies()];
mole_fractions[H2O] = 9.67E-4;
mole_fractions[H2] = 0.87;
mole_fractions[He] = 0.1;
mole_fractions[PH3] = 2.3E-6;

double Temp = 1000.;
double Kzz = 1E6;
gas->setState_TP(Temp, 300.0*OneAtm);
gas->setMoleFractions(mole_fractions);

double *flux_in = new double [gas->nSpecies()];
double *flux_out = new double [gas->nSpecies()];

//Domain
int nsp = gas->nSpecies();
int nSize = 100;
int nTime = 10000;
int iTemp = 0;
int iPress =  1;
int iKzz =  2;
int iAlt = 3;

VectorXd m_wdot(nsp);
VectorXd mole_frac(nsp);
VectorXd Un(nsp);
VectorXd Unext(nsp);
VectorXd Uprev(nsp);
VectorXd diff_flux(nsp);
VectorXd dQ(nsp);
double Keddy(nSize);
MatrixXd AtmData(4, nSize);
MatrixXd b(nsp, nSize);
//The extra dimensions stand for temperature, pressure and flux

//Initial Condition
//The input variables like T, P, Height, eddy diffusion, species can be included from a netCDF
//file
for (int i = 0; i < nSize; i++) {
      Temp = 500.0 + 100.0*i;
      Kzz = 1E12 + 1E11*i;
      AtmData(iAlt, i) = i;
      AtmData(iTemp,i) = Temp; 
      AtmData(iKzz, i) = Kzz;
      gas->getMoleFractions(&mole_frac[0]);
      b.col(i) = mole_frac;
    }


//Chemical evolution
float dt = 1e-8;
float dh = 0.1; //Km
int iPrev, iNext;

//Initiating Matrices for Time Evolution
Eigen::SparseMatrix<double>  m_wjac;
m_wjac.resize(nsp, nsp);
MatrixXd mat1(nsp, nsp);
mat1 = MatrixXd::Identity(nsp, nsp);
MatrixXd mat2(nsp, nsp);


//Iterating over each time step
//Iterating over each grid point
//Iterating over each species

for (int i = 0; i < nTime; i++) {
  for (int j = 1; j < nSize-1; j++) {
//Setting T, P, X for each grid point
    iPrev = j-1;
    iNext = j+1;
    Temp = AtmData(iTemp,j);
    gas->setState_TP(Temp, 300.0*OneAtm);
    mole_frac = b.col(j);
    gas->setMoleFractions(&mole_frac[0]);
    Keddy = AtmData(iKzz, j);
//Solving the net production for each species
    gas_kin->getNetProductionRates(&m_wdot[0]);
    m_wjac = gas_kin->netProductionRates_ddX();
    mat2 = ((mat1/dt) - m_wjac);
    mat2 = mat2.inverse();
    mat2 = mat2*m_wdot; 
//Diffusion terms
    Un = b.col(j);
    Unext = b.col(iNext);
    Uprev = b.col(iPrev);
    diff_flux = Keddy*(Unext - (2*Un) + Uprev)/(dh*dh);
//Integration for each species (Foward Euler simple case)
    dQ = mat2;
    b.col(j) = b.col(j) +  dQ;
    std::cout << std::endl;
    std::cout << "At time step t = ";
    std::cout << i*dt << " (s) ";
    std::cout << "mole fractions = ";
    std::cout << b.col(j).transpose() << " ";
    std::cout << std::endl;
  }


}


//Output

}

