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

// External library headers
#ifdef NETCDFOUTPUT
  #include <netcdf.h>
#endif

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace Cantera;

//Reading netCDF file for input parameters i.e. T, P, Kzz, Altitude
void NCFileRead(std::string fname)
{
#ifdef NETCDFOUTPUT
  int fileid, dimid, varid, err;
  char tname[80];
  nc_open(fname.c_str(), NC_NETCDF4, &fileid);

  nc_inq_dimid(fileid, "Wavenumber", &dimid);
  nc_inq_dimlen(fileid, dimid, (size_t*)len_);
  nc_inq_dimid(fileid, "Pressure", &dimid);
  nc_inq_dimlen(fileid, dimid, (size_t*)(len_ + 1));
  strcpy(tname, "T_");
  strcat(tname, name_.c_str());
  nc_inq_dimid(fileid, tname, &dimid);
  nc_inq_dimlen(fileid, dimid, (size_t*)(len_ + 2));

  axis_.resize(len_[0] + len_[1] + len_[2]);

  nc_inq_varid(fileid, "Wavenumber", &varid);
  nc_get_var_double(fileid, varid, axis_.data());
  err = nc_inq_varid(fileid, "Pressure", &varid);
  err = nc_get_var_double(fileid, varid, axis_.data() + len_[0]);
  nc_inq_varid(fileid, tname, &varid);
  nc_get_var_double(fileid, varid, axis_.data() + len_[0] + len_[1]);

  Real *temp = new Real[len_[1]];
  nc_inq_varid(fileid, "Temperature", &varid);
  nc_get_var_double(fileid, varid, temp);

#endif
}



//Function for Central Difference Scheme
MatrixXd diff_flux(int inxn, int inxNext, int inxPrev, MatrixXd Sp, int ns, double Kzz, double dx){
  VectorXd diff_flux(ns), Un(ns), Unext(ns), Uprev(ns); 
  Un = Sp.col(inxn);
  Unext = Sp.col(inxNext);
  Uprev = Sp.col(inxPrev);
  diff_flux = Kzz*(Unext - (2*Un) + Uprev)/(dx*dx);
return diff_flux;
}

//Index of the chemical species in the reaction
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

//Gas Properties, indices for I/O storage
int nsp = gas->nSpecies();
int nSize = 100;
int nTime = 10000;
int iTemp = 0;
int iPress =  1;
int iKzz =  2;
int iAlt = 3;

//Initial condition for mole fraction
VectorXd mole_fractions = VectorXd::Zero(nsp);
mole_fractions(H2O) = 9.67E-4;
mole_fractions(H2) = 0.87;
mole_fractions(He) = 0.1;
mole_fractions(PH3) = 2.3E-6;

double Temp = 1000.;
double Kzz = 1E6;
gas->setState_TP(Temp, 300.0*OneAtm);
gas->setMoleFractions_NoNorm(&mole_fractions[0]);

//Initial condition for boundary fluxes
VectorXd flux_in = VectorXd::Zero(nsp);
VectorXd flux_out = VectorXd::Zero(nsp);


//Domain

VectorXd m_wdot(nsp);
VectorXd mole_frac(nsp);
VectorXd Un(nsp);
VectorXd Unext(nsp);
VectorXd Uprev(nsp);
VectorXd flux1(nsp);
VectorXd flux2(nsp);
VectorXd flux3(nsp);
VectorXd flux4(nsp);
VectorXd dQ(nsp);
double Keddy(nSize);
MatrixXd AtmData(4, nSize);
MatrixXd b(nsp, nSize);
MatrixXd a(nsp, nSize);
//The extra dimensions stand for temperature, pressure and flux

//Initial Condition
//The input variables like T, P, Height, eddy diffusion, species can be included from a netCDF
//file
for (int i = 0; i < nSize; i++) {
      Temp = 500.0 + 100.0*i;
      Kzz = 1E2 + 1E1*i;
      AtmData(iAlt, i) = i;
      AtmData(iTemp,i) = Temp; 
      AtmData(iKzz, i) = Kzz;
      b.col(i) = mole_fractions;
    }


//Chemical evolution
double dt = 1e-8;
double dh = 1000; //m
int iPrev, iNext;

//Initiating Matrices for Time Evolution
Eigen::SparseMatrix<double>  m_wjac;
m_wjac.resize(nsp, nsp);
MatrixXd mat1(nsp, nsp);
mat1 = MatrixXd::Identity(nsp, nsp);
MatrixXd mat2(nsp, nsp);


for (int i = 0; i < nTime; i++) {
  for (int j = 0; j < nSize; j++) {
//Setting T, P, X for each grid point
    iPrev = j-1;
    iNext = j+1;
    Temp = AtmData(iTemp,j);
    gas->setState_TP(Temp, 300.0*OneAtm);
    mole_frac = b.col(j);
    gas->setMoleFractions_NoNorm(&mole_frac[0]);
    Keddy = AtmData(iKzz, j);
//Solving the net production for each species
    gas_kin->getNetProductionRates(&m_wdot[0]); //Extracting net production rates from Cantera
    m_wjac = gas_kin->netProductionRates_ddX(); //Extracting Jacobian from Cantera


//Backward Euler Scheme (Li and Chen, 2019)
    mat2 = ((mat1/dt) - m_wjac);
    mat2 = mat2.inverse();
    mat2 = mat2*m_wdot;
    a = b;
    if (j == 0){
    flux1 = flux_in;
    }
    if (j == nSize-1){
    flux1 = flux_out;
    }
    if ((j > 0) && (j < nSize-1)) {


//Diffusion terms (central difference scheme)
    flux1 = diff_flux(j,iNext,iPrev,a,nsp,Keddy,dh);
    a.col(j) = a.col(j) - (flux1*dt);
    flux2 = diff_flux(j,iNext,iPrev,a,nsp,Keddy,dh);
    a.col(j) = a.col(j) - (flux2*dt/2);
    flux3 = diff_flux(j,iNext,iPrev,a,nsp,Keddy,dh);
    a.col(j) = a.col(j) - (flux3*dt/2);
    flux4 = diff_flux(j,iNext,iPrev,a,nsp,Keddy,dh);
    flux1 = ((flux1 + (2*flux2) + (2*flux3) + flux4)/6);
    }
//Integration for each species (Foward Euler simple case)
    dQ = mat2;
    b.col(j) = dQ + b.col(j) - (flux1*dt);

//Printing the result at each step  
    std::cout << std::endl;
    std::cout << "At time step t = ";
    std::cout << i*dt << " (s) ";
    std::cout << "mole fractions = ";
    std::cout << b.col(j).transpose() << " ";
    std::cout << std::endl;
  
}

}

}

