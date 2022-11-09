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
#include <string>
#ifdef NETCDFOUTPUT
  #include <netcdf.h>
#endif

 

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace Cantera;
using namespace std;



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
  H = 0, H2O, OH, H2, O, CH, C, CH2, CH3, CH4, C2, C2H, C2H2, C2H3, C2H4, C2H5, C2H6, C4H2, CO, CO2, CH2OH, H2CO, HCO, CH3OH, CH3O, CH3CO, O2, H2CCO, HCCO
};

int main()
{


//Kinetics Network: VULCAN C-H-O
auto sol = newSolution("vulcan_Tsai2017.yaml");
auto gas = sol->thermo();
auto gas_kin = sol->kinetics();

//Gas Properties, indices for I/O storage
int nsp = gas->nSpecies();
fstream InFile;
int nSize = 0;
string data1, data2, data3;
InFile.open("/data4/ananyo/models/C3M/tests/VULCAN_atm_HD209_Kzz.txt"); 
getline(InFile, data1);
getline(InFile, data1);
while (getline(InFile, data1))
nSize++;
InFile.close();

int nTime = 10;
int iTemp = 0;
int iPress =  1;
int iKzz =  2;
int iAlt = 3;

//Initial condition for mole fraction
VectorXd mole_fractions = VectorXd::Zero(nsp);
mole_fractions(H2O) = 9.67E-6;
mole_fractions(H2) = 0.87;
mole_fractions(CO) = 1E-6;
mole_fractions(CH4) = 1E-5;
mole_fractions(O2) = 0;



double Temp, Press, Kzz;
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

//Input from txt file
int inx = 0;
InFile.open("/data4/ananyo/models/C3M/tests/VULCAN_atm_HD209_Kzz.txt"); 
getline(InFile, data1);
getline(InFile, data1);
while (InFile >> data1 >> data2 >> data3){
      AtmData(iPress, inx) = atof(data1.c_str())*0.1; //Pressure (dyn/cm^2) -> Pa
      AtmData(iTemp,inx) = atof(data2.c_str()); //Temperature (K)
      AtmData(iKzz, inx) = atof(data3.c_str())*1E-4; //Kzz (cm^2/s) -> m^2/s
inx++;
}
InFile.close();

//Initial condition for mole fractions
for (int i = 0; i < nSize; i++) {
      b.col(i) = mole_fractions;
    }


//Chemical evolution
double dt = 1e-7;
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
    Press = AtmData(iPress,j); 
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
//Integration for each species (RK4)
    dQ = mat2;
    b.col(j) = dQ + b.col(j) - (flux1*dt);

//Printing the result at each step  
    std::cout << std::endl;
    std::cout << "At time step t = ";
    std::cout << i*dt << " (s) ";
    std::cout << "mole fractions = ";
    std::cout << b.col(j).transpose() << " ";
    std::cout << std::endl;
  
}}

//Writing NetCDF output file
#ifdef NETCDFOUTPUT
#define NETCDFOUTPUT
stringstream msg;
int ifile;
string fname = "Vulcan_CHO_output.nc";
int status = nc_create(fname.c_str(), NC_NETCDF4, &ifile);
std::cout << status << std::endl;
#endif

/*
//Writing NetCDF output file
//Output dimensions: height, nSpecies
//Description: vertical distribution of the chemical species at the time step.
NcFile dataFile("sfc_pres_temp.nc", NcFile::Replace);
NcDim *altDim, *chemDim;
altDim = dataFile.add_dim("altitude", nSize);
chemDim = dataFile.add_dim("chemicals", nSpecies);
altitude = dataFile.add_var("Height", ncFloat, altDim);
altitude->add_att("units", "km");
NcVar ChemConc;
//Either one can manually select species that will be included in the netCDF file or the complete matrix is included
output = dataFile.add_var("mole_fraction", ncFloat, chemDim, altDim);
output->put(&b[0][0], nSpecies, nSize);
*/
}

