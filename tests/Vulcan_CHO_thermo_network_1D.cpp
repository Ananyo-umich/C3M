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
#include <configure.hpp>

#if NETCDFOUTPUT
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
double solar_O = 6.0618E-4;
double solar_C = 2.7761E-4;
double enrf_O = 1;
double enrf_C = 1;

mole_fractions(H2) = 0.87;
mole_fractions(H2O) = 0.793*solar_O*enrf_O*mole_fractions(H2);
mole_fractions(CH4) = solar_C*enrf_C*mole_fractions(H2);



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
double dt = 1e-5;
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
    gas->setState_TP(Temp, (Press/1.0132E5)*OneAtm);
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
    //dh = log(AtmData(iPress,j+1)/AtmData(iPress,j-1)); 
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
  
} 
    
}

//Storing the variables in intermediate variables
VectorXd cPress = AtmData.row(iPress);
VectorXd cH2O = b.row(H2O);
VectorXd cCO = b.row(CO);
VectorXd cCO2 = b.row(CO2);
VectorXd cCH3 = b.row(CH3);
VectorXd cCH4 = b.row(CH4);
VectorXd cC2H2 = b.row(C2H2);
VectorXd cC2H4 = b.row(C2H4);
VectorXd cC2H6 = b.row(C2H6);
VectorXd cCH3OH = b.row(CH3OH);
VectorXd cH = b.row(H);
VectorXd cO = b.row(O);
VectorXd cH2CO = b.row(H2CO);
//Writing NetCDF output file


#if NETCDFOUTPUT
stringstream msg;
int ifile;
string fname = "Vulcan_CHO_output.nc";
int status = nc_create(fname.c_str(), NC_NETCDF4, &ifile);
// Defining the dimensions
int alt, iPres, iH2O, iCO, iCO2, iCH3, iCH4, iC2H2, iC2H4, iC2H6, iCH3OH, iH, iO, iH2CO;
// Atmospheric Properties
nc_def_dim(ifile, "Pressure (Pa)", nSize, &alt);
nc_def_var(ifile, "Pressure (Pa)", NC_DOUBLE, 1, &alt, &iPres);
nc_put_var_double(ifile, iPres, &cPress[0]);

// Chemical Species
nc_def_var(ifile, "H2O", NC_DOUBLE, 1, &alt, &iH2O);
nc_put_var_double(ifile, iH2O, &cH2O[0]);

nc_def_var(ifile, "CO", NC_DOUBLE, 1, &alt, &iCO);
nc_put_var_double(ifile, iCO, &cCO[0]);

nc_def_var(ifile, "CO2", NC_DOUBLE, 1, &alt, &iCO2);
nc_put_var_double(ifile, iCO2, &cCO2[0]);

nc_def_var(ifile, "CH3", NC_DOUBLE, 1, &alt, &iCH3);
nc_put_var_double(ifile, iCH3, &cCH3[0]);

nc_def_var(ifile, "CH4", NC_DOUBLE, 1, &alt, &iCH4);
nc_put_var_double(ifile, iCH4, &cCH4[0]);

nc_def_var(ifile, "C2H2", NC_DOUBLE, 1, &alt, &iC2H2);
nc_put_var_double(ifile, iC2H2, &cC2H2[0]);

nc_def_var(ifile, "C2H4", NC_DOUBLE, 1, &alt, &iC2H4);
nc_put_var_double(ifile, iC2H4, &cC2H4[0]);

nc_def_var(ifile, "C2H6", NC_DOUBLE, 1, &alt, &iC2H6);
nc_put_var_double(ifile, iC2H6, &cC2H6[0]);

nc_def_var(ifile, "CH3OH", NC_DOUBLE, 1, &alt, &iCH3OH);
nc_put_var_double(ifile, iCH3OH, &cCH3OH[0]);

nc_def_var(ifile, "H", NC_DOUBLE, 1, &alt, &iH);
nc_put_var_double(ifile, iH, &cH[0]);

nc_def_var(ifile, "O", NC_DOUBLE, 1, &alt, &iO);
nc_put_var_double(ifile, iO, &cO[0]);

nc_def_var(ifile, "H2CO", NC_DOUBLE, 1, &alt, &iH2CO);
nc_put_var_double(ifile, iH2CO, &cH2CO[0]);

nc_close(ifile);
#endif

}

