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

int main(int argc, char **argv) {
  
//Reading input file
  IOWrapper infile;
  infile.Open("test_Venus.inp", IOWrapper::FileMode::read);

  ParameterInput *pinput = new ParameterInput();

  pinput->LoadFromFile(infile);
  infile.Close();

//Loading the input parameters for atmospheric profile and reaction network
  std::string atm_file = pinput->GetString("problem", "planet");
  std::string network_file = pinput->GetString("problem", "network");
  
//Reading the chemical kinetics network
  auto sol = newSolution(network_file);
  std::cout << "Hello" << std::endl;
  auto gas = sol->thermo();
  auto gas_kin = sol->kinetics();  
  int nsp = gas->nSpecies();
  int nrxn = gas_kin->nReactions();
  
//Initial condition for mole fraction
  VectorXd mole_fractions = VectorXd::Zero(nsp);
  
//Initial condition for boundary fluxes
  VectorXd flux_lower = VectorXd::Zero(nsp);
  VectorXd flux_upper = VectorXd::Zero(nsp);

//Loading the input parameters for initial condition
  std::string init_species_list = pinput->GetString("init", "species");
  std::regex pattern ("[a-zA-Z0-9_]+");
  std::smatch m;
  int species_inx;
  while (std::regex_search (init_species_list,m,pattern)) {
    for (auto x:m){
    std::string species_init_condition = pinput->GetString("init", x);
    species_inx = gas->speciesIndex(x);
    mole_fractions(species_inx) = atof(species_init_condition.c_str());
    }
    init_species_list = m.suffix().str();
  }


//Loading the input for upper boundary condition
  std::string ub_species_list = pinput->GetString("upperboundaryflux", "species");
  while (std::regex_search (ub_species_list,m,pattern)) {
    for (auto x:m){
    std::string species_upperboundary_condition = pinput->GetString("upperboundaryflux", x);
    species_inx = gas->speciesIndex(x);
    flux_upper(species_inx) = atof(species_upperboundary_condition.c_str());
    }
    ub_species_list = m.suffix().str();
  }

//Loading the input for lower boundary condition 
  std::string lb_species_list = pinput->GetString("lowerboundaryflux", "species");
  while (std::regex_search (lb_species_list,m,pattern)) {
    for (auto x:m){
    std::string species_lowerboundary_condition = pinput->GetString("lowerboundaryflux", x);
    species_inx = gas->speciesIndex(x);
    flux_lower(species_inx) = atof(species_lowerboundary_condition.c_str());
    }
    lb_species_list = m.suffix().str();
  }

//Integrator
  std::string time_step = pinput->GetString("integrator", "dt");
  std::string total_time_steps = pinput->GetString("integrator", "nSteps");
  double dt = atof(time_step.c_str());
  int nTime = stoi(total_time_steps.c_str());


//Photochemistry
  std::string photochem = pinput->GetString("problem", "photochem");
  if(photochem == "true")
  {
  
  
//Radiative transfer
//Reading the Stellar Irradiance Input
  std::string stellar_input_file = pinput->GetString("radtran", "solar");
  MatrixXd stellar_input = ReadStellarRadiationInput(stellar_input_file);  
  
  
//Reading absorber cross section for all absorbers
  std::string absorber_species_list = pinput->GetString("abscross", "absorbers");
  while (std::regex_search (absorber_species_list,m,pattern)) {
    for (auto x:m){
    std::string absorber_cross_info = pinput->GetString("abscross", x);
    std::string absorber_cross_database = absorber_cross_info.substr(0,absorber_cross_info.find(","));  
    std::string absorber_cross_file = absorber_cross_info.substr(absorber_cross_info.find(",")+1,absorber_cross_info.length()-1);
    
//Reading the absorption cross section for absorber as per the database structure
//Interpolating the cross sections to reference grid

   if(absorber_cross_database == "VULCAN"){
     MatrixXd cross_info = ReadVULCANPhotoAbsCrossSection(absorber_cross_file);
     MatrixXd cross_section_data = InterpolateCrossSection(stellar_input.row(0), cross_info.row(0), cross_info.row(1));
   }

/*    
   if(absorber_cross_database == "KINETICS"){
     std::cout << absorber_cross_database << std::endl;

    }
*/
    
   if(absorber_cross_database == "VPL"){
     MatrixXd cross_info = ReadAtmosCrossSection(absorber_cross_file);
   
   } 

    //species_inx = gas->speciesIndex(x);
    
    }
    absorber_species_list = m.suffix().str();
  }
  
//Reading cross section for all photochemical reactions
  for(int irxn = 0; irxn < nrxn; irxn++){
    std::string rxnEquation = gas_kin->reactionString(irxn);
    int pos = rxnEquation.find("=");
    rxnEquation.replace(pos-1, 3, "->");
    std::string photo_cross_info = pinput->GetOrAddString("photocross", rxnEquation, "nan");
    
    if(photo_cross_info != "nan"){
    std::string photo_cross_database = photo_cross_info.substr(0,photo_cross_info.find(","));  
    std::string photo_cross_file = photo_cross_info.substr(photo_cross_info.find(",")+1,photo_cross_info.length()-1);
    
    
    if(photo_cross_database == "VULCAN_Ion"){
     MatrixXd photoXsection_info = ReadVULCANPhotoIonCrossSection(photo_cross_file);
   
   }
    
    if(photo_cross_database == "VULCAN_Diss"){
     MatrixXd photoXsection_info = ReadVULCANPhotoDissCrossSection(photo_cross_file);
   
   }
   
    if(photo_cross_database == "VPL"){
     MatrixXd photoXsection_info = ReadAtmosCrossSection(photo_cross_file);
   
   }
   
   
    }
    
//If the string is not equal to nan and then proceed and load the photochemical cross sectikons
   }
   
  }


//Atmosphere Properties, indices for I/O storage
  fstream InFile;
  int nSize = 0;
  string data1, data2, data3;
  InFile.open(atm_file); 
  getline(InFile, data1);
  getline(InFile, data1);
  while (getline(InFile, data1))
    nSize++;
  InFile.close();
  
//Setting initial condition for chemical species
  MatrixXd ChemMoleFrac(nsp, nSize);
  MatrixXd a(nsp, nSize);
  for (int i = 0; i < nSize; i++) {
      ChemMoleFrac.col(i) = mole_fractions;
    }
  

//Atmospheric Profile Data
  MatrixXd AtmData(4, nSize);
  int iTemp = 0;
  int iPress =  1;
  int iKzz =  2;
  int iAlt = 3;
  
//Input from txt file
  int inx = 0;
  InFile.open(atm_file); 
  getline(InFile, data1);
  getline(InFile, data1);
  while (InFile >> data1 >> data2 >> data3){
      AtmData(iPress, inx) = atof(data1.c_str())*0.1; //Pressure (dyn/cm^2) -> Pa
      AtmData(iTemp,inx) = atof(data2.c_str()); //Temperature (K)
      AtmData(iKzz, inx) = atof(data3.c_str())*1E-4; //Kzz (cm^2/s) -> m^2/s
      inx++;
      }
  InFile.close();
 
//Chemical evolution
double dh = 1000; //m
int iPrev, iNext;

//Initiating Matrices for Time Evolution
Eigen::SparseMatrix<double>  m_wjac;
m_wjac.resize(nsp, nsp);
MatrixXd mat1(nsp, nsp);
mat1 = MatrixXd::Identity(nsp, nsp);
MatrixXd mat2(nsp, nsp);
VectorXd Un(nsp);
VectorXd Unext(nsp);
VectorXd Uprev(nsp);
VectorXd flux1(nsp);
VectorXd flux2(nsp);
VectorXd flux3(nsp);
VectorXd flux4(nsp);
VectorXd dQ(nsp);
double Keddy;
double Temp, Press, Kzz;
VectorXd m_wdot(nsp);
VectorXd mole_frac(nsp);


for (int i = 0; i < nTime; i++) {
  for (int j = 0; j < nSize; j++) {
//Setting T, P, X for each grid point
    iPrev = j-1;
    iNext = j+1;
    Temp = AtmData(iTemp,j);
    Press = AtmData(iPress,j); 
    gas->setState_TP(Temp, (Press/1.0132E5)*OneAtm);
    mole_frac = ChemMoleFrac.col(j);
    gas->setMoleFractions_NoNorm(&mole_frac[0]);
    Keddy = AtmData(iKzz, j);
//Solving the net production for each species
    gas_kin->getNetProductionRates(&m_wdot[0]); //Extracting net production rates from Cantera
    m_wjac = gas_kin->netProductionRates_ddX(); //Extracting Jacobian from Cantera


//Backward Euler Scheme (Li and Chen, 2019)
    mat2 = ((mat1/dt) - m_wjac);
    mat2 = mat2.inverse();
    mat2 = mat2*m_wdot;
    a = ChemMoleFrac;
    if (j == 0){
    flux1 = flux_upper;
    }
    if (j == nSize-1){
    flux1 = flux_lower;
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
    ChemMoleFrac.col(j) = dQ + ChemMoleFrac.col(j) - (flux1*dt);
  
} 
    std::cout << "Simulation completed at time step t = " << i << std::endl;
} 

  delete pinput;
}
