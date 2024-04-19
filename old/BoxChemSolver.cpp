// @sec3{Include files}
// Author: Ananyo Bhattacharya
// Affiliation: University of Michigan
// Email: ananyo@umich.edu
// The code computes the chemical evolution of a reaction network at a point


// Solution class describes a phase consists of a mixture of chemical species
#include <cantera/base/Solution.h>

// ThermoPhase object stores the thermodynamic state
#include <cantera/thermo.h>

// Kinetics object stores the chemical kinetics information
#include <cantera/kinetics.h>
#include <cantera/kinetics/Reaction.h>

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
  infile.Open("test_athena.inp", IOWrapper::FileMode::read);
  ParameterInput *pinput = new ParameterInput();

  pinput->LoadFromFile(infile);
  infile.Close();

//Loading the input parameters for atmospheric profile and reaction network
  std::string network_file = pinput->GetString("problem", "network");
  
//Reading the chemical kinetics network
  auto sol = newSolution(network_file);
  //auto sol = newSolution(network_file);
  auto gas = sol->thermo();
  auto gas_kin = sol->kinetics();  
  int nsp = gas->nSpecies();
  int nrxn = gas_kin->nReactions();
  std::cout << "Network imported" << std::endl;
//Initial condition for mole fraction
  VectorXd mole_fractions = VectorXd::Zero(nsp);
  

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
 std::cout << "Initial conditions loaded" << std::endl;
 
//Loading the atmospheric conditions
  std::string tp = pinput->GetString("problem", "temp");
  std::string ps = pinput->GetString("problem", "pres");
  double temp = atof(tp.c_str());
  double pres = atof(ps.c_str());

//Integrator
  std::string time_step = pinput->GetString("integrator", "dt");
  std::string max_time = pinput->GetString("integrator", "Tmax");
  double dt = atof(time_step.c_str());
  double Tmax = atof(max_time.c_str());
  double Ttot = 0.0;
  std::cout << "Time step: " << dt << std::endl;
  //int nTime = stoi(total_time_steps.c_str());


//Atmosphere Properties, indices for I/O storage

  

std::cout << "All inputs loaded into C3M " << std::endl;

//Photochemistry
  std::string photochem = pinput->GetString("problem", "photochem");
  if(photochem == "true")
  {
  
  std::cout << "Starting photochemistry calculations" << std::endl;
//Radiative transfer
//Reading the Stellar Irradiance Input
  std::string stellar_input_file = pinput->GetString("radtran", "solar");
  std::string radius = pinput->GetString("radtran", "radius");
  std::string reference =   pinput->GetString("radtran", "reference");
  
  double rad = atof(radius.c_str());
  double ref = atof(reference.c_str());
  std::cout << rad << " " << ref << std::endl;
  MatrixXd stellar_input = ReadStellarRadiationInput(stellar_input_file, rad, ref);
  std::cout << "Radiation Input Complete!" << std::endl;

  std::cout << "Initiating photochemistry!" << std::endl;
//Storing the photochemical cross section data
  int PhotoRxn = 0;
  
/*
 *auto& react = *(NetworkName->reaction(inumRxn));
  std::string Equation = react.equation();
 */
  
  for(int irxn = 0; irxn < nrxn; irxn++){
    auto& rxnObj = *(gas_kin->reaction(irxn));
    std::string rxnEquation = rxnObj.equation();
    int pos = rxnEquation.find("=");
    rxnEquation.replace(pos, 2, "->");
    std::string photo_cross_info = pinput->GetOrAddString("photocross", rxnEquation, "nan");
    if(photo_cross_info != "nan"){
    PhotoRxn++;
    }
    
    }

  Eigen::MatrixXd photo_cross_data(stellar_input.row(0).size(), PhotoRxn);
  Eigen::VectorXd RxnIndex(PhotoRxn);
  Eigen::VectorXd Jrate(PhotoRxn);
  
//Reading cross section for all photochemical reactions
  int ph_inx = 0;
  for(int irxn = 0; irxn < nrxn; irxn++){
    auto& rxnObj = *(gas_kin->reaction(irxn));
    std::string rxnEquation = rxnObj.equation();
    std::cout << rxnEquation << std::endl;
    int pos = rxnEquation.find("=");
    rxnEquation.replace(pos, 2, "->");
    std::string photo_cross_info = pinput->GetOrAddString("photocross", rxnEquation, "nan");
    std::cout << photo_cross_info << std::endl;    
    if(photo_cross_info != "nan"){
    RxnIndex(ph_inx) = irxn;
    std::string photo_cross_database = photo_cross_info.substr(0,photo_cross_info.find(","));  
    std::string photo_cross_file = photo_cross_info.substr(photo_cross_info.find(",")+1,photo_cross_info.length()-1);
    std::cout << "Identifying the reactions" << std::endl;
    
    if(photo_cross_database == "VULCAN_Ion"){
     MatrixXd photoXsection_info = ReadVULCANPhotoIonCrossSection(photo_cross_file);
     photo_cross_data.col(ph_inx) = InterpolateCrossSection(stellar_input.row(0), photoXsection_info.row(0), photoXsection_info.row(1)); 
   }
    
    if(photo_cross_database == "VULCAN_Diss"){
     MatrixXd photoXsection_info = ReadVULCANPhotoDissCrossSection(photo_cross_file);
     photo_cross_data.col(ph_inx) = InterpolateCrossSection(stellar_input.row(0), photoXsection_info.row(0), photoXsection_info.row(1));
   
   }
   
    if(photo_cross_database == "VPL"){
     MatrixXd photoXsection_info = ReadAtmosCrossSection(photo_cross_file);
     photo_cross_data.col(ph_inx) = InterpolateCrossSection(stellar_input.row(0), photoXsection_info.row(0), photoXsection_info.row(1));
   
   }

   if(photo_cross_database == "MPD"){
     MatrixXd photoXsection_info = ReadMPCrossSection(photo_cross_file);
     photo_cross_data.col(ph_inx) = InterpolateCrossSection(stellar_input.row(0), photoXsection_info.row(0), photoXsection_info.row(1));

   }
   
   ph_inx++;
    }
    
//If the string is not equal to nan and then proceed and load the photochemical cross sectikons
   }

//Chemical evolution
 std::cout << "Starting chemical evolution " << std::endl;
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

double Kb = 1.38E-23;
int i = 0;
double j_rate;
double j_O3;
VectorXd krate(nrxn);
//Calculating the photochemical reaction rate
std::cout << "Computing the photochemical reaction rates" << std::endl;
//Computing the photochemical reaction rate for each reaction
for(int rx = 0; rx < PhotoRxn; rx++){
  j_rate = PhotoChemRate(stellar_input.row(0),photo_cross_data.col(rx), stellar_input.row(1).transpose());
  std::cout << j_rate << std::endl;
//Setting the multiplier for each reaction
 gas_kin->setMultiplier(RxnIndex(rx), j_rate); 
 }

std::string ionkin_state = pinput->GetString("problem", "ionkinetics");
if(ionkin_state == "true"){
for(int irxn = 0; irxn < nrxn; irxn++){
    auto& rxnObj = *(gas_kin->reaction(irxn));
    std::string rxnEqn = rxnObj.equation();
    int pos = rxnEqn.find("=");
    rxnEqn.replace(pos, 2, "->");
    std::string ion_rate = pinput->GetOrAddString("ionkinetics_man", rxnEqn, "nan");
    if(ion_rate != "nan"){
    gas_kin->setMultiplier(irxn, atof(ion_rate.c_str())); 
    }
    
    }
    
}
std::cout << "Photochemical reaction rates computed!" << std::endl;
std::string OutFileName = pinput->GetString("output", "file");
std::ofstream outfile (OutFileName);

std::string time_step = pinput->GetOrAddString("problem", "time", "nan");
while(Ttot  < Tmax) {
//Switch off the precipitating source after certain time step
   if(ionkin_state == "true"){
      std::string source_time = pinput->GetString("integrator", "Tsource");
      double t_source = atof(source_time.c_str());
      if(Ttot >= t_source){
std::cout << "Source is switched off" << std::endl;
       for(int irxn = 0; irxn < nrxn; irxn++){
         auto& rxnObj = *(gas_kin->reaction(irxn));
         std::string rxnEqn = rxnObj.equation();
         int pos = rxnEqn.find("=");
         rxnEqn.replace(pos, 2, "->");
         std::string ion_rate = pinput->GetOrAddString("ionkinetics_man", rxnEqn, "nan");
         if(ion_rate != "nan"){
            gas_kin->setMultiplier(irxn, 0.0); 
            }
    
    }
      
      }
      
      
    }
//Setting T, P, X for each grid point
    gas->setState_TP(temp, (pres)*OneBar);
    gas->setMoleFractions(&mole_fractions[0]);
//Solving the net production for each species
    gas_kin->getNetProductionRates(&m_wdot[0]); //Extracting net production rates from Cantera
    m_wjac = gas_kin->netProductionRates_ddX(); //Extracting Jacobian from Cantera
    m_wjac = m_wjac/ gas->molarDensity(); 
//Integration for each species
//Backward Euler Scheme (Li and Chen, 2019)
    mat2 = ((mat1/dt) - m_wjac);
    mat2 = mat2.inverse();
    mat2 = mat2*(m_wdot / gas->molarDensity());
    dQ = mat2;
    mole_fractions = mole_fractions + dQ;
    //gas->setMoleFractions(&mole_fractions[0]);
    //gas->getMoleFractions(&mole_fractions[0]);   
    //std::cout << "Simulation completed at time t = " << Ttot << std::endl;
    Ttot = Ttot + dt;
    if(time_step == "log"){
      dt = dt*1.25;
    }
    gas_kin->getFwdRatesOfProgress(&krate[0]);
    std::cout << m_wdot.transpose() << std::endl;
    outfile << Ttot << " " << mole_fractions.transpose()  << std::endl;
} 
    std::cout << gas->report() << std::endl;
    //std::cout << j_rate << std::endl;
    outfile.close();
//Analytical solution
//double ratio = j_O3/(krate(1)*0.21*number_density);
  }

//This is where the photochemistry definition ends (Anything beyond is not defined in the scope of photochemistry) 


  if(photochem != "true"){
//If photochemistry is set equal to false
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
double dh;
int iPrev, iNext;

double Kb = 1.38E-23;
std::string ionkin_state = pinput->GetString("problem", "ionkinetics");
if(ionkin_state == "true"){
for(int irxn = 0; irxn < nrxn; irxn++){
    auto& rxnObj = *(gas_kin->reaction(irxn));
    std::string rxnEqn = rxnObj.equation();
    int pos = rxnEqn.find("=");
    rxnEqn.replace(pos, 2, "->");
    std::string ion_rate = pinput->GetOrAddString("ionkinetics_man", rxnEqn, "nan");
    if(ion_rate != "nan"){
    gas_kin->setMultiplier(irxn, atof(ion_rate.c_str())); 
    }
    
    }}


int i = 0;
std::string OutFileName = pinput->GetString("output", "file");
std::ofstream outfile (OutFileName);
std::string time_step = pinput->GetOrAddString("problem", "time", "nan");
while(Ttot  < Tmax) {
//Switch off the precipitating source after certain time step
   if(ionkin_state == "true"){
      std::string source_time = pinput->GetString("integrator", "Tsource");
      double t_source = atof(source_time.c_str());
      if(Ttot >= t_source){
std::cout << "Source is switched off" << std::endl;
       for(int irxn = 0; irxn < nrxn; irxn++){
         auto& rxnObj = *(gas_kin->reaction(irxn));
         std::string rxnEqn = rxnObj.equation();
         int pos = rxnEqn.find("=");
         rxnEqn.replace(pos, 2, "->");
         std::string ion_rate = pinput->GetOrAddString("ionkinetics_man", rxnEqn, "nan");
         if(ion_rate != "nan"){
            gas_kin->setMultiplier(irxn, 0.0); 
            }
    
    }
      
      }
      
      
    }

//Setting T, P, X for each grid point
    gas->setState_TP(temp, (pres)*OneBar);
    gas->setMoleFractions(&mole_fractions[0]);
//Solving the net production for each species
    gas_kin->getNetProductionRates(&m_wdot[0]); //Extracting net production rates from Cantera
    gas_kin->getFwdRatesOfProgress(&m_wdot[0]);
    std::cout << 1E3*6.022E23*1E-6*m_wdot.transpose() << std::endl;
    m_wjac = gas_kin->netProductionRates_ddX(); //Extracting Jacobian from Cantera
    m_wjac = m_wjac/ gas->molarDensity(); 
//Integration for each species
//Backward Euler Scheme (Li and Chen, 2019)
    mat2 = ((mat1/dt) - m_wjac);
    mat2 = mat2.inverse();
    mat2 = mat2*(m_wdot / gas->molarDensity());
    dQ = mat2;
    mole_fractions = mole_fractions + dQ;
    //gas->setMoleFractions(&mole_fractions[0]);
    //gas->getMoleFractions(&mole_fractions[0]);   
    std::cout << "Simulation completed at time t = " << Ttot << std::endl;
    Ttot = Ttot + dt;
    if(time_step == "log"){
      dt = dt*1.25;
    }
    outfile << Ttot << " " << mole_fractions.transpose()  << std::endl;
}

    std::cout << gas->report() << std::endl;
    outfile.close();
}
//Writing output into NetCDF file



delete pinput;
}
