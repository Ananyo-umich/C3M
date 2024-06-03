// @sec3{Include files}
// Author: Ananyo Bhattacharya
// Affiliation: University of Michigan
// Email: ananyo@umich.edu
// The code computes the chemical evolution of a reaction network at a point

// Solution class describes a phase consists of a mixture of chemical species
#include <cantera/base/Solution.h>

// ThermoPhase object stores the thermodynamic state
#include <cantera/thermo.h>
#include <cantera/thermo/Phase.h>
// Kinetics object stores the chemical kinetics information
#include <cantera/kinetics.h>
#include <cantera/kinetics/Reaction.h>

// output stream
#include <cantera/numerics/Integrator.h>
#include <cantera/numerics/eigen_dense.h>
#include <cantera/numerics/eigen_sparse.h>
#include <cantera/zeroD/IdealGasConstPressureReactor.h>
#include <cantera/zeroD/Reactor.h>
#include <cantera/zeroD/ReactorNet.h>

#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

// C3M header
#include <c3m/PhotoChemistry.hpp>
#include <c3m/RadTran.hpp>
#include <c3m/actinic_flux.hpp>
#include <c3m/atm_chemistry.hpp>
#include <c3m/atm_chemistry_simulator.hpp>
#include <c3m/boundary.hpp>


// application
#include <application/application.hpp>
//#include <gtest/gtest.h>

// NetCDF Output
#if NETCDFOUTPUT
#include <netcdf.h>
#endif

// YAML file IO
#include "cantera/base/ct_defs.h"
#include "cantera/ext/yaml-cpp/yaml.h"


int main(int argc, char** argv) {
//Reading the input file, and creating Cantera objects
   Application::Start(argc, argv);
   std::string FileName = "SOx_box.yaml";
   auto sol = Cantera::newSolution(FileName);
   auto gas = sol->thermo();
   auto gas_kin = sol->kinetics();
   auto atm = std::make_shared<AtmChemistry>("atm", sol);
   int nsp = gas->nSpecies();
   int nrxn = gas_kin->nReactions();
   


//Using YAML to extract components of atmospheric structure
  YAML::Node config  = YAML::LoadFile(FileName);
  YAML::Node AtmFile = config["problem"][0]["planet"]; 

  int nz = 0;
  string data1, data2, data3, data4, data5;
  fstream InFile;
  InFile.open(FileName);
  getline(InFile, data1);
  getline(InFile, data1);
  while (getline(InFile, data1))
    nz++;
  InFile.close();

  MatrixXd AtmData(5, nz);
  int iTemp = 0;
  int iPress =  1;
  int iKzz =  2;
  int iAlt = 3;
  int iNd = 4;

  int inx = 0;
  InFile.open(FileName);
  getline(InFile, data1);
  getline(InFile, data1);
  while (InFile >> data1 >> data2 >> data3 >> data4 >> data5){
      AtmData(iPress, inx) = atof(data1.c_str())*1E2; //Pressure (mbar) -> Pa
      AtmData(iTemp,inx) = atof(data2.c_str()); //Temperature (K)
      AtmData(iKzz, inx) = atof(data3.c_str())*1E-4; //Kzz (cm^2/s) -> m^2/s
      AtmData(iAlt, inx) = atof(data4.c_str())*1E3; //Altitude (m) 
      AtmData(iNd, inx) = atof(data5.c_str())*1E6/(6.022E23*1E3); //Concentration (1/cm^3) -> kmol/m^3
      inx++;
      }
  InFile.close();

//Set up the grid
  double U0 = 0.;

// surface (lower boundary)
  auto surface = std::make_shared<SurfaceBoundary>("surface", sol);

// space  (upper boundary)
  auto space = std::make_shared<SpaceBoundary>("space", sol);


// set up simulation
  AtmChemistrySimulator *pchem;
  pchem = new AtmChemistrySimulator({surface, atm, space});

// resize function is called inside setupGrid
  atm->setupGrid(nz, AtmData.row(iAlt).data());
  std::shared_ptr<std::vector<double>> hydro;
  hydro = std::make_shared<std::vector<double>>(atm->nPoints() * 3);
  for (size_t n = 0; n < atm->nPoints(); ++n) {
      hydro->at(n * 3) = AtmData(iTemp,n);
      hydro->at(n * 3 + 1) = AtmData(iPress, n);
      hydro->at(n * 3 + 2) = U0;
    }

  atm->setHydro(hydro, 3);

//Homogeneous input condition
  YAML::Node species_list = config["init"][0]["species"];
   std::string init_species_list = species_list.as<std::string>();
   std::regex pattern ("[a-zA-Z0-9_]+");
   std::smatch m;
   int species_inx;
   inx = 1;
   while (std::regex_search (init_species_list,m,pattern)) {
      for (auto x:m){
        YAML::Node list_sp = config["init"][inx][x.str()];
        std::string species_init_condition = list_sp.as<std::string>();
        species_inx = gas->speciesIndex(x);
//Setting initial condition using flat profile function
        pchem->setFlatProfile(atm, species_inx, atof(species_init_condition.c_str()));
    }
    inx++;
    init_species_list = m.suffix().str();
  }
//Initial condition based on input profile (TBA)


//Integrator
  YAML::Node T_max = config["integrator"][1]["Tmax"];
  double Tmax = atof(T_max.as<std::string>().c_str());

//Extract the TOA stellar radiation, and scale it with orbit radius


//Boundary conditions are added here
  YAML::Node Surface_BC = config["problem"][1]["Tmax"]; 
  std::string SurfaceBC = Surface_BC.as<std::string>();

  YAML::Node Space_BC = config["problem"][1]["Tmax"];
  std::string init_species_list = species_list.as<std::string>();



  delete pchem;
  Application::Destroy();
}











