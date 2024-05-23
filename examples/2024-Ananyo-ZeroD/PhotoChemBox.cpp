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

// application
#include <application/application.hpp>
//#include <interpolation.hpp>

// NetCDF Output
#if NETCDFOUTPUT
#include <netcdf.h>
#endif

// YAML file IO
#include "cantera/base/ct_defs.h"
#include "cantera/ext/yaml-cpp/yaml.h"


int main(int argc, char** argv) {
//Reading the input file, and creating Cantera objects
   std::string FileName = "SOx_box.yaml";
   auto sol = Cantera::newSolution(FileName);
   auto gas = sol->thermo();
   auto gas_kin = sol->kinetics();
   int nsp = gas->nSpecies();
   int nrxn = gas_kin->nReactions();

//Using YAML to extract problem parameters
   YAML::Node config  = YAML::LoadFile("SOx_box.yaml");
   YAML::Node pressNode = config["problem"][0]["pres"];
   double Pres = atof(pressNode.as<std::string>().c_str());
   if (pressNode) {
            std::cout << "Found press: " << pressNode.as<std::string>() << std::endl;
        } else {
            std::cout << "press not found." << std::endl;
        }

   YAML::Node tempNode = config["problem"][1]["temp"];
   double Temp = atof(tempNode.as<std::string>().c_str());
   if (tempNode) {
            std::cout << "Found temp: " << tempNode.as<std::string>() << std::endl;
        } else {
            std::cout << "temp not found." << std::endl;
	} 
 
//Creating solver objects, and matrices
   VectorXd mole_fractions = VectorXd::Zero(nsp);
   VectorXd Output = VectorXd::Zero(nsp);
   YAML::Node species_list = config["init"][0]["species"];
   std::string init_species_list = species_list.as<std::string>();
   std::regex pattern ("[a-zA-Z0-9_]+");
   std::smatch m;
   int species_inx;
   int inx = 1;
   while (std::regex_search (init_species_list,m,pattern)) {
      for (auto x:m){
        YAML::Node list_sp = config["init"][inx][x.str()];
        std::string species_init_condition = list_sp.as<std::string>();
        species_inx = gas->speciesIndex(x);
        mole_fractions(species_inx) = atof(species_init_condition.c_str());
    }
    inx++;
    init_species_list = m.suffix().str();
  }
   
   double h = 6.626e-34; //Planck's constant
   double c = 3e8;  //Speed of light 
   auto app = Application::GetInstance();
   auto stellar_input_file = app->FindResource("stellar/sun.ir");
   auto stellar_input = ReadStellarRadiationInput(stellar_input_file, 1., 1.);
   std::cout << "Radiation Input Complete!" << std::endl;

   double factor = 1/(h*c);
   Eigen::VectorXd wav = Eigen::Map<Eigen::VectorXd>(stellar_input.first.data(), stellar_input.first.size());;
  Eigen::VectorXd irr = Eigen::Map<Eigen::VectorXd>(stellar_input.second.data(), stellar_input.second.size());;
  Eigen::VectorXd actinicFlux = (wav.array()*irr.array()).matrix()*factor;

   std::shared_ptr<ActinicFlux> aflux = std::make_shared<ActinicFlux>();
  aflux->setWavelength(stellar_input.first);
  std::vector<double> actinicFluxVec(actinicFlux.data(), actinicFlux.data() + actinicFlux.size());

// Call the setTOAFlux function to set the actinic flux
   aflux->setTOAFlux(actinicFluxVec);
   aflux->initialize();
   sol->kinetics()->handleActinicFlux(aflux);

   YAML::Node time_max = config["integrator"][1]["Tmax"];
   double Tmax = atof(time_max.as<std::string>().c_str());
   gas->setMoleFractions(&mole_fractions[0]);
   gas->setState_TP(Temp, (Pres)*Cantera::OneBar);
   Cantera::ReactorNet sim;
   Cantera::IdealGasConstPressureReactor reactor;
   reactor.setSolution(sol);
   reactor.setEnergy(false);
   //reactor.initialize(0.0);
   double rtol = 1.0e-15;  // relative tolerance
   double atol = 1.0e-30; // absolute tolerance
   sim.addReactor(reactor);
   sim.setTolerances(rtol, atol);
   sim.setSensitivityTolerances(rtol, atol);
   reactor.initialize();
   sim.initialize();
   sim.advance(Tmax);
   gas->getMoleFractions(&Output[0]);
   std::cout << "<<<<<<<<< Mole Fraction >>>>>>>>>>" << std::endl;
   std::cout << Output.transpose() << std::endl;
   std::cout << "Wrapped up simulation" << std::endl;
}
