// external
#include <gtest/gtest.h>

// application
#include <application/application.hpp>

// Solution class describes a phase consists of a mixture of chemical species
#include <cantera/base/Solution.h>

// ThermoPhase object stores the thermodynamic state
#include <cantera/thermo.h>

// Kinetics object stores the chemical kinetics information
#include <cantera/kinetics.h>
#include <cantera/kinetics/Reaction.h>

// Reactor class, and numerics are responsible for solving the system
#include <cantera/zeroD/Reactor.h>
#include <cantera/zeroD/IdealGasConstPressureReactor.h>
#include <cantera/zeroD/ReactorNet.h>
#include <cantera/numerics/Integrator.h>
//#include <cantera/numerics/CVodesIntegrator.h>

// ZeroDim object stores the reactor information
#include <cantera/zerodim.h>

// c3m
#include <c3m/RadTran.hpp>

// YAML file IO
#include "cantera/base/ct_defs.h"
#include "cantera/ext/yaml-cpp/yaml.h"

#include <sstream>
#include <regex>

TEST(Input, ChapmanCycle) {

// Reading the chemical kinetics network
  auto sol = Cantera::newSolution("chapman_box.yaml");
  auto gas = sol->thermo();
  auto gas_kin = sol->kinetics();
  int nsp = gas->nSpecies();
  int nrxn = gas_kin->nReactions();
  YAML::Node config  = YAML::LoadFile("chapman_box.yaml");
// Reading the temperature, and pressure for box input file
// Find the node corresponding to "press"
  YAML::Node pressNode = config["problem"][0]["pres"];
  double Pres = atof(pressNode.as<std::string>().c_str());
  if (pressNode) {
            std::cout << "Found press: " << pressNode.as<std::string>() << std::endl;
        } else {
            std::cout << "press not found." << std::endl;
        }
  ASSERT_EQ(Pres, 2.5E-3);
// Find the node corresponding to "temp"
  YAML::Node tempNode = config["problem"][1]["temp"];
  double Temp = atof(tempNode.as<std::string>().c_str());       
  if (tempNode) {
            std::cout << "Found temp: " << tempNode.as<std::string>() << std::endl;
        } else {
            std::cout << "temp not found." << std::endl;
        }
  ASSERT_EQ(Temp, 251);
// Reading the time scale for integration
  
// Initiating solver objects and matrices
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
  
  //sim.initialize(); 
  double Tmax = 1E7;
  
// Test 1: Earth's Ozone chemistry at 30 km
  double P1 = 12E-3;
  double T1 = 227;
  double jO2_c1 = 6e-11;
  double jO3_c1 = 1.2e-3;

// Setting the multiplier for photodissociation reactions
  gas->setMoleFractions(&mole_fractions[0]);
  gas->setState_TP(T1, (P1)*Cantera::OneBar);
  gas_kin->setMultiplier(0, jO2_c1);
  gas_kin->setMultiplier(2, jO3_c1);
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
  int nmax = 1E9;
  Cantera::Integrator& integrator = sim.integrator();
  integrator.setMaxSteps(nmax);
  reactor.initialize();
  sim.initialize();
  sim.advance(Tmax);
  gas->getMoleFractions(&Output[0]);
  VectorXd krate(nrxn);
  gas_kin->getFwdRateConstants(&krate[0]);
  std::cout << krate.transpose() << std::endl; 

  std::cout << mole_fractions.transpose() << std::endl;
  std::cout << "Results for test at 30 km" << std::endl;
  std::cout << Output.transpose() << std::endl;
// Check against analytical solution
// ASSERT_NEAR(kfwd[0], 0.00127603, 1.0e-6);  
  ASSERT_NEAR(Output(0), 0.78, 0.01);
  ASSERT_NEAR(Output(1), 1.101859475627416e-09, 1E-12);
  ASSERT_NEAR(Output(2), 0.22, 0.01);
  ASSERT_NEAR(Output(3), 3.413346708728647e-05, 1E-8);
 
// Test 2: Earth's Ozone chemistry at 40 km
  double P2 = 2.5E-3;
  double T2 = 251;
  double jO2_c2 = 5e-10;
  double jO3_c2 = 1.9e-3;

// Setting the multiplier for photodissociation reactions
  gas->setMoleFractions(&mole_fractions[0]);
  gas->setState_TP(T2, (P2)*Cantera::OneBar);
  gas_kin->setMultiplier(0, jO2_c2);
  gas_kin->setMultiplier(2, jO3_c2);
  Cantera::ReactorNet sim2;
  Cantera::IdealGasConstPressureReactor reactor2;
  reactor2.setSolution(sol);
  reactor2.setEnergy(false);
  sim2.addReactor(reactor2);
  sim2.setTolerances(rtol, atol);
  sim2.setSensitivityTolerances(rtol, atol);
  Cantera::Integrator& integrator2 = sim2.integrator();
  integrator2.setMaxSteps(nmax);
  reactor2.initialize();
  sim2.initialize();
  sim2.advance(Tmax);
  gas->getMoleFractions(&Output[0]);
  gas_kin->getFwdRateConstants(&krate[0]);
  std::cout << krate.transpose() << std::endl;

  std::cout << "Results for test at 40 km" << std::endl;
  std::cout << Output.transpose() << std::endl;
// Check against analytical solution
  ASSERT_NEAR(Output(0), 0.78, 0.01);
  ASSERT_NEAR(Output(1), 3.5730286419165366e-08, 1E-9);
  ASSERT_NEAR(Output(2), 0.22, 0.01);
  ASSERT_NEAR(Output(3), 1.9446510471788222e-05, 1E-7);



}


int main(int argc, char **argv) {
  Application::Start(argc, argv);

  testing::InitGoogleTest(&argc, argv);

  int result = RUN_ALL_TESTS();

  Application::Destroy();

  return result;
}
