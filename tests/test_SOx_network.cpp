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

class SOxNet : public testing::Test {
 public:
  // data
  shared_ptr<Cantera::ThermoPhase> gas;
  shared_ptr<Cantera::Kinetics> gas_kin;

  // constructor
  SOxNet() {
    auto sol = Cantera::newSolution("SOx_box.yaml");
    gas = sol->thermo();
    gas_kin = sol->kinetics();
    int nsp = gas->nSpecies();
    int nrxn = gas_kin->nReactions();
    YAML::Node config  = YAML::LoadFile("SOx_box.yaml");
 
  // Find the node corresponding to "press"
    YAML::Node pressNode = config["problem"][0]["pres"];
    double Pres = atof(pressNode.as<std::string>().c_str());
    if (pressNode) {
            std::cout << "Found press: " << pressNode.as<std::string>() << std::endl;
        } else {
            std::cout << "press not found." << std::endl;
        }
  // Find the node corresponding to "temp"
    YAML::Node tempNode = config["problem"][1]["temp"];
    double Temp = atof(tempNode.as<std::string>().c_str());
    if (tempNode) {
            std::cout << "Found temp: " << tempNode.as<std::string>() << std::endl;
        } else {
            std::cout << "temp not found." << std::endl;
        }


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

  YAML::Node time_max = config["integrator"][1]["Tmax"];
  double Tmax = atof(time_max.as<std::string>().c_str());

 //Test SOx reactions for Venus atmosphere at 68 km. Photochemical reaction rates based on Zhang et al., (2012)
  double J_SO2 = 1.92E-4;
  double J_SO3 = 3.94E-5;

 // Setting the multiplier for photodissociation reactions
  gas->setMoleFractions(&mole_fractions[0]);
  gas->setState_TP(Temp, (Pres)*Cantera::OneBar);
  gas_kin->setMultiplier(0, J_SO2);
  gas_kin->setMultiplier(1, J_SO3);
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

  }
};







int main(int argc, char **argv) {
  Application::Start(argc, argv);

  testing::InitGoogleTest(&argc, argv);

  int result = RUN_ALL_TESTS();

  Application::Destroy();

  return result;
}

