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

// ZeroDim object stores the reactor information
#include <cantera/zerodim.h>

// c3m
#include <c3m/RadTran.hpp>

// YAML file IO
#include "cantera/base/ct_defs.h"
#include "cantera/ext/yaml-cpp/yaml.h"

TEST(Input, ChapmanCycle) {

// Reading the chemical kinetics network
  auto sol = Cantera::newSolution("chapman_box.yaml");
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
  

// Test 1: Earth's Ozone chemistry at 30 km
  P1 = 2.5E-3
  T1 = 251


// Test 2: Earth's Ozone chemistry at 40 km
  P2 = 2
  T2 = 2 
  
}


int main(int argc, char **argv) {
  Application::Start(argc, argv);

  testing::InitGoogleTest(&argc, argv);

  int result = RUN_ALL_TESTS();

  Application::Destroy();

  return result;
}
