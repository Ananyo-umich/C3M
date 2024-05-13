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


// Reading the time scale for integration


// Test 1: Earth's Ozone chemistry at 30 km



// Test 2: Earth's Ozone chemistry at 40 km

  
}


int main(int argc, char **argv) {
  Application::Start(argc, argv);

  testing::InitGoogleTest(&argc, argv);

  int result = RUN_ALL_TESTS();

  Application::Destroy();

  return result;
}
