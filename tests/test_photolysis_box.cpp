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

TEST(ZeroDim, OxygenPhotolysisBox) {
  auto app = Application::GetInstance();

  // Reading the chemical kinetics network
  auto sol = Cantera::newSolution("test_photolysis_box.yaml");

  // Initial condition for mole fraction
  sol->thermo()->setState_TPX(250., 0.1 * Cantera::OneAtm, "O2:0.21, N2:0.78");

  // Calculating photochemical reaction rate
  auto stellar_input_file = app->FindResource("stellar/sun.ir");
  auto stellar_input = ReadStellarRadiationInput(stellar_input_file, 1., 1.);
  std::cout << "Radiation Input Complete!" << std::endl;

  // Updating the actinic flux within yaml file [All in SI units]
  sol->kinetics()->setWavelength(stellar_input.row(0).data(),
                                 stellar_input.row(0).size());
  sol->kinetics()->updateActinicFlux(stellar_input.row(1).data());

  // Reactor
  Cantera::Reactor reactor(sol);
  reactor.initialize();

  // Reactor Network
  Cantera::ReactorNet network;
  network.addReactor(reactor);
  network.initialize();

  double time_step = 1.e-5;
  double max_time = 1.e-2;

  double time = 0.;
  while (network.time() < max_time) {
    time = network.time() + time_step;
    network.advance(time);
  }
}

int main(int argc, char **argv) {
  Application::Start(argc, argv);

  testing::InitGoogleTest(&argc, argv);

  int result = RUN_ALL_TESTS();

  Application::Destroy();

  return result;
}
