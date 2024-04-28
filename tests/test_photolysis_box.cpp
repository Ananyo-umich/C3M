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

class PhotolysisCH4 : public testing::Test {
 public:
  // data
  shared_ptr<Cantera::ThermoPhase> phase;
  shared_ptr<Cantera::Kinetics> kin;

  // constructor
  PhotolysisCH4() {
    phase = Cantera::newThermo("photolysis_ch4.yaml");
    kin = Cantera::newKinetics({phase}, "photolysis_ch4.yaml");

    // set the initial state
    std::string X = "CH4:0.02 N2:0.98";
    phase->setState_TPX(200.0, Cantera::OneAtm, X);

    // set wavelength
    std::vector<double> wavelength(10);
    std::vector<double> actinic_flux(10);

    for (int i = 0; i < 10; i++) {
      wavelength[i] = 20.0e-9 + i * 20.0e-9;
      actinic_flux[i] = 1.e18;
    }

    kin->setWavelength(wavelength.data(), wavelength.size());
    kin->updateActinicFlux(actinic_flux.data());
  }
};

TEST_F(PhotolysisCH4, check_phase) {
  ASSERT_EQ(phase->nElements(), 3);
  ASSERT_EQ(phase->nSpecies(), 8);
}

TEST_F(PhotolysisCH4, check_kinetics) {
  ASSERT_EQ(kin->nReactions(), 2);
  ASSERT_EQ(kin->nTotalSpecies(), 8);
  ASSERT_EQ(kin->nPhases(), 1);
  ASSERT_EQ(kin->nWavelengths(), 10);
}

TEST_F(PhotolysisCH4, check_fwd_rate_constants) {
  std::vector<double> kfwd(kin->nReactions());

  kin->getFwdRateConstants(kfwd.data());

  ASSERT_NEAR(kfwd[0], 1.1178762941869657e-10, 1.0e-14);
  ASSERT_NEAR(kfwd[1], 0., 1.0e-18);

  int iCH4 = kin->kineticsSpeciesIndex("CH4");
  int iCH3 = kin->kineticsSpeciesIndex("CH3");
  int i1CH2 = kin->kineticsSpeciesIndex("(1)CH2");
  int i3CH2 = kin->kineticsSpeciesIndex("(3)CH2");
  int iCH = kin->kineticsSpeciesIndex("CH");
  int iH2 = kin->kineticsSpeciesIndex("H2");
  int iH = kin->kineticsSpeciesIndex("H");

  ASSERT_EQ(iCH4, 0);
  ASSERT_EQ(iCH3, 1);
  ASSERT_EQ(i1CH2, 2);
  ASSERT_EQ(i3CH2, 3);
  ASSERT_EQ(iCH, 4);
  ASSERT_EQ(iH2, 5);
  ASSERT_EQ(iH, 6);

  double kCH3 = kin->productStoichCoeff(iCH3, 0);
  ASSERT_NEAR(kCH3, 0.390204, 1.0e-4);

  double k1CH2 = kin->productStoichCoeff(i1CH2, 0);
  ASSERT_NEAR(k1CH2, 0.268438, 1.0e-4);

  double k3CH2 = kin->productStoichCoeff(i3CH2, 0);
  ASSERT_NEAR(k3CH2, 0.103706, 1.0e-4);

  double kCH = kin->productStoichCoeff(iCH, 0);
  ASSERT_NEAR(kCH, 0.237652, 1.0e-4);

  double kH2 = kin->productStoichCoeff(iH2, 0);
  ASSERT_NEAR(kH2, 0.50609, 1.0e-4);

  double kH = kin->productStoichCoeff(iH, 0);
  ASSERT_NEAR(kH, 0.835268, 1.0e-4);

  ASSERT_NEAR(kCH3 + k1CH2 + k3CH2 + kCH, 1.0, 1.0e-14);
  ASSERT_NEAR(3 * kCH3 + 2 * k1CH2 + 2 * k3CH2 + kCH + 2 * kH2 + kH, 4.0,
              1.0e-14);
}

/*TEST(ZeroDim, PhotolysisO2) {
  auto app = Application::GetInstance();

  // Reading the chemical kinetics network
  auto sol = Cantera::newSolution("photolysis_o2.yaml");

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
  Cantera::IdealGasReactor reactor(sol);
  reactor.setEnergy(false);
  reactor.initialize();

  std::cout << "T = " << reactor.temperature() << std::endl;
  std::cout << "rho = " << reactor.density() << std::endl;
  std::cout << "mass fractions = ";
  for (size_t i = 0; i < sol->thermo()->nSpecies(); i++) {
    std::cout << sol->thermo()->speciesName(i) << ":"
              << sol->thermo()->massFraction(i) << " ";
  }
  std::cout << std::endl;

  // Reactor Network
  Cantera::ReactorNet network;
  network.addReactor(reactor);
  network.initialize();

  double time_step = 10.;
  double max_time = 1000.;

  double time = 0.;
  while (network.time() < max_time) {
    time = network.time() + time_step;
    network.advance(time);
  }
  std::cout << "Time = " << time << std::endl;

  std::cout << "T = " << reactor.temperature() << std::endl;
  std::cout << "rho = " << reactor.density() << std::endl;
  std::cout << "mass fractions = ";
  for (size_t i = 0; i < sol->thermo()->nSpecies(); i++) {
    std::cout << sol->thermo()->speciesName(i) << ":"
              << sol->thermo()->massFraction(i) << " ";
  }
  std::cout << std::endl;
}*/

int main(int argc, char **argv) {
  Application::Start(argc, argv);

  testing::InitGoogleTest(&argc, argv);

  int result = RUN_ALL_TESTS();

  Application::Destroy();

  return result;
}
