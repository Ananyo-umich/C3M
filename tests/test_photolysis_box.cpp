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

class MethanePhotolysis : public testing::Test {
 public:
  // data
  shared_ptr<Cantera::ThermoPhase> phase;
  shared_ptr<Cantera::Kinetics> kin;

  // constructor
  MethanePhotolysis() {
    phase = Cantera::newThermo("../data/ch4_photolysis.yaml");
    kin = Cantera::newKinetics({phase}, "../data/ch4_photolysis.yaml");

    // set the initial state
    std::string X = "CH4:0.02 N2:0.98";
    phase->setState_TPX(200.0, Cantera::OneAtm, X);

    // set wavelength
    std::vector<double> wavelength(10);
    std::vector<double> actinic_flux(10);

    for (int i = 0; i < 10; i++) {
      wavelength[i] = 20.0 + i * 20.0;
      actinic_flux[i] = 1.0;
    }

    kin->setWavelength(wavelength.data(), wavelength.size());
    kin->updateActinicFlux(actinic_flux.data());
  }
};

TEST_F(MethanePhotolysis, check_phase) {
  ASSERT_EQ(phase->nElements(), 3);
  ASSERT_EQ(phase->nSpecies(), 8);
}

TEST_F(MethanePhotolysis, check_kinetics) {
  ASSERT_EQ(kin->nReactions(), 2);
  ASSERT_EQ(kin->nTotalSpecies(), 8);
  ASSERT_EQ(kin->nPhases(), 1);
  ASSERT_EQ(kin->nWavelengths(), 10);
}

TEST_F(MethanePhotolysis, check_fwd_rate_constants) {
  std::vector<double> kfwd(kin->nReactions());

  kin->getFwdRateConstants(kfwd.data());

  ASSERT_NEAR(kfwd[0], 3.06820e-14, 1.0e-18);
  ASSERT_NEAR(kfwd[1], 3.2e-16, 1.0e-18);

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

  double kCH4 = kin->productStoichCoeff(iCH4, 0);
  ASSERT_NEAR(kCH4, 0.635657, 1.0e-4);

  double kCH3 = kin->productStoichCoeff(iCH3, 0);
  ASSERT_NEAR(kCH3, 0.142168, 1.0e-4);

  double k1CH2 = kin->productStoichCoeff(i1CH2, 0);
  ASSERT_NEAR(k1CH2, 0.0978033, 1.0e-4);

  double k3CH2 = kin->productStoichCoeff(i3CH2, 0);
  ASSERT_NEAR(k3CH2, 0.0377844, 1.0e-4);

  double kCH = kin->productStoichCoeff(iCH, 0);
  ASSERT_NEAR(kCH, 0.0865869, 1.0e-4);

  double kH2 = kin->productStoichCoeff(iH2, 0);
  ASSERT_NEAR(kH2, 0.18439, 1.0e-4);

  double kH = kin->productStoichCoeff(iH, 0);
  ASSERT_NEAR(kH, 0.304324, 1.0e-4);

  ASSERT_NEAR(kCH4 + kCH3 + k1CH2 + k3CH2 + kCH, 1.0, 1.0e-14);
  ASSERT_NEAR(4 * kCH4 + 3 * kCH3 + 2 * k1CH2 + 2 * k3CH2 + kCH + 2 * kH2 + kH,
              4.0, 1.0e-14);
}

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
