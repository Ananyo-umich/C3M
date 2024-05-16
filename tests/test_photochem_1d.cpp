// external
#include <gtest/gtest.h>

// application
#include <application/application.hpp>

// cantera
#include <cantera/base/Solution.h>

// c3m
#include <c3m/atm_chemistry.hpp>
#include <c3m/atm_chemistry_simulator.hpp>
#include <c3m/boundary.hpp>

class TestPhotochem1D : public testing::Test {
 public:
  // data
  AtmChemistrySimulator *pchem;
  int iN2, iO, iO2, iO3;
  std::shared_ptr<std::vector<double>> hydro;

  // constructor
  TestPhotochem1D() {
    // atmosphere
    auto mech = Cantera::newSolution("photolysis_o2.yaml");
    auto atm = std::make_shared<AtmChemistry>("atm", mech);

    iN2 = atm->componentIndex("N2");
    iO = atm->componentIndex("O");
    iO2 = atm->componentIndex("O2");
    iO3 = atm->componentIndex("O3");

    // grid
    double T0 = 300;
    double P0 = 0.01;
    double U0 = 0.;

    int nz = 10;
    double height = 100.E3;
    std::vector<double> z(nz);
    double dz = height / nz;
    for (int iz = 0; iz < nz; iz++) {
      z[iz] = dz / 2. + iz * dz;
    }

    // resize function is called inside setupGrid
    atm->setupGrid(nz, z.data());

    hydro = std::make_shared<std::vector<double>>(atm->nPoints() * 3);
    for (size_t n = 0; n < atm->nPoints(); ++n) {
      hydro->at(n * 3) = T0;
      hydro->at(n * 3 + 1) = P0;
      hydro->at(n * 3 + 2) = U0;
    }

    atm->setHydro(hydro, 3);

    // surface
    auto surface = std::make_shared<SurfaceBoundary>("surface", mech);

    // space
    auto space = std::make_shared<SpaceBoundary>("space", mech);

    // set up simulation
    pchem = new AtmChemistrySimulator({surface, atm, space});
    // pchem = new AtmChemistrySimulator({atm});
    pchem->initFromFile("stellar/sun.ir");

    double vN2 = 0.70;
    double vO = 0.1;
    double vO2 = 0.21;
    double vO3 = 0.1;

    pchem->setFlatProfile(atm, iN2, vN2);
    pchem->setFlatProfile(atm, iO, vO);
    pchem->setFlatProfile(atm, iO2, vO2);
    pchem->setFlatProfile(atm, iO3, vO3);

    std::string X = "O2:0.21 N2:0.79";
    pchem->find<Connector>("surface")->setSpeciesDirichlet(X);
    pchem->find<Connector>("space")->setSpeciesDirichlet(X);
  }

  ~TestPhotochem1D() { delete pchem; }
};

TEST_F(TestPhotochem1D, check_domain_index) {
  ASSERT_EQ(iN2, 0);
  ASSERT_EQ(iO, 1);
  ASSERT_EQ(iO2, 3);
  ASSERT_EQ(iO3, 4);

  auto atm = pchem->find<>("atm");
  int ncomp = atm->nComponents();
  ASSERT_EQ(ncomp, 5);
  int npoints = atm->nPoints();
  ASSERT_EQ(npoints, 10);
  int size = pchem->size();
  ASSERT_EQ(size, (npoints + 2) * ncomp);
}

TEST_F(TestPhotochem1D, check_profile) {
  auto atm = pchem->find<AtmChemistry>("atm");

  double vN2 = 0.70;
  double vO = 0.1;
  double vO2 = 0.21;
  double vO3 = 0.1;

  for (int j = 0; j < atm->nPoints(); j++) {
    ASSERT_DOUBLE_EQ(pchem->value(atm, iN2, j), vN2);
    ASSERT_DOUBLE_EQ(pchem->value(atm, iO, j), vO);
    ASSERT_DOUBLE_EQ(pchem->value(atm, iO2, j), vO2);
    ASSERT_DOUBLE_EQ(pchem->value(atm, iO3, j), vO3);
  }
}

TEST_F(TestPhotochem1D, check_time_step) {
  auto atm = pchem->find<AtmChemistry>("atm");

  pchem->setMaxTimeStep(1.E9);
  pchem->show();

  int nsteps = 10;
  double dt = 1.0;

  pchem->timeStep(200, dt, 8);
  pchem->show();
}

int main(int argc, char **argv) {
  Application::Start(argc, argv);

  testing::InitGoogleTest(&argc, argv);

  int result = RUN_ALL_TESTS();

  Application::Destroy();

  return result;
}
