// external
#include <gtest/gtest.h>

// application
#include <application/application.hpp>

// c3m
#include <c3m/atm_chemistry.hpp>
#include <c3m/atm_chemistry_simulator.hpp>
#include <c3m/boundary.hpp>

class TestAdvectionDiffusion : public testing::Test {
 public:
  // data
  AtmChemistrySimulator *pchem;
  std::shared_ptr<std::vector<double>> hydro;
  int iN2, iO;

  // constructor
  TestAdvectionDiffusion() {
    auto mech = Cantera::newSolution("advection_diffusion.yaml");
    auto atm = std::make_shared<AtmChemistry>("atm", mech);

    iN2 = atm->componentIndex("N2");
    iO = atm->componentIndex("O");

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
    pchem->initFromFile("stellar/sun.ir");

    double vN2 = 0.9;
    double vO = 0.1;
    pchem->setFlatProfile(atm, iN2, vN2);
    pchem->setFlatProfile(atm, iO, vO);

    std::string X = "O:0. N2:1.";
    pchem->find<Connector>("surface")->setSpeciesDirichlet(X);

    X = "O:0.9 N2:0.1";
    pchem->find<Connector>("space")->setSpeciesDirichlet(X);
  }

  ~TestAdvectionDiffusion() { delete pchem; }
};

TEST_F(TestAdvectionDiffusion, check_index) {
  ASSERT_EQ(iN2, 0);
  ASSERT_EQ(iO, 1);

  auto atm = pchem->find<>("atm");
  int ncomp = atm->nComponents();
  ASSERT_EQ(ncomp, 2);
  int npoints = atm->nPoints();
  ASSERT_EQ(npoints, 10);
  int size = pchem->size();
  ASSERT_EQ(size, (npoints + 2) * ncomp);
}

TEST_F(TestAdvectionDiffusion, check_diffusion) {
  auto atm = pchem->find<>("atm");

  int nstep = 499;
  double dt = 1.0;

  pchem->timeStep(nstep, dt, 8);
  pchem->show();
}

TEST_F(TestAdvectionDiffusion, check_neg_advection) {
  auto atm = pchem->find<>("atm");

  for (size_t n = 0; n < atm->nPoints(); ++n) {
    hydro->at(n * 3 + 2) = -10.;
  }

  int nstep = 100;
  double dt = 1.0;

  pchem->timeStep(nstep, dt, 8);
  pchem->show();
}

TEST_F(TestAdvectionDiffusion, check_pos_advection) {
  auto atm = pchem->find<>("atm");

  for (size_t n = 0; n < atm->nPoints(); ++n) {
    hydro->at(n * 3 + 2) = 10.;
  }

  int nstep = 100;
  double dt = 1.0;

  pchem->timeStep(nstep, dt, 8);
  pchem->show();
}

int main(int argc, char **argv) {
  Application::Start(argc, argv);

  testing::InitGoogleTest(&argc, argv);

  int result = RUN_ALL_TESTS();

  Application::Destroy();

  return result;
}
