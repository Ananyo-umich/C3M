// external
#include <gtest/gtest.h>

// application
#include <application/application.hpp>

// cantera
#include <cantera/base/Solution.h>
#include <cantera/oneD/DomainFactory.h>
#include <cantera/onedim.h>

// c3m
#include <c3m/atm_chemistry.hpp>
#include <c3m/atm_chemistry_simulator.hpp>
#include <c3m/boundary.hpp>

class TestPhotochem1D : public testing::Test {
 public:
  // data
  AtmChemistrySimulator *pchem;
  int iT, iP, iN2, iO, iO2, iO3;

  // constructor
  TestPhotochem1D() {
    // atmosphere
    auto mech = Cantera::newSolution("photolysis_o2.yaml");
    auto atm = std::make_shared<AtmChemistry>("atm", mech);

    iT = atm->componentIndex("T");
    iP = atm->componentIndex("P");
    iN2 = atm->componentIndex("N2");
    iO = atm->componentIndex("O");
    iO2 = atm->componentIndex("O2");
    iO3 = atm->componentIndex("O3");

    // grid
    double uin = 0.;
    double T = 300;
    double P = 101325;

    int nz = 10;
    double height = 100.E3;
    std::vector<double> z(nz);
    double dz = height / nz;
    for (int iz = 0; iz < nz; iz++) {
      z[iz] = dz / 2. + iz * dz;
    }

    // resize function is called inside setupGrid
    atm->setupGrid(nz, z.data());

    // surface
    std::string X = "O2:0.21 N2:0.79";
    auto surface = std::make_shared<SurfaceBoundary>("surface", mech);
    surface->setTemperature(T);

    // space
    auto space = std::make_shared<SpaceBoundary>("space", mech);
    surface->setTemperature(T);

    // set up simulation
    // pchem = new AtmChemistrySimulator({surface, atm, space});
    pchem = new AtmChemistrySimulator({atm});
    pchem->initFromFile("stellar/sun.ir");
  }

  ~TestPhotochem1D() { delete pchem; }
};

TEST_F(TestPhotochem1D, check_domain_index) {
  ASSERT_EQ(iT, 0);
  ASSERT_EQ(iP, 1);
  ASSERT_EQ(iN2, 3);
  ASSERT_EQ(iO, 4);
  ASSERT_EQ(iO2, 6);
  ASSERT_EQ(iO3, 7);

  // int dom = static_cast<int>(pchem->domainIndex("surface"));
  // ASSERT_EQ(dom, 0);
  int iatm = pchem->domainIndex("atm");
  ASSERT_EQ(iatm, 0);
  int ncomp = pchem->domain(iatm).nComponents();
  ASSERT_EQ(ncomp, 8);
  int npoints = pchem->domain(iatm).nPoints();
  ASSERT_EQ(npoints, 10);
  int size = pchem->size();
  ASSERT_EQ(size, npoints * ncomp);
  // dom = static_cast<int>(pchem->domainIndex("space"));
  // ASSERT_EQ(dom, 2);
}

TEST_F(TestPhotochem1D, check_profile) {
  int iatm = pchem->domainIndex("atm");
  double vN2 = 0.70;
  double vO = 0.1;
  double vO2 = 0.21;
  double vO3 = 0.1;

  pchem->setFlatProfile(iatm, iN2, vN2);
  pchem->setFlatProfile(iatm, iO, vO);
  pchem->setFlatProfile(iatm, iO2, vO2);
  pchem->setFlatProfile(iatm, iO3, vO3);

  for (int j = 0; j < pchem->domain(iatm).nPoints(); j++) {
    ASSERT_DOUBLE_EQ(pchem->value(iatm, iN2, j), vN2);
    ASSERT_DOUBLE_EQ(pchem->value(iatm, iO, j), vO);
    ASSERT_DOUBLE_EQ(pchem->value(iatm, iO2, j), vO2);
    ASSERT_DOUBLE_EQ(pchem->value(iatm, iO3, j), vO3);
  }
}

TEST_F(TestPhotochem1D, check_attenuation) {
  int iatm = pchem->domainIndex("atm");
  double T0 = 300.;
  double P0 = 101325.;

  pchem->setFlatProfile(iatm, iT, T0);
  pchem->setFlatProfile(iatm, iP, P0);
}

TEST_F(TestPhotochem1D, check_time_step) {
  int iatm = pchem->domainIndex("atm");
  double vN2 = 0.79;
  double vO2 = 0.21;
  double T0 = 300.;
  double P0 = 0.01;

  pchem->setFlatProfile(iatm, iN2, vN2);
  pchem->setFlatProfile(iatm, iO2, vO2);
  pchem->setFlatProfile(iatm, iT, T0);
  pchem->setFlatProfile(iatm, iP, P0);

  pchem->show();

  int nsteps = 10;
  double dt = 1.0;

  pchem->timeStep(1, dt, 8);
}

int main(int argc, char **argv) {
  Application::Start(argc, argv);

  testing::InitGoogleTest(&argc, argv);

  int result = RUN_ALL_TESTS();

  Application::Destroy();

  return result;
}
