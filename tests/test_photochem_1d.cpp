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
  int iN2, iO, iO2, iO3;

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
    double uin = .1;
    double T = 300;
    double P = 101325;

    int nz = 50;
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
    surface->setMoleFractions(X);

    double rho_in = surface->solution()->thermo()->density();
    surface->setMdot(uin * rho_in);
    surface->setTemperature(T);

    // space
    auto space = std::make_shared<SpaceBoundary>("space", mech);
    double rho_out = space->solution()->thermo()->density();
    double uout = surface->mdot() / rho_out;

    // set up simulation
    // pchem = new Cantera::OneDim({surface, atm, space});
    pchem = new AtmChemistrySimulator({atm});
  }

  ~TestPhotochem1D() { delete pchem; }
};

TEST_F(TestPhotochem1D, check_domain_index) {
  ASSERT_EQ(iN2, 3);
  ASSERT_EQ(iO, 4);
  ASSERT_EQ(iO2, 6);
  ASSERT_EQ(iO3, 7);

  // int dom = static_cast<int>(pchem->domainIndex("surface"));
  // ASSERT_EQ(dom, 0);
  int dom = pchem->domainIndex("atm");
  ASSERT_EQ(dom, 0);
  int ncomp = pchem->domain(dom).nComponents();
  ASSERT_EQ(ncomp, 8);
  int npoints = pchem->domain(dom).nPoints();
  ASSERT_EQ(npoints, 50);
  int size = pchem->size();
  ASSERT_EQ(size, npoints * ncomp);
  // dom = static_cast<int>(pchem->domainIndex("space"));
  // ASSERT_EQ(dom, 2);
}

TEST_F(TestPhotochem1D, check_profile) {
  int datm = pchem->domainIndex("atm");
  double vN2 = 0.70;
  double vO = 0.1;
  double vO2 = 0.21;
  double vO3 = 0.1;

  pchem->setFlatProfile(datm, iN2, vN2);
  pchem->setFlatProfile(datm, iO, vO);
  pchem->setFlatProfile(datm, iO2, vO2);
  pchem->setFlatProfile(datm, iO3, vO3);

  for (int j = 0; j < pchem->domain(datm).nPoints(); j++) {
    ASSERT_DOUBLE_EQ(pchem->value(datm, iN2, j), vN2);
    ASSERT_DOUBLE_EQ(pchem->value(datm, iO, j), vO);
    ASSERT_DOUBLE_EQ(pchem->value(datm, iO2, j), vO2);
    ASSERT_DOUBLE_EQ(pchem->value(datm, iO3, j), vO3);
  }

  pchem->show();
}

TEST_F(TestPhotochem1D, check_time_step) {
  int datm = pchem->domainIndex("atm");
  double vN2 = 0.79;
  double vO2 = 0.21;

  pchem->setFlatProfile(datm, iN2, vN2);
  pchem->setFlatProfile(datm, iO2, vO2);

  int nsteps = 10;
  double dt = 1.0;

  // pchem->timeStep(1, dt, 8);
}

int main(int argc, char **argv) {
  Application::Start(argc, argv);

  testing::InitGoogleTest(&argc, argv);

  int result = RUN_ALL_TESTS();

  Application::Destroy();

  return result;
}
