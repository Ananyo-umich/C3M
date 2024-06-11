// external
#include <gtest/gtest.h>

// application
#include <application/application.hpp>

// cantera
#include <cantera/base/Solution.h>

// output stream
#include <cantera/numerics/eigen_dense.h>
#include <cantera/numerics/eigen_sparse.h>

// c3m
#include <c3m/atm_chemistry.hpp>
#include <c3m/atm_chemistry_simulator.hpp>
#include <c3m/boundary.hpp>
#include <c3m/CustomTransport.hpp>


#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>


// YAML file IO
#include "cantera/base/ct_defs.h"
#include "cantera/ext/yaml-cpp/yaml.h"

class TestEarthChem : public testing::Test {
 public:
  // data
  AtmChemistrySimulator *pchem;
  int iCO2, iN2, iH2O, iCH4, iNO2, iNO, iO2;
  std::shared_ptr<std::vector<double>> hydro;

  // constructor
  TestEarthChem() {
    // atmosphere
    auto mech = Cantera::newSolution("Earth_CHNOS.yaml");
    auto atm = std::make_shared<AtmChemistry>("atm", mech);
    auto gas = mech->thermo();
    auto gas_kin = mech->kinetics();
    int nsp = gas->nSpecies();
    int nrxn = gas_kin->nReactions();
    Cantera::ThermoPhase* gasThermo = gas.get();
    iCO2 = atm->componentIndex("CO2");
    iN2 = atm->componentIndex("N2");
    iO2 = atm->componentIndex("O2");
    iH2O = atm->componentIndex("H2O");
    iCH4 = atm->componentIndex("CH4");
    iNO2 = atm->componentIndex("NO2");
    iNO = atm->componentIndex("NO");
    std::cout << "Created a mechanism" << std::endl;
    //Using YAML to extract components of atmospheric structure
    YAML::Node config  = YAML::LoadFile("Earth_CHNOS.yaml");
    YAML::Node Atm_File = config["problem"][0]["planet"];
    YAML::Node Planet_Name = config["problem"][1]["name"];
    std::string AtmFile = Atm_File.as<std::string>();
    std::string PlanetName = Planet_Name.as<std::string>();
    Eigen::VectorXd mWt(nsp);  // Molecular weight
    for (int insp = 0; insp < nsp; insp++) {
          mWt(insp) = gas->molecularWeight(insp);  // Kg/kml
        }
  

    int nz = 0;
    std::string data1, data2, data3, data4, data5;
    std::fstream InFile;
    InFile.open(AtmFile);
    getline(InFile, data1);
    getline(InFile, data1);
    while (getline(InFile, data1))
    nz++;
    InFile.close();

    Eigen::MatrixXd AtmData(5, nz);
    int iTemp = 0;
    int iPress =  1;
    int iKzz =  2;
    int iAlt = 3;
    int iNd = 4;

    int inx = 0;
    std::vector<double> z(nz);
    InFile.open(AtmFile);
    getline(InFile, data1);
    getline(InFile, data1);
    while (InFile >> data1 >> data2 >> data3 >> data4 >> data5){
      AtmData(iPress, inx) = atof(data1.c_str())*1E2; //Pressure (mbar) -> Pa
      AtmData(iTemp,inx) = atof(data2.c_str()); //Temperature (K)
      AtmData(iKzz, inx) = atof(data3.c_str())*1E-4; //Kzz (cm^2/s) -> m^2/s
      z[inx] = atof(data4.c_str())*1E3; //Altitude (m)
      AtmData(iAlt, inx) = atof(data4.c_str())*1E3; //Altitude (m) 
      AtmData(iNd, inx) = atof(data5.c_str())*1E6/(6.022E23*1E3); //Concentration (1/cm^3) -> kmol/m^3
      inx++;
      }
  InFile.close();
     std::cout << "Extracted input from model Earth atmosphere" << std::endl;
    // grid
    double T0 = 300;
    double P0 = 0.01;
    double U0 = 0.;

    Eigen::MatrixXd Dzz = Eigen::MatrixXd::Zero(nsp,nsp);
    // resize function is called inside setupGrid
    atm->setupGrid(nz, z.data());

    hydro = std::make_shared<std::vector<double>>(atm->nPoints() * 3);
    for (size_t n = 0; n < atm->nPoints(); ++n) {
      hydro->at(n * 3) = AtmData(iTemp,n);
      hydro->at(n * 3 + 1) = AtmData(iPress, n);
      hydro->at(n * 3 + 2) = U0;
      atm->setEddyDiffusionCoeff(0.0, n);
      std::cout << "Pres: " << AtmData(iPress, n) << " & " << AtmData(iTemp,n) << std::endl;
      //Eigen::VectorXd Diag  = handleCustomMolecularDiffusion(PlanetName, gasThermo, AtmData(iPress, n), AtmData(iTemp,n), mWt);
      
      //Dzz.diagonal() = Diag;
      atm->setBinaryDiffusionCoeff(Dzz, n); 
    }
    
    atm->setHydro(hydro, 3);
    std::cout << "Hydro done!" << std::endl;
    // surface
    auto surface = std::make_shared<SurfaceBoundary>("surface", mech);
    std::cout << "Surface BC done!" << std::endl;
    // space
    auto space = std::make_shared<SpaceBoundary>("space", mech);
    std::cout << "Space BC done!" << std::endl;

    // set up simulation
    pchem = new AtmChemistrySimulator({surface, atm, space});
    pchem->initFromFile("stellar/sun.ir");
    std::cout << "Set Stellar conditions" << std::endl;

    // Setting up initial profiles

    pchem->setFlatProfile(atm, iO2, 0.22);
    pchem->setFlatProfile(atm, iN2, 0.78);
    atm->setGravity(-9.8);
   
    //Initializing BC for surface and space domains
    Eigen::VectorXd SurfaceBC = Eigen::VectorXd::Ones(nsp)*1e-30;
    Eigen::VectorXd SpaceBC = Eigen::VectorXd::Ones(nsp)*1e-30;
    
    //Applying special BC for space and surface
    //Surface
    SurfaceBC(iN2) = 0.78; 
    SurfaceBC(iO2) = 0.21;
    SurfaceBC(iH2O) = 1E-6;
    SurfaceBC(iCO2) = 400E-6;

    //Space
    SpaceBC(iN2) = 0.78;
    SpaceBC(iO2) = 0.21;

    //pchem->find<Connector>("surface")->setDirichletBC(SurfaceBC);
    //pchem->find<Connector>("space")->setDirichletBC(SpaceBC);

    std::string X = "O2:1.0";
    pchem->find<Connector>("surface")->setSpeciesDirichlet(X);
    pchem->find<Connector>("space")->setSpeciesDirichlet(X);
    //std::string X = "N2:0.78 O2:0.22 H2O:1e-6 CO2:400e-6";
   // pchem->find<Connector>("surface")->setSpeciesDirichlet(X);
   // std::string X_space = "N2:0.78 O2:0.22";
   // pchem->find<Connector>("space")->setSpeciesDirichlet(X_space);

    }

  ~TestEarthChem() { delete pchem; }
};


TEST_F(TestEarthChem, check_basic) {
  auto atm = pchem->find<>("atm");

}

TEST_F(TestEarthChem, check_abundance) {
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






