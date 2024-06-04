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

class TestJupiter : public testing::Test {
 public:
  // data
  AtmChemistrySimulator *pchem;
  int iH2, iHe;
  std::shared_ptr<std::vector<double>> hydro;

  // constructor
  TestJupiter() {
    // atmosphere
    auto mech = Cantera::newSolution("Jupiter_diffusion.yaml");
    auto atm = std::make_shared<AtmChemistry>("atm", mech);
    auto gas = mech->thermo();
    auto gas_kin = mech->kinetics();
    int nsp = gas->nSpecies();
    int nrxn = gas_kin->nReactions();
    Cantera::ThermoPhase* gasThermo = gas.get();
    iH2 = atm->componentIndex("H2");
    iHe = atm->componentIndex("He");
   
    //Using YAML to extract components of atmospheric structure
    YAML::Node config  = YAML::LoadFile("Jupiter_diffusion.yaml");
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
      atm->setEddyDiffusionCoeff(AtmData(iKzz, n), n);
      std::cout << "Pres: " << AtmData(iPress, n) << " & " << AtmData(iTemp,n) << std::endl;
      Eigen::VectorXd Diag  = handleCustomMolecularDiffusion(PlanetName, gasThermo, AtmData(iPress, n), AtmData(iTemp,n), mWt);
      
      Dzz.diagonal() = Diag;
      atm->setBinaryDiffusionCoeff(Dzz, n); 
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


    pchem->setFlatProfile(atm, iH2, 0.8);
    pchem->setFlatProfile(atm, iHe, 0.2);
    atm->setGravity(-24);

    std::string X = "H2:0.86 He:0.14";
    pchem->find<Connector>("surface")->setSpeciesDirichlet(X);
    std::string X_space = "H2:0.9 He:0.1";
    pchem->find<Connector>("space")->setSpeciesDirichlet(X_space);
  }

  ~TestJupiter() { delete pchem; }
};

/*
class TestJupiter : public testing::Test{
 public:
  // data
  AtmChemistrySimulator *pchem;
  std::shared_ptr<std::vector<double>> hydro;
  int iH2, iHe;

  // constructor
  TestJupiter() {
   // atmosphere 
    auto sol = Cantera::newSolution("Jupiter_diffusion.yaml");
    //auto gas = sol->thermo();
   // auto gas_kin = sol->kinetics();
    auto atm = std::make_shared<AtmChemistry>("atm", sol);  
   // int nsp = gas->nSpecies();
   // int nrxn = gas_kin->nReactions();
   // Cantera::ThermoPhase* gasThermo = gas.get();
    
    iH2 = atm->componentIndex("H2");
    iHe = atm->componentIndex("He");
    
   //Using YAML to extract components of atmospheric structure
    YAML::Node config  = YAML::LoadFile("Jupiter_diffusion.yaml");
//    YAML::Node Atm_File = config["problem"][0]["planet"];
//    YAML::Node Planet_Name = config["problem"][1]["name"];
//    std::string AtmFile = Atm_File.as<std::string>();
//    std::string PlanetName = Planet_Name.as<std::string>();
//    Eigen::VectorXd mWt(nsp);  // Molecular weight

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
    std::cout << "Grid set up done!" << std::endl;
    hydro = std::make_shared<std::vector<double>>(atm->nPoints() * 3);
    for (size_t n = 0; n < atm->nPoints(); ++n) {
      hydro->at(n * 3) = T0;
      hydro->at(n * 3 + 1) = P0;
      hydro->at(n * 3 + 2) = U0;
    }

    atm->setHydro(hydro, 3);
    std::cout << "Hydro set up done!" << std::endl;
// surface (lower boundary)
  auto surface = std::make_shared<SurfaceBoundary>("surface", sol);
  std::cout << "Surface set up done!" << std::endl;
// space  (upper boundary)
  auto space = std::make_shared<SpaceBoundary>("space", sol);
  std::cout << "Space set up done!" << std::endl;
// set up simulation
  AtmChemistrySimulator *pchem;
  pchem = new AtmChemistrySimulator({surface, atm, space});
  pchem->initFromFile("stellar/sun.ir");


    pchem->setFlatProfile(atm, iH2, 1.0);
    pchem->setFlatProfile(atm, iHe, 0.0);

    std::string X = "H2:0.86 He:0.14";
    pchem->find<Connector>("surface")->setSpeciesDirichlet(X);
    pchem->find<Connector>("space")->setSpeciesDirichlet(X);

  } };
*/

/*
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
    InFile.open(AtmFile);
    getline(InFile, data1);
    getline(InFile, data1);
    while (InFile >> data1 >> data2 >> data3 >> data4 >> data5){
      AtmData(iPress, inx) = atof(data1.c_str())*1E2; //Pressure (mbar) -> Pa
      AtmData(iTemp,inx) = atof(data2.c_str()); //Temperature (K)
      AtmData(iKzz, inx) = atof(data3.c_str())*1E-4; //Kzz (cm^2/s) -> m^2/s
      AtmData(iAlt, inx) = atof(data4.c_str())*1E3; //Altitude (m) 
      AtmData(iNd, inx) = atof(data5.c_str())*1E6/(6.022E23*1E3); //Concentration (1/cm^3) -> kmol/m^3
      inx++;
      }
  InFile.close();



// resize function is called inside setupGrid
  atm->setupGrid(nz, AtmData.row(iAlt).data());
  std::shared_ptr<std::vector<double>> hydro;
  hydro = std::make_shared<std::vector<double>>(atm->nPoints() * 3);
  for (int insp = 0; insp < nsp; insp++) {
          mWt(insp) = gas->molecularWeight(insp);  // Kg/kml
        }
        
  for (size_t n = 0; n < atm->nPoints(); ++n) {
      hydro->at(n * 3) = AtmData(iTemp,n);
      hydro->at(n * 3 + 1) = AtmData(iPress, n);
      hydro->at(n * 3 + 2) = U0;
      std::cout << "seems good" << AtmData(iTemp,n)  << std::endl;
      //atm->setEddyDiffusionCoeff(AtmData(iKzz, n), n);
      Eigen::MatrixXd Dzz = Eigen::MatrixXd::Zero(n,n);
      Dzz.diagonal() = handleCustomMolecularDiffusion(PlanetName, gasThermo, AtmData(iPress, n), AtmData(iTemp,n), mWt);
      //atm->setBinaryDiffusionCoeff(Dzz, n); 
    }

  atm->setHydro(hydro, 3);
  std::cout << "Loaded all inputs" << std::endl;
//Homogeneous input condition
  pchem->setFlatProfile(atm, iH2, 1.0);
  std::cout << "Set up H2 mixing ratio profile " << std::endl;
  pchem->setFlatProfile(atm, iHe, 0.0);

  std::cout << "Set up input profiles" << std::endl;

  std::string X_bot = "H2:0.86 He:0.14";
  pchem->find<Connector>("surface")->setSpeciesDirichlet(X_bot);
  
  std::string X_top = "H2:0.0 He:0.0";
  pchem->find<Connector>("space")->setSpeciesNeumann(X_top);
  */









TEST_F(TestJupiter, check_basic) {
  auto atm = pchem->find<>("atm");

}

TEST_F(TestJupiter, check_abundance) {
  auto atm = pchem->find<>("atm");
 
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






