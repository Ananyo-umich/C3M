// external
#include <gtest/gtest.h>

// application
#include <application/application.hpp>

// cantera
#include "cantera/base/ct_defs.h"
#include <cantera/base/Solution.h>
#include <cantera/kinetics/Reaction.h>
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
    auto mech = Cantera::newSolution("photolysis_o2.yaml");
    auto atm = std::make_shared<AtmChemistry>("atm", mech);
    auto gas = mech->thermo();
    auto gas_kin = mech->kinetics();
    int nsp = gas->nSpecies();
    int nrxn = gas_kin->nReactions();
    Cantera::ThermoPhase* gasThermo = gas.get();
    iN2 = atm->componentIndex("N2");
    iO2 = atm->componentIndex("O2");
    std::cout << "Created a mechanism" << std::endl;
    //Using YAML to extract components of atmospheric structure
    YAML::Node config  = YAML::LoadFile("photolysis_o2.yaml");
    YAML::Node Atm_File = config["problem"][0]["planet"];
    YAML::Node Planet_Name = config["problem"][1]["name"];
    std::string AtmFile = Atm_File.as<std::string>();
    std::string PlanetName = Planet_Name.as<std::string>();
    for (int irxn = 0; irxn < nrxn; irxn++){
          auto& rxnObj = *(gas_kin->reaction(irxn));
          std::string rxnEquation = rxnObj.equation();
          std::cout << rxnEquation << std::endl;
	}

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
      std::cout << "Pres: " << AtmData(iPress, inx) << ", Temp: " << AtmData(iTemp,inx) << std::endl;
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
    atm->setGraivty(9.8);

    hydro = std::make_shared<std::vector<double>>(atm->nPoints() * 3);
    for (size_t n = 0; n < atm->nPoints(); ++n) {
      hydro->at(n * 3) = AtmData(iTemp,n); 
      hydro->at(n * 3 + 1) = AtmData(iPress,n);
      hydro->at(n * 3 + 2) = U0;
      atm->setEddyDiffusionCoeff(AtmData(iKzz, n), n);
    }
    
    atm->setGravity(-9.8);
    atm->setHydro(hydro, 3);
    std::cout << "Hydro done!" << std::endl;

    for (size_t n = 0; n < atm->nPoints(); ++n) {
    std::cout << n << ": " << atm->getP(n) << std::endl; }
    // surface
    auto surface = std::make_shared<SurfaceBoundary>("surface", mech);
    std::cout << "Surface BC done!" << std::endl;
    // space
    auto space = std::make_shared<SpaceBoundary>("space", mech);
    std::cout << "Space BC done!" << std::endl;

    // set up simulation
    pchem = new AtmChemistrySimulator({space, atm, surface});
    pchem->initFromFile("stellar/sun.ir");
    std::cout << "Set Stellar conditions" << std::endl;

    // Setting up initial profiles

    pchem->setFlatProfile(atm, iO2, 0.2);
    pchem->setFlatProfile(atm, iN2, 0.8);
   
    //Initializing BC for surface and space domains
    Eigen::VectorXd SurfaceBC = Eigen::VectorXd::Ones(nsp)*0.0;
    Eigen::VectorXd SpaceBC = Eigen::VectorXd::Ones(nsp)*0.0;
    

    std::string X = "O2:0.21 N2:0.79";
    pchem->find<Connector>("surface")->setSpeciesDirichlet(X);
    std::string X2 = "O2:0.0 N2:0.0";
    pchem->find<Connector>("space")->setSpeciesNeumann(X2);
    
   
    //pchem->find<Connector>("surface")->setDirichletBC(SurfaceBC);
    //pchem->find<Connector>("space")->setDirichletBC(SpaceBC);
    double dt = 1E-10;
    
   // pchem->advance(10);
    pchem->timeStep(100, dt, 10);
    pchem->show();

    }

  ~TestEarthChem() { delete pchem; }
};


TEST_F(TestEarthChem, check_abundance) {
  auto atm = pchem->find<AtmChemistry>("atm");
 
 // pchem->setMaxTimeStep(1.E9);
 // pchem->show();

 // int nsteps = 10;
 // double dt = 1.0;



}


int main(int argc, char **argv) {
  Application::Start(argc, argv);

  testing::InitGoogleTest(&argc, argv);

  int result = RUN_ALL_TESTS();

  Application::Destroy();

  return result;
}






