// @sec3{Include files}
// Author: Ananyo Bhattacharya
// Affiliation: University of Michigan
// Email: ananyo@umich.edu
// The code computes the chemical evolution of a reaction network at a point


// Solution class describes a phase consists of a mixture of chemical species
#include <cantera/base/Solution.h>

// ThermoPhase object stores the thermodynamic state
#include <cantera/thermo.h>
#include <cantera/thermo/Phase.h>
// Kinetics object stores the chemical kinetics information
#include <cantera/kinetics.h>
#include <cantera/kinetics/Reaction.h>

// output stream
#include <iostream>
#include <fstream>
#include <cantera/numerics/eigen_dense.h>
#include <cantera/numerics/eigen_sparse.h>
#include <cantera/zeroD/Reactor.h>
#include <cantera/zeroD/IdealGasConstPressureReactor.h>
#include <cantera/zeroD/ReactorNet.h>
#include <cantera/numerics/Integrator.h>
#include <cantera/numerics/CVodesIntegrator.h>
#include <vector>
#include <string>
#include <sstream>
#include <regex>


// Athena++ header
#include <parameter_input.hpp>

// C3M header
#include <PhotoChemistry.hpp>
#include <RadTran.hpp>
#include <interpolation.hpp>


// NetCDF Output
#if NETCDFOUTPUT
  #include <netcdf.h>
#endif

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace Cantera;
using namespace std;


int main(int argc, char **argv) {

//Reading input file
  IOWrapper infile;
  infile.Open("test_athena.inp", IOWrapper::FileMode::read);
  ParameterInput *pinput = new ParameterInput();

  pinput->LoadFromFile(infile);
  infile.Close();

//Loading the input parameters for atmospheric profile and reaction network
  std::string network_file = pinput->GetString("problem", "network");

//Reading the chemical kinetics network
  auto sol = newSolution(network_file);
  //auto sol = newSolution(network_file);
  auto gas = sol->thermo();
  auto gas_kin = sol->kinetics();
  int nsp = gas->nSpecies();
  int nrxn = gas_kin->nReactions();
  std::cout << "Network imported" << std::endl;
//Initial condition for mole fraction
  VectorXd mole_fractions = VectorXd::Zero(nsp);
  VectorXd ChemConc = VectorXd::Zero(nsp);
  VectorXd m_wdot = VectorXd::Zero(nsp);

//Loading the input parameters for initial condition
  std::string init_species_list = pinput->GetString("init", "species");
  std::regex pattern ("[a-zA-Z0-9_]+");
  std::smatch m;
  int species_inx;
  while (std::regex_search (init_species_list,m,pattern)) {
    for (auto x:m){
    std::string species_init_condition = pinput->GetString("init", x);
    species_inx = gas->speciesIndex(x);
    mole_fractions(species_inx) = atof(species_init_condition.c_str());
    }
    init_species_list = m.suffix().str();
  }
 std::cout << "Initial conditions loaded" << std::endl;

//Loading the atmospheric conditions
  std::string tp = pinput->GetString("problem", "temp");
  std::string ps = pinput->GetString("problem", "pres");
  double temp = atof(tp.c_str());
  double pres = atof(ps.c_str());

//Integrator
  std::string time_step = pinput->GetString("integrator", "dt");
  std::string max_time = pinput->GetString("integrator", "Tmax");
  double dt = atof(time_step.c_str());
  double Tmax = atof(max_time.c_str());
  double Ttot = 0.0;
  std::cout << "Time step: " << dt << std::endl;

  std::cout << "All inputs loaded into C3M " << std::endl;
  //Photochemistry
  std::string photochem = pinput->GetString("problem", "photochem");

  std::cout << "Starting photochemistry calculations" << std::endl;
//Radiative transfer
//Reading the Stellar Irradiance Input
  std::string stellar_input_file = pinput->GetString("radtran", "solar");
  std::string radius = pinput->GetString("radtran", "radius");
  std::string reference =   pinput->GetString("radtran", "reference");

  double rad = atof(radius.c_str());
  double ref = atof(reference.c_str());
  std::cout << rad << " " << ref << std::endl;
  MatrixXd stellar_input = ReadStellarRadiationInput(stellar_input_file, rad, ref);
  std::cout << "Radiation Input Complete!" << std::endl;

  std::cout << "Starting chemical evolution " << std::endl;
//Setting up Cantera reactor object
  gas->setMoleFractions(&mole_fractions[0]);
  gas->setState_TP(temp, (pres)*OneBar);
  double P = gas->pressure();
  //std::cout << mole_fractions.transpose() << std::endl;
  //gas->setMoleFractions(&mole_fractions[0]);
  P = gas->pressure();
 // std::cout << "Pressure: " << P/1E5 << ", " << pres << std::endl;
 double totalDensity = gas->molarDensity();
 // std::cout << totalDensity*1e3*6.022e23*1e-6 << std::endl;
  //gas->setMoleFractions(&mole_fractions[0]);
  ReactorNet sim;
  IdealGasConstPressureReactor reactor;
  reactor.setThermoMgr(*gas);
  reactor.setKineticsMgr(*gas_kin);
  reactor.setEnergy(false);
  reactor.initialize(0.0);
  double rtol = 1.0e-15;  // relative tolerance
  double atol = 1.0e-30; // absolute tolerance
  sim.addReactor(reactor);
  sim.setTolerances(rtol, atol);
  sim.setSensitivityTolerances(rtol, atol);
  int nmax = 1E9;
  Integrator& integrator = sim.integrator();
  integrator.setMaxSteps(nmax);
  sim.initialize();

//Setting photochemical reaction rates
  if(photochem == "true"){
  for(int irxn = 0; irxn < nrxn; irxn++){
    auto& rxnObj = *(gas_kin->reaction(irxn));
    std::string rxnEqn = rxnObj.equation();
    std::cout << rxnEqn << std::endl;
    int pos = rxnEqn.find("=");
    rxnEqn.replace(pos, 2, "->");
    std::string jrate = pinput->GetOrAddString("Jrates", rxnEqn, "nan");
    if(jrate != "nan"){
    gas_kin->setMultiplier(irxn, atof(jrate.c_str()));
    std::cout << rxnEqn << " " << jrate << std::endl;
    }

    } }

    Ttot = dt;

  std::string OutFileName = pinput->GetString("output", "file");
  std::ofstream outfile (OutFileName);
//Integrating the reaction rates

/*
reactor.initialize();
sim.initialize();
totalDensity = gas->molarDensity();
//std::cout << totalDensity*1e3*6.022e23*1e-6 << std::endl;
sim.advance(Tmax);
gas->getMoleFractions(&mole_fractions[0]);
gas->getConcentrations(&ChemConc[0]);
gas_kin->getNetProductionRates(&m_wdot[0]);
std::cout << mole_fractions.transpose() << std::endl;
//std::cout << ChemConc.transpose()*1e3*6.022e23*1e-6 << std::endl;
//std::cout << m_wdot.transpose() << std::endl;
double mm  = gas->meanMolecularWeight();
P = gas->pressure();
double T = gas->temperature();
std::cout << mm << std::endl;
std::cout << P/1E5 << std::endl;
std::cout << T << std::endl;
*/

while(Ttot  < Tmax) {
    init_species_list = pinput->GetString("init", "species");
    while (std::regex_search (init_species_list,m,pattern)) {
      for (auto x:m){
        std::string species_init_condition = pinput->GetString("init", x);
        species_inx = gas->speciesIndex(x);
        mole_fractions(species_inx) = atof(species_init_condition.c_str());
    }
    init_species_list = m.suffix().str();
  }

    gas->setMoleFractions(&mole_fractions[0]);
    gas->setState_TP(temp, (pres)*OneBar);
    reactor.initialize();
    sim.initialize();
    sim.advance(dt);
    double rho_new = gas->density();
    cout << "Time: " << sim.time() << " s" << endl;
    dt = dt*1.25;
    gas->getMoleFractions(&mole_fractions[0]);
    //std::cout << mole_fractions.transpose() << std::endl;
    VectorXd krate(nrxn);
    gas_kin->getFwdRateConstants(&krate[0]);
    std::cout << krate.transpose() << std::endl;
    std::cout << mole_fractions.transpose() << std::endl;
    outfile << Ttot << " " << mole_fractions.transpose()  << std::endl;
    Ttot = Ttot + dt;
   }
 //  double R = GasConstant;
//   double NA = Avogadro;
//   double n_total = (pres*1e5/(R*temp))*NA;
//  totalDensity = gas->density();
//  std::cout << R << " " << NA << " " << pres << " " << temp <<  std::endl;
//  std::cout << n_total*1e-6 << std::endl;
//  std::cout << totalDensity << std::endl;
  std::cout << "Simulation Complete!" << std::endl;
}


