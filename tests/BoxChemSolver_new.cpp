// @sec3{Include files}
// Author: Ananyo Bhattacharya
// Affiliation: University of Michigan
// Email: ananyo@umich.edu
// The code computes the chemical evolution of a reaction network at a point

// Solution class describes a phase consists of a mixture of chemical species
#include <cantera/base/Solution.h>

// ThermoPhase object stores the thermodynamic state
#include <cantera/thermo.h>

// Kinetics object stores the chemical kinetics information
#include <cantera/kinetics.h>
#include <cantera/kinetics/Reaction.h>

// output stream
#include <cantera/kinetics/KineticsFactory.h>
#include <cantera/numerics/eigen_dense.h>
#include <cantera/numerics/eigen_sparse.h>
#include <cantera/thermo/ThermoFactory.h>

#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

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
  // Reading input file
  IOWrapper infile;
  infile.Open("test_athena.inp", IOWrapper::FileMode::read);
  ParameterInput *pinput = new ParameterInput();

  pinput->LoadFromFile(infile);
  infile.Close();

  // Loading the input parameters for atmospheric profile and reaction network
  std::string network_file = pinput->GetString("problem", "network");
  std::cout << network_file << std::endl;
  // Reading the chemical kinetics network
  auto sol = newSolution("ch4_photolysis.yaml");
  // auto sol = newSolution(network_file);
  auto gas = sol->thermo();
  auto gas_kin = sol->kinetics();
  int nsp = gas->nSpecies();
  int nrxn = gas_kin->nReactions();
  std::cout << "Network imported" << std::endl;
  // Initial condition for mole fraction
  VectorXd mole_fractions = VectorXd::Zero(nsp);

  // Loading the input parameters for initial condition
  std::string init_species_list = pinput->GetString("init", "species");
  std::regex pattern("[a-zA-Z0-9_]+");
  std::smatch m;
  int species_inx;

  while (std::regex_search(init_species_list, m, pattern)) {
    for (auto x : m) {
      std::string species_init_condition = pinput->GetString("init", x);
      species_inx = gas->speciesIndex(x);
      mole_fractions(species_inx) = atof(species_init_condition.c_str());
    }
    init_species_list = m.suffix().str();
  }
  std::cout << "Initial conditions loaded" << std::endl;

  // Loading the atmospheric conditions
  std::string tp = pinput->GetString("problem", "temp");
  std::string ps = pinput->GetString("problem", "pres");
  double temp = atof(tp.c_str());
  double pres = atof(ps.c_str());

  // Integrator
  std::string time_step = pinput->GetString("integrator", "dt");
  std::string max_time = pinput->GetString("integrator", "Tmax");
  double dt = atof(time_step.c_str());
  double Tmax = atof(max_time.c_str());
  double Ttot = 0.0;
  std::cout << "Time step: " << dt << std::endl;

  std::cout << "All inputs loaded into C3M " << std::endl;
  VectorXd dQ(nsp);
  Eigen::SparseMatrix<double> m_wjac;
  m_wjac.resize(nsp, nsp);
  MatrixXd mat1(nsp, nsp);
  mat1 = MatrixXd::Identity(nsp, nsp);
  MatrixXd mat2(nsp, nsp);
  VectorXd m_wdot(nsp);
  // Calculating photochemical reaction rate
  std::string stellar_input_file = pinput->GetString("radtran", "solar");
  std::string radius = pinput->GetString("radtran", "radius");
  std::string reference = pinput->GetString("radtran", "reference");

  double rad = atof(radius.c_str());
  double ref = atof(reference.c_str());
  std::cout << rad << " " << ref << std::endl;
  MatrixXd stellar_input =
      ReadStellarRadiationInput(stellar_input_file, rad, ref);
  std::cout << "Radiation Input Complete!" << std::endl;

  // Updating the actinic flux within yaml file [All in SI units]
  gas_kin->setWavelength(stellar_input.row(0).data(),
                         stellar_input.row(0).size());
  gas_kin->updateActinicFlux(stellar_input.row(1).data());

  // Solving the chemical reactions
  while (Ttot < Tmax) {
    // Setting T, P, X for each grid point
    gas->setState_TP(temp, (pres)*OneBar);
    gas->setMoleFractions(&mole_fractions[0]);
    // Solving the net production for each species
    gas_kin->getNetProductionRates(
        &m_wdot[0]);  // Extracting net production rates from Cantera
    m_wjac =
        gas_kin->netProductionRates_ddX();  // Extracting Jacobian from Cantera
    m_wjac = m_wjac / gas->molarDensity();
    // Integration for each species
    // Backward Euler Scheme (Li and Chen, 2019)
    mat2 = ((mat1) - (m_wjac * dt / gas->molarDensity()));
    mat2 = mat2.inverse();
    mat2 = dt * mat2 * ((m_wdot));
    dQ = mat2 / gas->molarDensity();
    Ttot = Ttot + dt;
    mole_fractions = mole_fractions + dQ;
    dt = dt * 1.25;
    std::cout << "Simulation completed at time: " << Ttot << std::endl;
  }

  std::cout << "*** Simulation complete! ***" << std::endl;
  delete pinput;
}
