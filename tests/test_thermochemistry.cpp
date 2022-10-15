// @sec3{Include files}
// Solution class describes a phase consists of a mixture of chemical species
#include <cantera/base/Solution.h>


// ThermoPhase object stores the thermodynamic state
#include <cantera/thermo.h>

// Kinetics object stores the chemical kinetics information
#include <cantera/kinetics.h>

// output stream
#include <iostream>
#include <cantera/numerics/eigen_dense.h>
#include <cantera/numerics/eigen_sparse.h>
#include <vector>

using Eigen::MatrixXd;
using Eigen::VectorXd;
// defines Real
//#include <defs.hpp>

using namespace Cantera;

// ordering of species in the enum class must be the same as those
// in the chemical mechanism file
enum {
  PH3 = 0, H3PO4, HOPO, H, PH2, H2O, H2, OH, PO, PO2, HPO, HOPO2,
    HPOH, H2POH, He
};

// @sec3{Main program}
int main(int argc, char **argv)
{
  // Read in chemical reactions. The first argument is the input chemistry
  // mechanism file. The second argument is the name of the phase.
  // newSolution returns std::shared_ptr<Solution>, which can be automatically
  // deduced using the "auto" keyword from C++11
  auto sol = newSolution("P_Reaction_Set_Simplified.yaml", "p_reaction_set");

  // obtain the underlying ThermoPhase object. ThermoPhase derives from
  // Phase, and augments it with thermodynamic properties of the substances
  // in addition to the ones that define the thermodynamic state such as
  // temperature, pressure and composition
  auto gas = sol->thermo();
  
  auto gas_kin = sol->kinetics();

  // check out number of species
  std::cout << "Number of species = "
            << gas->nSpecies() << std::endl;

  // check out species names given its index
  std::cout << "Species names = ";
  for (int i = 0; i < gas->nSpecies(); ++i) {
    std::cout << gas->speciesName(i) << " ";
  }
  std::cout << std::endl;

  // check out index of species given its name
  std::cout << "Index of H2O = "
            << gas->speciesIndex("H2O") << std::endl;

  double *mole_fractions = new double [gas->nSpecies()];
  mole_fractions[H2O] = 9.67E-4;
  mole_fractions[H2] = 0.87;
  mole_fractions[He] = 0.1;
  mole_fractions[PH3] = 2.3E-6;

  gas->setState_TP(1000., 300.0*OneAtm);
  gas->setMoleFractions(mole_fractions);
  VectorXd m_wdot(gas->nSpecies());
  float dt = 0.1;
  //gas_kin->getEquilibriumConstants();
  Eigen::SparseMatrix<double>  m_wjac;
  m_wjac.resize(gas->nSpecies(), gas->nSpecies());
  gas_kin->getNetProductionRates(&m_wdot[0]);
  m_wjac = gas_kin->netProductionRates_ddX();
  MatrixXd mat1(gas->nSpecies(), gas->nSpecies());
  MatrixXd mat2(gas->nSpecies(), gas->nSpecies());
  mat1 = MatrixXd::Identity(gas->nSpecies(), gas->nSpecies());
  mat2 = ((mat1/dt) - m_wjac);
  mat2 = mat2.inverse();
  mat2 = mat2*m_wdot;


  // Calculate chemical equilibrium under constant temperature
  // and pressure
  // Net production rates of Chemicals
  //std::cout << "Net Production Rates = ";         
  //for (int i = 0; i < m_wdot.size(); i++) {
  //  std::cout << m_wdot.at(i) << " ";
  //}
  
  
  std::cout << std::endl;
  std::cout << "net production = ";         
  std::cout << mat2 << " ";
  std::cout << std::endl;
            
  gas->equilibrate("TP");
  
  

  // print a report of chemical mixture
  std::cout << gas->report() << std::endl;

  delete[] mole_fractions;
}
