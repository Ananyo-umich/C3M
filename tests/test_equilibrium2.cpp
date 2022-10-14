// @sec3{Include files}
// Solution class describes a phase consists of a mixture of chemical species
#include <cantera/base/Solution.h>

// ThermoPhase object stores the thermodynamic state
#include <cantera/thermo.h>

// output stream
#include <iostream>

using namespace Cantera;

// ordering of species in the enum class must be the same as those
// in the chemical mechanism file
enum {
  iH2 = 0, iH = 1, iO = 2, iO2 = 3, iOH = 4, iH2O = 5,
  iHO2 = 6, iH2O2 = 7, iAR = 8, iN2 = 9
};

void equil_demo()
{
    // Create a new Solution object
    auto sol = newSolution("h2o2.yaml", "ohmech", "None");
    auto gas = sol->thermo();

    gas->setState_TPX(1500.0, 2.0*OneAtm, "O2:1.0, H2:3.0, AR:1.0");
    gas->equilibrate("TP");
    std::cout << gas->report() << std::endl;
}

// @sec3{Main program}
int main(int argc, char **argv)
{
  // Read in chemical reactions. The first argument is the input chemistry
  // mechanism file. The second argument is the name of the phase.
  // newSolution returns std::shared_ptr<Solution>, which can be automatically
  // deduced using the "auto" keyword from C++11
  auto sol = newSolution("h2o2.yaml", "ohmech");

  // obtain the underlying ThermoPhase object. ThermoPhase derives from
  // Phase, and augments it with thermodynamic properties of the substances
  // in addition to the ones that define the thermodynamic state such as
  // temperature, pressure and composition
  auto gas = sol->thermo();

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
  mole_fractions[iO2] = 0.2;
  mole_fractions[iH2] = 0.6;
  mole_fractions[iAR] = 0.2;

  gas->setState_TP(1500., 20.0*OneAtm);
  gas->setMoleFractions(mole_fractions);

  // Calculate chemical equilibrium under constant temperature
  // and pressure
  gas->equilibrate("TP");

  // print a report of chemical mixture
  std::cout << gas->report() << std::endl;

  delete[] mole_fractions;
}
