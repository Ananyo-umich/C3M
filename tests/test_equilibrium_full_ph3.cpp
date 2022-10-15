// @sec3{Include files}
// Solution class describes a phase consists of a mixture of chemical species
#include <cantera/base/Solution.h>

// ThermoPhase object stores the thermodynamic state
#include <cantera/thermo.h>

// output stream
#include <iostream>

// defines Real
//#include <defs.hpp>

using namespace Cantera;

// ordering of species in the enum class must be the same as those
// in the chemical mechanism file
enum {
  PH3 = 0, H3PO4, HOPO, H, PH2, PH, H2O, H2, O2, O, OH, HO2, PO, PO2,
    PO3, HPO, HOPO2, P2O3, P, P2, P4, P2O, P2O2, HPOH, H2POH, He
};

// @sec3{Main program}
int main(int argc, char **argv)
{
  // Read in chemical reactions. The first argument is the input chemistry
  // mechanism file. The second argument is the name of the phase.
  // newSolution returns std::shared_ptr<Solution>, which can be automatically
  // deduced using the "auto" keyword from C++11
  auto sol = newSolution("P_Reaction_Set.yaml", "p_reaction_set");

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
  mole_fractions[H2O] = 4.805E-4;
  mole_fractions[H2] = 0.864;
  mole_fractions[He] = 0.157;
  mole_fractions[P] = 7.0E-7;//2.3E-6;
  
  gas->setMoleFractions(mole_fractions);

  // gas->setState_TP(204.73933649289103, 201691.45547303313);
  gas->setState_TP(86.25592417061611, 10181.517217181829);
  gas->equilibrate("TP");
  std::cout << gas->report() << std::endl;
  gas->setState_TP(150.23696682464453,68538.95838650088 );
  gas->equilibrate("TP");
  std::cout << gas->report() << std::endl;
  gas->setState_TP(204.73933649289103,201691.45547303313);
  gas->equilibrate("TP");
  std::cout << gas->report() << std::endl;
  gas->setState_TP(301.8957345971564, 791475.5439411161 );
  gas->equilibrate("TP");
  std::cout << gas->report() << std::endl;
  gas->setState_TP(399.0521327014218, 2090800.0412787204);
  gas->equilibrate("TP");
  std::cout << gas->report() << std::endl;
  gas->setState_TP(607.5829383886255, 9305720.40929699 );
  gas->equilibrate("TP");
  std::cout << gas->report() << std::endl;
  gas->setState_TP(989.0995260663508, 51396968.00771518 );
  gas->equilibrate("TP");
  std::cout << gas->report() << std::endl;
  gas->setState_TP(1197.6303317535544,99999999.9999999 );

  // Calculate chemical equilibrium under constant temperature
  // and pressure
  gas->equilibrate("TP");

  // print a report of chemical mixture
  std::cout << gas->report() << std::endl;

  delete[] mole_fractions;
}
