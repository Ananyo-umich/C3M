// @sec3{Include files}
// Solution class describes a phase consists of a mixture of chemical species
#include <cantera/base/Solution.h>

// ThermoPhase object stores the thermodynamic state
#include <cantera/thermo.h>

// output stream
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>

// defines Real
//#include <defs.hpp>
using namespace std;
using namespace Cantera;

// ordering of species in the enum class must be the same as those
// in the chemical mechanism file
enum {
  CO = 0, H2, H2O, O2, H2O2, CH4, H2CO, CH3OH, CO2, CH3OOH, C2H2,
    C2H4, C2H6, CH2CO, CH3CHO, C2H5OH, C2H5OOH, CH3COOOH, cC2H4O, C, CH,
    _1CH2, _3CH2, O3P, H, OH, OOH, CH3, HCO, CH2OH, CH3O, CH3OO, C2H, C2H3,
    C2H5, CHCO, CH2CHO, CH3CO, C2H5O, C2H4OOH, C2H5OO, CH3COOO, CH3OCO,
    CO2H, _1C2H4OH, _2C2H4OH, C3H8, C4H8Y, C4H10, C2H5CHO, C3H7OH, C3H7O,
    C4H9O, C2H6CO, C3H8CO, C2H3CHO, _1C3H7, _2C3H7, _1C4H9, _2C4H9, O1D, N2D,
    NO3, HONO2, CH3ONO, CH3NO2, HNO2, CH3NO, NO2, HONO, HCNN, HCNO, N2O,
    NCO, HNO, HOCN, NNH, H2CN, N4S, CN, HNCO, NO, NH, NH2, HCN, NH3, N2,
    N2O4, N2O3, N2H2, N2H3, N2H4, HNNO, HNOH, HNO3, NH2OH, H2NO, CNN, H2CNO,
    C2N2, HCNH, HON, NCN, He, CH3NH2, CH3NH, CH2NH2, CH2NH
};

// @sec3{Main program}
int main(int argc, char **argv)
{

  //Read the temperature and pressure from ECCM (Romani and Atreya) input file
  //col1: Z (km), col2: T (K), col3: P(bar), col4: NH3-H2O solution, col5: H2O ice
  //col6: NH4SH solid, col7: NH3 ice, col8: H2S ice, col9: CH4 ice, col10: Ar ice
  string Z, T, P, a4, a5, a6, a7, a8, a9, a10;
  double Temp[1560], Press[1560];
  int num = 0;
  fstream inFile;
  inFile.open("/home/ananyo/models/C3M/build/bin/tp_profile.txt");
   if(inFile.fail()) // checks to see if file opended 
    { 
      cout << "error" << endl; 
      return 1; // no point continuing if the file didn't open...
    } 
  while(inFile >> Z >> T >> P >> a4 >> a5 >> a6 >> a7 >> a8 >> a9 >> a10 ) //Reading till the end of input file
  {
    Temp[num] = atof(T.c_str());
    Press[num] = atof(P.c_str());
    num++;
         }
  inFile.close();
  
  // Read in chemical reactions. The first argument is the input chemistry
  // mechanism file. The second argument is the name of the phase.
  // newSolution returns std::shared_ptr<Solution>, which can be automatically
  // deduced using the "auto" keyword from C++11
  auto sol = newSolution("CNOH_reaction_network.yaml", "reaction network");

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

  //Initial condition for mole fraction of the system
  double *mole_fractions = new double [gas->nSpecies()];
  mole_fractions[H2O] = 4.805E-4;
  mole_fractions[H2] = 0.864;
  mole_fractions[He] = 0.157;
  mole_fractions[C] = 2.37E-3;
  

  //Iterating over each vertical grid point using ECCM input T, P profile
  for (int i = 0; i < num; i++){
  std::cout << "T, P at level: " << i<<std::endl;
  gas->setMoleFractions(mole_fractions); //Setting mole fraction
  gas->setState_TP(Temp[i], Press[i]*OneAtm/1.0132); //Setting the T, P from ECCM at point i
  gas->equilibrate("TP"); //Thermochemical equilibrium from Cantera
}

  //The output from thermochemical equilibrium can be extracted from the mole fraction
  //and stored in a vector or matrix format.

  delete[] mole_fractions;
}
