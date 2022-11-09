// @sec3{Include files}
// Solution class describes a phase consists of a mixture of chemical species
#include <cantera/base/Solution.h>
#include <PhotoChemistry.hpp>


// ThermoPhase object stores the thermodynamic state
#include <cantera/thermo.h>

// Kinetics object stores the chemical kinetics information
#include <cantera/kinetics.h>
#include <cantera/kinetics/Reaction.h>

// output stream
#include <iostream>
#include <cantera/numerics/eigen_dense.h>
#include <cantera/numerics/eigen_sparse.h>
#include <vector>
#include <string>
#include <stdlib.h>
#include <fstream>
using Eigen::MatrixXd;
using Eigen::VectorXd;
// defines Real
//#include <defs.hpp>

using namespace Cantera;
// ordering of species in the enum class must be the same as those
// in the chemical mechanism file
enum {
  O2 = 0, O, O3, O_1D
};

MatrixXd Actinic_Flux(MatrixXd wav){
  VectorXd acf(wav.size());
  for (int i = 0; i < wav.size() ; i++) {
     acf(i) = 100; 
  }

  return acf;
}


// @sec3{Main program}
int main(int argc, char **argv)
{
  // Read in chemical reactions. The first argument is the input chemistry
  // mechanism file. The second argument is the name of the phase.
  // newSolution returns std::shared_ptr<Solution>, which can be automatically
  // deduced using the "auto" keyword from C++11
  auto sol = newSolution("chapman.yaml");
  // obtain the underlying ThermoPhase object. ThermoPhase derives from
  // Phase, and augments it with thermodynamic properties of the substances
  // in addition to the ones that define the thermodynamic state such as
  // temperature, pressure and composition
  auto gas = sol->thermo();
  auto gas_kin = sol->kinetics();
  int nsp = gas->nSpecies();
  
  // check out species names given its index
  std::cout << "Species names = ";
  for (int i = 0; i < gas->nSpecies(); ++i) {
    std::cout << gas->speciesName(i) << " ";
  }
  std::cout << std::endl;
  
  //Initial condition for mole fraction
  double *mole_fractions = new double [gas->nSpecies()];
  mole_fractions[O2] = 0.21;
  mole_fractions[O] = 1E-10;
  mole_fractions[O3] = 1E-6;
  mole_fractions[O_1D] = 1E-10;
  
  VectorXd mole_frac(gas->nSpecies());
  //Setting T, P and mole fraction
  gas->setState_TP(255, OneAtm);
  gas->setMoleFractions_NoNorm(mole_fractions);
  gas->getMoleFractions(&mole_frac[0]);  
  
  
  //Printing the Initial condition
  std::cout << std::endl;
  std::cout << "Initial condition = ";
  for (int i = 0 ; i < gas->nSpecies(); i++){
  std::cout << mole_fractions[i] << " ";}
  std::cout << std::endl;

  
  //Reading the cross sections from photochemical database
  MatrixXd Rxn0 = ReadVULCANCrossSection("/data4/ananyo/models/C3M/data/VULCAN/O2/O2_cross.csv");
  MatrixXd Rxn2 = ReadVULCANCrossSection("/data4/ananyo/models/C3M/data/VULCAN/O3/O3_cross.csv");
  MatrixXd Rxn3 = ReadVULCANCrossSection("/data4/ananyo/models/C3M/data/VULCAN/O3/O3_cross.csv");
  
  
  
  //gas_kin->setMultiplier(0, rate);
  //gas_kin->setMultiplier(5, rate2);
  //VectorXd mole_frac(gas->nSpecies());
  float dt = 1e-1; //Time step (s)
  int nSteps = 10000; //Number of time steps 
  
  //Initiating Matrices for Time Evolution
  VectorXd m_wdot(gas->nSpecies()); //Net Production Rate
  Eigen::SparseMatrix<double>  m_wjac; //Jacobian
  m_wjac.resize(gas->nSpecies(), gas->nSpecies()); 
  MatrixXd mat1(gas->nSpecies(), gas->nSpecies());
  MatrixXd mat2(gas->nSpecies(), gas->nSpecies());
  mat1 = MatrixXd::Identity(gas->nSpecies(), gas->nSpecies()); 
  MatrixXd Act_Flux;
  double rxnRate0, rxnRate2, rxnRate3;

  //Iterating over time steps
  for (int i = 0; i < nSteps; i++) {
    //Calculating photochemical reaction rate for all the photochemical reactions  
    Act_Flux = Actinic_Flux(Rxn0.row(0));
    rxnRate0 = PhotoChemRate(Rxn0.row(0)*1E-9, Rxn0.row(2)*1E-4, Act_Flux);
    Act_Flux = Actinic_Flux(Rxn2.row(0));
    rxnRate2 = PhotoChemRate(Rxn2.row(0)*1E-9, Rxn2.row(2)*1E-4, Act_Flux);
    Act_Flux = Actinic_Flux(Rxn3.row(0));
    rxnRate3 = PhotoChemRate(Rxn3.row(0)*1E-9, Rxn3.row(2)*1E-4, Act_Flux);
  
    // Setting the multiplier for each photochemical reaction
    gas_kin->setMultiplier(0, rxnRate0);
    gas_kin->setMultiplier(2, rxnRate2);
    gas_kin->setMultiplier(3, rxnRate3);
    
    //Extracting Net Production Rates from Cantera
    gas_kin->getNetProductionRates(&m_wdot[0]);
    //Extracting Jacobian from Cantera
    m_wjac = gas_kin->netProductionRates_ddX();

    //Backward Euler Scheme (Li and Chen, 2019)
    mat2 = ((mat1/dt) - m_wjac);
    mat2 = mat2.inverse();
    mat2 = mat2*m_wdot;
    //Updating the mole fractions
    gas->getMoleFractions(&mole_frac[0]);
    mole_frac = mole_frac + mat2;
    gas->setMoleFractions_NoNorm(&mole_frac[0]);
    std::cout << std::endl;
    std::cout << "Mole Fraction at t = " << (i+1)*dt << " (s) \n" ;
    for (int i = 0 ; i < gas->nSpecies(); i++){
    std::cout << mole_frac(i) << " ";}
    std::cout << std::endl;
  }

 
}
