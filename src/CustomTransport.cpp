// @sec3{Include files}
// Author: Ananyo Bhattacharya
// Affiliation: University of Michigan
// Email: ananyo@umich.edu
// C/C++ headers
#include <vector>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <fstream>

// Athena++ header
#include <parameter_input.hpp>

// Cantera headers
#include <cantera/kinetics/Kinetics.h>
#include <cantera/base/ct_defs.h>
#include <cantera/kinetics/ReactionData.h>
#include <cantera/kinetics/ReactionRate.h>
#include <cantera/kinetics/MultiRate.h>
#include <cantera/base/Solution.h>
#include <cantera/thermo.h>

// C3M headers
#include <configure.hpp>
#include "CustomTransport.hpp"
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;


void handleCustomTransport(string PlanetName, Cantera::ThermoPhase* NetworkName, double Pres, double Temp, double diff){

  if(PlanetName == "JupiterAurora")
    {JupiterMolDiff(NetworkName, Pres, Temp, diff);
    }

}


void JupiterMolDiff(Cantera::ThermoPhase* NetworkName, double Pres, double Temp, double diff){
//Molecular diffusion in Hydrogen atmosphere [Egert et al., 2017]
 double Ndensity =  NetworkName->molarDensity()*1E3*6.022E23/1E6; // kmol/m^3 -> #/cm^3 conversion
 double MolDiff = 1.98E20*pow(Temp, 0.51)/(pow(log(5.34E6/Temp), 2)*Ndensity); //cm^2/s
 MolDiff = MolDiff*1E-4; //cm^2/s -> m^2/s
 diff = diff + MolDiff; 
}

