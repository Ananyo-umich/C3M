#ifndef SRC_CUSTOMRATE_HPP_
#define SRC_CUSTOMRATE_HPP_

// Author: Ananyo Bhattacharya
// Affiliation: University of Michigan
// Email: ananyo@umich.edu

// C/C++ headers
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Cantera headers
#include <cantera/base/Solution.h>
#include <cantera/base/ct_defs.h>
#include <cantera/kinetics/Kinetics.h>
#include <cantera/kinetics/MultiRate.h>
#include <cantera/kinetics/Reaction.h>
#include <cantera/kinetics/ReactionData.h>
#include <cantera/kinetics/ReactionRate.h>
#include <cantera/thermo.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

void handleCustomChemistry(string PlanetName, Cantera::Kinetics* NetworkName,
                           double Pres, double Temp);
void JupiterAuroralChemistry_VT(Cantera::Kinetics* NetworkName, double Pres,
                                double Temp);

#endif  // SRC_CUSTOMRATE_HPP_
