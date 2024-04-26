#ifndef SRC_CUSTOMTRANSPORT_HPP_
#define SRC_CUSTOMTRANSPORT_HPP_

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
#include <cantera/kinetics/ReactionData.h>
#include <cantera/kinetics/ReactionRate.h>
#include <cantera/thermo.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

VectorXd handleCustomMolecularDiffusion(string PlanetName,
                                        Cantera::ThermoPhase* NetworkName,
                                        double Pres, double Temp, VectorXd mWt);

VectorXd JupiterMolDiff(Cantera::ThermoPhase* NetworkName, double Pres,
                        double Temp, VectorXd mWt);

#endif  // SRC_CUSTOMTRANSPORT_HPP_
