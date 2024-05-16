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

void DiffusiveFlux(Eigen::MatrixXd Atm, Eigen::MatrixXd Xf,
                   Eigen::MatrixXd flux, Eigen::VectorXd alpha, double grav,
                   std::string PlanetName, std::string FileName);

#endif  // SRC_CUSTOMRATE_HPP_
