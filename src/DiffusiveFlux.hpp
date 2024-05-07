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

void DiffusiveFlux(const std::vector<std::vector<double>>& Atm, const std::vector<std::vector<double>>& Ni, std::vector<std::vector<double>>& flux);


#endif  // SRC_CUSTOMRATE_HPP_
