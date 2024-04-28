#ifndef SRC_RADTRAN_HPP_
#define SRC_RADTRAN_HPP_

// C/C++ headers
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

Eigen::MatrixXd ReadVULCANPhotoAbsCrossSection(string Cross_Section_File);
std::pair<std::vector<double>, std::vector<double>> ReadStellarRadiationInput(
    string Solar_Input_File, double rad, double ref);

#endif  // SRC_RADTRAN_HPP_
