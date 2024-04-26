#ifndef SRC_INTERPOLATION_HPP_
#define SRC_INTERPOLATION_HPP_

// C/C++ headers
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// C3M headers
#include <configure.hpp>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

Eigen::MatrixXd InterpolateCrossSection(Eigen::MatrixXd reference_wavelength,
                                        Eigen::MatrixXd input_wavelength,
                                        Eigen::MatrixXd input_cross_section);
Eigen::MatrixXd InterpolateQYield(Eigen::MatrixXd reference_wavelength,
                                  Eigen::MatrixXd input_wavelength,
                                  Eigen::MatrixXd qyield);

#endif  // SRC_INTERPOLATION_HPP_
