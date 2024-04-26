#ifndef SRC_VAPORS_HPP_
#define SRC_VAPORS_HPP_

// The header contains the saturation vapor pressure corresponding to different
// vapors
// C/C++ headers
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Cantera headers
#include <cantera/base/ct_defs.h>
#include <cantera/kinetics/Kinetics.h>
#include <cantera/kinetics/MultiRate.h>
#include <cantera/kinetics/ReactionData.h>
#include <cantera/kinetics/ReactionRate.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

#endif  // SRC_VAPORS_HPP_
