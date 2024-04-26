// The header contains the saturation vapor pressure corresponding to different
// vapors C/C++ headers
#include <stdlib.h>

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

// C3M headers
#include <configure.hpp>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
