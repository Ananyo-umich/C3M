// C/C++ headers
#include <stdlib.h>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// C3M headers
#include <configure.hpp>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

Eigen::MatrixXd ReadVULCANPhotoAbsCrossSection(string Cross_Section_File);
Eigen::MatrixXd ReadStellarRadiationInput(string Solar_Input_File, double rad,
                                          double ref);
