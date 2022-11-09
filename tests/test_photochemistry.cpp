//Photochemistry Header File
#include <PhotoChemistry.hpp>

// Standard C++ Header Files
#include <string>
#include <stdlib.h>
#include <iostream>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

MatrixXd Actinic_Flux(MatrixXd wav){
  VectorXd acf(wav.size());
  for (int i = 0; i < wav.size() ; i++) {
     acf(i) = 100; 
  }

  return acf;
}

int main(int arg, char **argc) {
  //Photochemical cross section file name
  string FileName = "/data4/ananyo/models/C3M/data/ATMOS/SO/SO.XS.dat";
  //Photochemical data
  MatrixXd Data = ReadAtmosCrossSection(FileName);
  MatrixXd wavelength = Data.row(0);
  MatrixXd Xsection = Data.row(1);
  MatrixXd Act_Flux = Actinic_Flux(wavelength);
  double Rate;
  Rate = PhotoChemRate(wavelength*1E-10, Xsection*1E-4, Act_Flux);
  std::cout << "The photochemical reaction rate is: " << Rate << std::endl;
  
}
