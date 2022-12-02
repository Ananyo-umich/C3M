// C/C++ headers
#include <iostream>
#include <fstream>

// C3M headers
#include "PhotoChemistry.hpp"
using Eigen::MatrixXd;
using Eigen::VectorXd;


//Print wavelengths
void printWavelength(Eigen::VectorXd wavelengths_)
{
  for (int i = 0; i < wavelengths_.size(); ++i) {
    std::cout << wavelengths_[i] << std::endl;
  }
}

//Print Cross Sections
void printCrossSection(Eigen::VectorXd crossSection_)
{

for (int i = 0; i < crossSection_.size(); ++i) {
    std::cout << crossSection_[i] << std::endl;
  }
}

// Calculate the photochemical reaction rate based on wavelength, cross section and actinic flux
double PhotoChemRate(Eigen::MatrixXd wavelengths_,Eigen::MatrixXd crossSection_, Eigen::MatrixXd actinicFlux_)
{double rate = 0.0;
for (int i = 1; i < crossSection_.size() ; ++i) {
    rate = rate + ((wavelengths_(i) - wavelengths_(i-1))*crossSection_(i)*actinicFlux_(i));
  }
  return rate;
}


// Function to read the photochemical cross sections from VULCAN database
// The output file will contain wavelength (nm), photoabsorption cross section (cm^2), photodissociation cross section (cm^2) and photoionization cross section (cm^2)
Eigen::MatrixXd  ReadVULCANCrossSection(string VULCAN_ID){
  fstream InFile;
  InFile.open(VULCAN_ID); 
  string wavlength;
  string photoabs;
  string photodiss;
  string photoion; 
  int num = 0;
  int rows = 0;
  while (getline(InFile, wavlength))
  rows++;
  InFile.close();
  
  std::cout << "rows = " << rows << std::endl;
  Eigen::MatrixXd Output(4, rows-1);
  InFile.open(VULCAN_ID); 
  getline(InFile,wavlength);
  while(getline(InFile,wavlength, ',')){
  getline(InFile,photoabs,',');
  getline(InFile,photodiss,',');
  getline(InFile,photoion,'\n');
  
  
  Output(0, num) = atof(wavlength.c_str());
  Output(1, num) = atof(photoabs.c_str());
  Output(2, num) = atof(photodiss.c_str());
  Output(3, num) = atof(photoion.c_str());
  
  num++;
  }
  InFile.close();
  return Output;
}

// Function to read cross section from VPL photochemical database
// The output file contains the wavelength (Angstroms) and cross section (cm^2)
Eigen::MatrixXd ReadAtmosCrossSection(string AtmosFileName){
  fstream InFile;
  InFile.open(AtmosFileName); 
  string wavlength;
  string photoabs;
  int num = 0;
  int rows = 0;
  while (getline(InFile, wavlength))
  rows++;
  InFile.close();
  
  
  Eigen::MatrixXd Output(2, rows-4);
  InFile.open(AtmosFileName); 
  getline(InFile,wavlength);
  getline(InFile,wavlength);
  getline(InFile,wavlength);
  getline(InFile,wavlength);
  while(InFile >> wavlength >> photoabs){
  
  Output(0, num) = atof(wavlength.c_str());
  Output(1, num) = atof(photoabs.c_str());

  
  num++;
  }
  InFile.close();
  return Output;
  
  
}




