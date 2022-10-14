// C/C++ headers
#include <iostream>

// C3M headers
#include "PhotoChemistry.hpp"

PhotoChemistry::PhotoChemistry()
{
  std::cout << "I'm PhotoChemistry constructor" << std::endl;
  data = new Real[10];

  wavelengths_.push_back(1.);
  wavelengths_.push_back(2.);
}

PhotoChemistry::~PhotoChemistry()
{
  std::cout << "I'm PhotoChemistry destructor" << std::endl;
  delete[] data;
}

void PhotoChemistry::printWavelength()
{
  for (size_t i = 0; i < wavelengths_.size(); ++i) {
    std::cout << wavelengths_[i] << std::endl;
  }
}


