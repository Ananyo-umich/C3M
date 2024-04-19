#include <iostream>
#include <string>

// Athena++ header
#include <parameter_input.hpp>

int main(int argc, char **argv) {

  IOWrapper infile;
  infile.Open("test_athena.inp", IOWrapper::FileMode::read);

  ParameterInput *pinput = new ParameterInput();

  pinput->LoadFromFile(infile);
  infile.Close();

  std::string atm_file = pinput->GetString("problem", "planet");
  std::string network_file = pinput->GetString("problem", "network");

  std::cout << atm_file << std::endl;
  std::cout << network_file << std::endl;

  delete pinput;
}
