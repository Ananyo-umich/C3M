// C/C++ headers
#include <vector>

// Cantera headers
#include <cantera/kinetics/Kinetics.h>

// C3M headers
#include <configure.hpp>

class PhotoChemistry : public Cantera::Kinetics {
public:
  // function
  PhotoChemistry();
  ~PhotoChemistry();

  void printWavelength();

protected:
  // data
  std::vector<Real> wavelengths_;
  Real *data;
};
