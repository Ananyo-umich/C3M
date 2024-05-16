#ifndef SRC_BOUNDARY_HPP_
#define SRC_BOUNDARY_HPP_

// C/C++
#include <vector>

// cantera
#include <cantera/base/Solution.h>
#include <cantera/oneD/Domain1D.h>

struct DirichletBC {
  size_t id;
  double value;
};

struct NeumannBC {
  size_t id;
  double value;
};

class Connector : public Cantera::Domain1D {
 public:
  Connector(size_t nv = 1, size_t points = 1) : Domain1D(nv, points, 0.) {}

  std::string domainType() const override { return "connector"; }

  bool isConnector() override { return true; }

  //! Initial value of solution component @e n at grid point @e j.
  double initialValue(size_t n, size_t j) override;

  //! Name of the nth component.
  std::string componentName(size_t n) const override;

  //! Set the mole fractions by specifying a string.
  void setSpeciesDirichlet(const std::string& xin);

  //! Set the mole fractions by specifying an array.
  void setSpeciesNeumann(const std::string& xin);

  //! The eval function of a connector modifies the tri-diagonal matrices
  void eval(size_t jGlobal, double* yGlobal, double* rsdGlobal,
            integer* diagGlobal, double rdt) override;

 protected:
  std::vector<DirichletBC> dirichlet;
  std::vector<NeumannBC> neumann;
};

class SurfaceBoundary : public Connector {
 public:
  SurfaceBoundary() = default;

  SurfaceBoundary(std::string const& id,
                  std::shared_ptr<Cantera::Solution> sol);

  std::string domainType() const override { return "surface"; }
};

class SpaceBoundary : public Connector {
 public:
  SpaceBoundary() = default;

  SpaceBoundary(std::string const& id, std::shared_ptr<Cantera::Solution> sol);

  std::string domainType() const override { return "space"; }
};

#endif  // SRC_BOUNDARY_HPP_
