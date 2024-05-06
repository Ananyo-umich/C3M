#ifndef SRC_BOUNDARY_HPP_
#define SRC_BOUNDARY_HPP_

// cantera
#include <cantera/oneD/Boundary1D.h>

class SurfaceBoundary : public Cantera::Inlet1D {
 public:
  SurfaceBoundary() = default;

  SurfaceBoundary(std::string const& id, std::shared_ptr<Cantera::Solution> sol)
      : Cantera::Inlet1D(sol, id) {}

  std::string domainType() const override { return "surface"; }
};

class SpaceBoundary : public Cantera::Outlet1D {
 public:
  SpaceBoundary() = default;

  SpaceBoundary(std::string const& id, std::shared_ptr<Cantera::Solution> sol)
      : Cantera::Outlet1D(sol, id) {}

  std::string domainType() const override { return "space"; }
};

#endif  // SRC_BOUNDARY_HPP_
