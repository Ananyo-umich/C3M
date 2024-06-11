// C/C++
#include <iostream>

// cantera
#include <cantera/base/stringUtils.h>
#include <cantera/thermo/ThermoPhase.h>

// c3m
#include "boundary.hpp"

double Connector::initialValue(size_t n, size_t j) { return 0.; }

std::string Connector::componentName(size_t n) const {
  return solution()->thermo()->speciesName(n);
}


void Connector::setDirichletBC(Eigen::VectorXd D_bc){
  auto mf = solution()->thermo();
  size_t ns = mf->nSpecies();
  for(size_t ins = 0; ins < ns; ins++){
    DirichletBC bc;
    bc.id = ins;
    bc.value = D_bc(ins);
    dirichlet.push_back(bc);
  }

}

  
void Connector::setNeumannBC(Eigen::VectorXd N_bc){
  auto mf = solution()->thermo();
  size_t ns = mf->nSpecies();
  for(size_t ins = 0; ins < ns; ins++){
    NeumannBC bc;
    bc.id = ins;
    bc.value = N_bc(ins);
    neumann.push_back(bc);
  } 


}

void Connector::setSpeciesDirichlet(const std::string& xin) {
  auto xmap = Cantera::parseCompString(xin);
  // will throw an exception if the species are not in the phase
  auto mf = solution()->thermo()->getCompositionFromMap(xmap);

  for (auto x : xmap) {
    DirichletBC bc;
    bc.id = componentIndex(x.first);
    bc.value = x.second;
    dirichlet.push_back(bc);
  }
}

void Connector::setSpeciesNeumann(const std::string& xin) {
  auto xmap = Cantera::parseCompString(xin);
  // will throw an exception if the species are not in the phase
  auto mf = solution()->thermo()->getCompositionFromMap(xmap);

  for (auto x : xmap) {
    NeumannBC bc;
    bc.id = componentIndex(x.first);
    bc.value = x.second;
    neumann.push_back(bc);
  }
}

void Connector::eval(size_t jGlobal, double* xGlobal, double* rsdGlobal,
                     integer* diagGlobal, double rdt) {
  // If evaluating a Jacobian, and the global point is outside the domain of
  // influence for this domain, then skip evaluating the residual
  if (jGlobal != Cantera::npos &&
      (jGlobal + 1 < firstPoint() || jGlobal > lastPoint() + 1)) {
    return;
  }

  if (left() != nullptr) {
    if (left()->isConnector()) {
      throw Cantera::CanteraError("Connector::eval",
                                  "Left domain cannot be a connector");
    }
  }

  if (right() != nullptr) {
    if (right()->isConnector()) {
      throw Cantera::CanteraError("Connector::eval",
                                  "Right domain cannot be a connector");
    }
  }

  Eigen::VectorXd vleft(nComponents());
  vleft.setZero();
  Eigen::VectorXd vright(nComponents());
  vright.setZero();

  // start of local part of global arrays
  double* x = xGlobal + loc();
  double* rsd = rsdGlobal + loc();

  // apply dirichlet boundary conditions
  for (auto bc : dirichlet) {
    x[bc.id] = bc.value;
    rsd[bc.id] = 0.;

    if (left() != nullptr) {
      size_t nv = left()->nComponents();
      double* xleft = xGlobal + left()->loc() + left()->size() - nv;
      xleft[bc.id] = bc.value;
      vleft(bc.id) = bc.value;
    }

    if (right() != nullptr) {
      double* xright = xGlobal + right()->loc();
      xright[bc.id] = bc.value;
      vright(bc.id) = bc.value;
    }
  }

  // apply neumann boundary conditions
  for (auto bc : neumann) {
    double *xleft, *xright;
    if (left() != nullptr) {
      size_t nv = left()->nComponents();
      size_t points = right()->nPoints();
      double dz = left()->z(points - 1) - left()->z(points - 2);

      xleft = xGlobal + left()->loc() + left()->size() - nv;
      xleft[bc.id] = xleft[bc.id + nv] - bc.value * dz;
      vleft(bc.id) = bc.value * dz;
      left()->B(points - 2) += left()->C(points - 2);
    }

    if (right() != nullptr) {
      size_t nv = right()->nComponents();
      double dz = right()->z(1) - right()->z(0);

      xright = xGlobal + right()->loc();
      xright[bc.id] = xright[bc.id + nv] - bc.value * dz;
      vright(bc.id) = -bc.value * dz;
      right()->B(1) += right()->A(1);
    }

    if (left() != nullptr && right() != nullptr) {
      x[bc.id] = 0.5 * (xleft[bc.id] + xright[bc.id]);
      rsd[bc.id] = xright[bc.id] - xleft[bc.id];
    } else if (left() != nullptr) {
      x[bc.id] = xleft[bc.id];
      rsd[bc.id] = 0.;
    } else if (right() != nullptr) {
      x[bc.id] = xright[bc.id];
      rsd[bc.id] = 0.;
    } else {
      throw Cantera::CanteraError("Connector::eval",
                                  "Neumann boundary condition without domain");
    }
  }

  if (left()) {
    size_t points = left()->nPoints();
    Eigen::VectorXd vD = left()->C(points - 2) * vleft;

    size_t start = left()->loc() + left()->size() - 2 * left()->nComponents();
    // for (size_t n = 0; n < nComponents(); ++n) {
    //   rsdGlobal[start + n] -= vD(n);
    // }
    left()->C(points - 2).setZero();
  }

  if (right()) {
    Eigen::VectorXd vD = right()->A(1) * vright;
    size_t start = right()->loc() + right()->nComponents();
    // for (size_t n = 0; n < nComponents(); ++n) {
    //   rsdGlobal[start + n] -= vD(n);
    // }
    right()->A(1).setZero();
  }
}

SurfaceBoundary::SurfaceBoundary(std::string const& id,
                                 std::shared_ptr<Cantera::Solution> sol)
    : Connector(sol->thermo()->nSpecies(), 1) {
  setSolution(sol);
  setID(id);
}

SpaceBoundary::SpaceBoundary(std::string const& id,
                             std::shared_ptr<Cantera::Solution> sol)
    : Connector(sol->thermo()->nSpecies(), 1) {
  setSolution(sol);
  setID(id);
}
