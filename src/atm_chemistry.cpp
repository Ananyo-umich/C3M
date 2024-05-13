// cantera
#include <cantera/base/Solution.h>
#include <cantera/kinetics/Kinetics.h>
#include <cantera/thermo/ThermoPhase.h>
#include <cantera/transport/Transport.h>

// c3m
#include "actinic_flux.hpp"
#include "atm_chemistry.hpp"

template <typename T>
inline T relu(T x) {
  return x > 0. ? x : 0;
}

// Domain1D constructor will resize the solution vector to the correct size
AtmChemistry::AtmChemistry(std::string const& id,
                           std::shared_ptr<Cantera::Solution> sol,
                           size_t npoints, size_t stride)
    : Cantera::Domain1D(
          stride == Cantera::npos ? sol->thermo()->nSpecies() : stride,
          npoints) {
  setSolution(sol);
  setID(id);
}

std::string AtmChemistry::domainType() const { return "AtmChemistry"; }

void AtmChemistry::resetBadValues(double* xg) {
  auto thermo = solution()->thermo();

  // loc() returns the offset of the local solution vector in the global vector
  for (size_t j = 0; j < nPoints(); j++) {
    double* X = xg + loc() + stride() * j;
    // sets internal state to X, get rid of bad values
    thermo->setMoleFractions(X);
    // retrieve the mass fractions
    thermo->getMoleFractions(X);
  }
}

size_t AtmChemistry::size() const { return nSpecies() * nPoints(); }

void AtmChemistry::resize(size_t stride, size_t points) {
  // debug
  int nsp = nSpecies();
  Cantera::Domain1D::resize(stride, points);

  m_do_species.resize(nsp, true);
  m_delta_z.resize(points, 0.);

  m_wtm.resize(points, 0.);
  m_hk.resize(nsp, points, 0.);

  m_diff_flux.resize(nsp, points, 0.);
  m_diff_t.resize(nsp, points, 0.);

  m_Xmid.resize(nsp, 0.);

  m_wdot.resize(nsp, points, 0.);
  m_conv_t.resize(nsp, points, 0.);

  // tri-diagnoal matrices
  m_A.resize(points);
  m_B.resize(points);
  m_C.resize(points);
  m_D.resize(points);

  for (int j = 0; j < points; ++j) {
    m_A[j].resize(nsp, nsp);
    m_B[j].resize(nsp, nsp);
    m_C[j].resize(nsp, nsp);
    m_D[j].resize(nsp);
    m_D[j].setZero();
  }

  // half-grid
  m_Xmid.resize(points);
  m_Tmid.resize(points);
  m_Pmid.resize(points);
  m_dwtm.resize(nsp, points);
}

void AtmChemistry::setupGrid(size_t n, const double* z) {
  Domain1D::setupGrid(n, z);
  m_delta_z[0] = m_z[1] - m_z[0];
  size_t points = nPoints();
  for (int j = 1; j < points - 1; j++) {
    m_delta_z[j] = (m_z[j + 1] - m_z[j - 1]) / 2.;
  }
  m_delta_z[points - 1] = m_z[points - 1] - m_z[points - 2];
}

double AtmChemistry::initialValue(size_t n, size_t j) { return 0.; }

std::string AtmChemistry::componentName(size_t n) const {
  return solution()->thermo()->speciesName(n);
}

void AtmChemistry::setMidpoint(const double* x, size_t j) {
  auto thermo = solution()->thermo();

  m_Tmid[j] = 0.5 * (getT(j) + getT(j + 1));
  m_Pmid[j] = sqrt(getP(j) * getP(j + 1));
  thermo->setTemperature(m_Tmid[j]);
  thermo->setPressure(m_Pmid[j]);

  const double* Xj = x + stride() * j;
  const double* Xjp = x + stride() * (j + 1);
  for (size_t k = 0; k < nSpecies(); k++) {
    m_Xmid[k] = 0.5 * (Xj[k] + Xjp[k]);
  }

  thermo->setMoleFractions_NoNorm(m_Xmid.data());
}

void AtmChemistry::update(double const* x) {
  auto thermo = solution()->thermo();

  for (size_t j = 0; j < nPoints(); j++) {
    thermo->setTemperature(getT(j));
    thermo->setPressure(getP(j));

    const double* X = x + stride() * j;
    thermo->setMoleFractions_NoNorm(X);
    m_wtm[j] = thermo->meanMolecularWeight();

    // debug
    // std::cout << getT(j) << ", " << getP(j) << std::endl;

    thermo->getPartialMolarEnthalpies(&m_hk(0, j));
  }
}

void AtmChemistry::eval(size_t jGlobal, double* xGlobal, double* rsdGlobal,
                        integer* diagGlobal, double rdt) {
  // If evaluating a Jacobian, and the global point is outside the domain of
  // influence for this domain, then skip evaluating the residual
  if (jGlobal != Cantera::npos &&
      (jGlobal + 1 < firstPoint() || jGlobal > lastPoint() + 1)) {
    return;
  }

  // start of local part of global arrays
  double* x = xGlobal + loc();
  double* rsd = rsdGlobal + loc();
  integer* diag = diagGlobal + loc();

  // jmin and jmax does not include boundary points
  size_t jmin, jmax;
  if (jGlobal == Cantera::npos) {  // evaluate all interior points
    jmin = 1;
    jmax = nPoints() - 2;
  } else {  // evaluate points for Jacobian
    size_t jpt = (jGlobal == 0) ? 0 : jGlobal - firstPoint();
    jmin = std::max<size_t>(jpt, 2) - 1;
    jmax = std::min(jpt + 1, nPoints() - 2);
  }

  updateProperties(x, rdt, jmin, jmax);
  evalResidual(x, rsd, diag, rdt, jmin, jmax);
}

void AtmChemistry::updateProperties(double const* x, double rdt, size_t j0,
                                    size_t j1) {
  // properties are computed for grid points from j0 to j1
  updateThermo(x, rdt, j0, j1);
  updateReaction(x, j0, j1);
  updateConvection(x, j0, j1);
  // updateDiffusion(y, j0, j1);

  for (int j = j0; j <= j1; ++j) {
    m_A[j].makeCompressed();
    m_B[j].makeCompressed();
    m_C[j].makeCompressed();
  }
}

void AtmChemistry::updateThermo(double const* x, double rdt, size_t j0,
                                size_t j1) {
  auto thermo = solution()->thermo();
  size_t nsp = nSpecies();

  for (size_t j = j0; j <= j1; j++) {
    thermo->setTemperature(getT(j));
    thermo->setPressure(getP(j));

    const double* X = x + stride() * j;
    thermo->setMoleFractions_NoNorm(X);

    m_wtm[j] = thermo->meanMolecularWeight();
    thermo->getPartialMolarEnthalpies(&m_hk(0, j));

    m_B[j].setZero();
    if (rdt > 0.) {
      for (size_t k = 0; k < nsp; ++k) m_B[j].insert(k, k) = rdt;
    } else {
      for (size_t k = 0; k < nsp; ++k) m_B[j].insert(k, k) = 1.;
    }
  }
}

void AtmChemistry::updateReaction(double const* x, size_t j0, size_t j1) {
  auto thermo = solution()->thermo();
  auto kinetics = solution()->kinetics();

  for (size_t j = j0; j <= j1; j++) {
    kinetics->setActinicFluxLevel(j);

    const double* X = x + stride() * j;
    thermo->setMoleFractions_NoNorm(X);
    kinetics->getNetProductionRates(&m_wdot(0, j));

    m_B[j] -= kinetics->netProductionRates_ddX();
  }
}

void AtmChemistry::updateConvection(double const* x, size_t j0, size_t j1) {
  auto thermo = solution()->thermo();
  int nsp = thermo->nSpecies();

  // upwind convective transport
  for (size_t j = j0; j <= j1; j++) {
    double uj = getU(j);
    int usgn = uj > 0. ? 1 : -1;
    double dz = m_z[j - usgn] - m_z[j];

    m_A[j].setZero();
    m_C[j].setZero();
    for (size_t k = 0; k < nsp; ++k) {
      m_conv_t(k, j) = -uj * (getX(x, k, j - usgn) - getX(x, k, j)) / dz;
      m_A[j].insert(k, k) = relu(uj) / dz;
      m_C[j].insert(k, k) = -relu(-uj) / dz;
      m_B[j].coeffRef(k, k) -= uj / dz;
    }
  }
}

/*void AtmChemistry::updateDiffusion(const double* y, size_t j0, size_t j1) {
  std::cout << "I'm calling updateDiffusion" << std::endl;
  auto thermo = solution()->thermo();
  int nsp = thermo->nSpecies();

  // eddy and molecular transport
  for (size_t j = j0; j < j1; j++) {
    double dz = z(j + 1) - z(j);
    for (size_t k = 0; k < nsp; k++) {
      m_diff_flux(k, j) = getDiffusionFlux(y, j);
      m_diff_flux(k, j + 1) = getDiffusionFlux(y, j + 1);
      m_diff(k, j) = (m_diff_flux(k, j + 1) - m_diff_flux(k, j)) / m_delta_z[j];
      m_Keddy[j] = getEddyDiffusionCoeff(y, j);
      m_Kbinary[j] = getBinaryDiffusionCoeff(y, j);
    }

  for (size_t j = j0; j < j1; j++) {
    auto Kjm = (m_Keddy[j] + m_Kbinary[j]) / (m_delta_z[j] * (m_z[j] - m_z[j -
1])); auto Kjp = (m_Keddy[j + 1] + m_Kbinary[j + 1]) / (m_delta_z[j] * (m_z[j +
1] - m_z[j]));

    auto Djm = 0.5 * m_KBinary[j] * m_dwtm[j] * m_grav
                   / (m_Tmid[j] * Cantera::GasConstant * m_delta_z[j]);
    auto Djp = 0.5 * kBinary[j + 1] * m_dwtm[j + 1] * m_grav
                   / (m_Tmid[j + 1] * Cantera::GasConstant * m_delta_z[j]);

    m_A[j] += Kjm - Djm;
    m_C[j] += Kjp + Djp;
    m_B[j] += - Kjm - Kjp - Djm + Djp;
  }
}*/

void AtmChemistry::evalResidual(double* x, double* rsd, int* diag, double rdt,
                                size_t j0, size_t j1) {
  for (size_t j = j0; j <= j1; j++) {
    for (size_t k = 0; k < nSpecies(); k++) {
      rsd[index(k, j)] = m_wdot(k, j) + m_diff_t(k, j) - m_conv_t(k, j);
      diag[index(k, j)] = m_do_species[k];
    }
  }
}
