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
  int nsp = nSpecies();
  Cantera::Domain1D::resize(stride, points);

  m_do_species.resize(nsp, true);
  m_delta_z.resize(points, 0.);

  m_wtm.resize(points, 0.);
  m_hk.resize(nsp, points, 0.);

  m_Keddy.resize(points, 0.);
  m_Kbinary.resize(points);
  m_diff_flux.resize(nsp, points, 0.);
  m_diff_t.resize(nsp, points, 0.);
  for (size_t j = 0; j < points; j++) {
    m_Kbinary[j].resize(nsp, nsp);
    m_Kbinary[j].setZero();
  }

  m_Xmid.resize(nsp, 0.);

  m_wdot.resize(nsp, points, 0.);
  m_conv_t.resize(nsp, points, 0.);

  // tri-diagnoal matrices
  m_A.resize(points);
  m_B.resize(points);
  m_C.resize(points);
  m_D.resize(points);

  for (size_t j = 0; j < points; ++j) {
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
  m_dwtm.resize(points);

  for (size_t j = 0; j < points; ++j) {
    m_dwtm[j].resize(nsp);
  }
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

  for (size_t j = 1; j < nPoints(); j++) {
    m_Tmid[j] = 0.5 * (getT(j - 1) + getT(j));
    m_Pmid[j] = sqrt(getP(j - 1) * getP(j));
    thermo->setTemperature(m_Tmid[j]);
    thermo->setPressure(m_Pmid[j]);

    const double* Xjm = x + stride() * (j - 1);
    const double* Xj = x + stride() * j;
    for (size_t k = 0; k < nSpecies(); ++k) {
      m_Xmid[k] = 0.5 * (Xjm[k] + Xj[k]);
    }

    thermo->setMoleFractions_NoNorm(m_Xmid.data());
    for (size_t k = 0; k < nSpecies(); ++k) {
      m_dwtm[j](k) = thermo->molecularWeight(k) - thermo->meanMolecularWeight();
    }
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

  // udpate m_wdot, m_B
  updateReaction(rdt, x, jmin, jmax);

  // update m_conv_t, m_A, m_B, m_C
  updateConvection(x, jmin, jmax);

  // update m_diff_t, m_A, m_B, m_C
  updateDiffusion(x, jmin, jmax);

  // make compressed matrices and evaluate residual
  for (size_t j = jmin; j <= jmax; ++j) {
    m_A[j].makeCompressed();
    m_B[j].makeCompressed();
    m_C[j].makeCompressed();

    for (size_t k = 0; k < nSpecies(); ++k) {
      rsd[index(k, j)] = m_wdot(k, j) + m_diff_t(k, j) - m_conv_t(k, j);
      diag[index(k, j)] = m_do_species[k];
    }
  }

  // disable species at the boundaries
  for (size_t k = 0; k < nSpecies(); ++k) {
    for (size_t j = 0; j < jmin; ++j) {
      rsd[index(k, j)] = 0.;
      diag[index(k, j)] = 0;
    }

    for (size_t j = jmax + 1; j < nPoints(); ++j) {
      rsd[index(k, j)] = 0.;
      diag[index(k, j)] = 0;
    }
  }
}

void AtmChemistry::updateReaction(double rdt, double const* x, size_t j0,
                                  size_t j1) {
  auto thermo = solution()->thermo();
  auto kinetics = solution()->kinetics();
  int nsp = thermo->nSpecies();

  for (size_t j = j0; j <= j1; j++) {
    m_B[j].setZero();
    if (rdt > 0.) {
      for (size_t k = 0; k < nsp; ++k) m_B[j].insert(k, k) = rdt;
    } else {
      for (size_t k = 0; k < nsp; ++k) m_B[j].insert(k, k) = 1.E-10;
    }

    kinetics->setActinicFluxLevel(j);

    const double* X = x + stride() * j;
    thermo->setMoleFractions_NoNorm(X);
    kinetics->getNetProductionRates(&m_wdot(0, j));

    m_B[j] -= kinetics->netProductionRates_ddX();
  }
}

void setEddyDiffusionCoeff(double Kzz, size_t j){
   m_Keddy[j] = Kzz; 
}

double getEddyDiffusionCoeff(size_t j){
   double Kzz = m_Keddy[j];
   return Kzz
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
      m_conv_t(k, j) = uj * (getX(x, k, j - usgn) - getX(x, k, j)) / dz;
      m_A[j].insert(k, k) = relu(uj) / dz;
      m_C[j].insert(k, k) = -relu(-uj) / dz;
      m_B[j].coeffRef(k, k) -= uj / dz;
    }
  }
}

void AtmChemistry::updateDiffusion(const double* x, size_t j0, size_t j1) {
  auto thermo = solution()->thermo();
  int nsp = thermo->nSpecies();

  // eddy and molecular diffusion transport
  Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(nsp, nsp);

  for (size_t j = j0; j <= j1 + 1; j++) {
      m_Keddy[j] = getEddyDiffusionCoeff(j);
      m_Kbinary[j] = getBinaryDiffusionCoeff(x, j);

    Eigen::VectorXd Kflux = (m_Keddy[j] * identity + m_Kbinary[j]) * dXdz(x, j);
    for (size_t k = 0; k < nsp; k++) {
      m_diff_flux(k, j) = Kflux(k);
    }
  }

  for (size_t j = j0; j <= j1; j++) {
    for (size_t k = 0; k < nsp; k++) {
      m_diff_t(k, j) =
          (m_diff_flux(k, j + 1) - m_diff_flux(k, j)) / m_delta_z[j];

      double Kjm = (m_Keddy[j] + m_Kbinary[j](k, k)) /
                   (m_delta_z[j] * (m_z[j] - m_z[j - 1]));
      double Kjp = (m_Keddy[j + 1] + m_Kbinary[j + 1](k, k)) /
                   (m_delta_z[j] * (m_z[j + 1] - m_z[j]));

      Eigen::VectorXd Djm = 0.5 * m_Kbinary[j] * m_dwtm[j] * m_grav /
                            (m_Tmid[j] * Cantera::GasConstant * m_delta_z[j]);
      Eigen::VectorXd Djp =
          0.5 * m_Kbinary[j + 1] * m_dwtm[j + 1] * m_grav /
          (m_Tmid[j + 1] * Cantera::GasConstant * m_delta_z[j]);

      m_A[j].coeffRef(k, k) += -Kjm + Djm(k);
      m_C[j].coeffRef(k, k) += -Kjp - Djp(k);
      m_B[j].coeffRef(k, k) += Kjm + Kjp + Djm(k) - Djp(k);
    }
  }
}
