// cantera
#include <cantera/base/Solution.h>
#include <cantera/kinetics/Kinetics.h>
#include <cantera/thermo/ThermoPhase.h>
#include <cantera/transport/Transport.h>

// c3m
#include "actinic_flux.hpp"
#include "atm_chemistry.hpp"

// Domain1D constructor will resize the solution vector to the correct size
AtmChemistry::AtmChemistry(std::string const& id,
                           std::shared_ptr<Cantera::Solution> sol,
                           size_t npoints)
    : Cantera::Domain1D(Offset::Y + sol->thermo()->nSpecies(), npoints) {
  setSolution(sol);
  setID(id);
}

std::string AtmChemistry::domainType() const { return "AtmChemistry"; }

std::string AtmChemistry::componentName(size_t n) const {
  int nsp = solution()->thermo()->nSpecies();
  switch (n) {
    case Offset::T:
      return "T";
    case Offset::P:
      return "P";
    case Offset::U:
      return "U";
    default:
      if (n >= Offset::Y && n < (Offset::Y + nsp)) {
        return solution()->thermo()->speciesName(n - Offset::Y);
      } else {
        return "<unknown>";
      }
  }
}

// m_nv is the number of components in the solution vector
void AtmChemistry::resetBadValues(double* xg) {
  auto thermo = solution()->thermo();

  // loc() returns the offset of the local solution vector in the global vector
  double* x = xg + loc();
  for (size_t j = 0; j < m_points; j++) {
    double* Y = x + m_nv * j + Offset::Y;
    // sets internal state to Y, get rid of bad values
    thermo->setMassFractions(Y);
    // retrieve the mass fractions
    thermo->getMassFractions(Y);
  }
}

void AtmChemistry::resize(size_t components, size_t points) {
  // debug
  std::cout << "I'm resizing the domain" << std::endl;
  std::cout << "components: " << components << std::endl;
  std::cout << "points: " << points << std::endl;
  Cantera::Domain1D::resize(components, points);
  m_rho.resize(points, 0.);
  m_wtm.resize(points, 0.);
  m_cp.resize(points, 0.);
  m_tcon.resize(points, 0.);
  m_fixedtemp.resize(points, 0.);

  int nsp = solution()->thermo()->nSpecies();
  m_ymid.resize(nsp, 0.);
  m_diff.resize(nsp * points, 0.);
  m_do_species.resize(components, true);

  m_hk.resize(nsp, points, 0.);
  m_wdot.resize(nsp, points, 0.);
  m_flux.resize(nsp, points, 0.);
}

void AtmChemistry::setMidpoint(const double* x, size_t j) {
  auto thermo = solution()->thermo();
  int nsp = thermo->nSpecies();

  thermo->setTemperature(0.5 * (getT(x, j) + getT(x, j + 1)));
  thermo->setPressure(sqrt(getP(x, j) * getP(x, j + 1)));

  const double* yyj = x + m_nv * j + Offset::Y;
  const double* yyjp = x + m_nv * (j + 1) + Offset::Y;
  for (size_t k = 0; k < nsp; k++) {
    m_ymid[k] = sqrt(yyj[k] * yyjp[k]);
  }

  thermo->setMassFractions_NoNorm(m_ymid.data());
}

void AtmChemistry::eval(size_t jGlobal, double* xGlobal, double* rsdGlobal,
                        integer* diagGlobal, double rdt) {
  std::cout << "I'm calling eval" << std::endl;
  std::cout << "jGlobal: " << jGlobal << std::endl;
  std::cout << "rdt: " << rdt << std::endl;

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

  size_t jmin, jmax;
  if (jGlobal == Cantera::npos) {  // evaluate all points
    jmin = 0;
    jmax = m_points - 1;
  } else {  // evaluate points for Jacobian
    size_t jpt = (jGlobal == 0) ? 0 : jGlobal - firstPoint();
    jmin = std::max<size_t>(jpt, 1) - 1;
    jmax = std::min(jpt + 1, m_points - 1);
  }

  std::cout << "jmin: " << jmin << std::endl;
  std::cout << "jmax: " << jmax << std::endl;

  updateProperties(jGlobal, x, jmin, jmax);

  evalContinuity(x, rsd, diag, rdt, jmin, jmax);

  if (m_do_energy) {
    evalEnergy(x, rsd, diag, rdt, jmin, jmax);
  } else {
    size_t j0 = std::max<size_t>(jmin, 1);
    size_t j1 = std::min(jmax, m_points - 2);
    for (size_t j = j0; j <= j1; j++) {
      rsd[index(Offset::T, j)] = 0.;
      rsd[index(Offset::P, j)] = 0.;
      diag[index(Offset::T, j)] = 0;
      diag[index(Offset::P, j)] = 0;
    }
  }

  evalSpecies(x, rsd, diag, rdt, jmin, jmax);
}

void AtmChemistry::updateProperties(size_t jg, double* x, size_t jmin,
                                    size_t jmax) {
  // properties are computed for grid points from j0 to j1
  size_t j0 = std::max<size_t>(jmin, 1) - 1;
  size_t j1 = std::min(jmax + 1, m_points - 1);

  updateThermo(x, j0, j1);
  updateReactionRates(x, j0, j1);

  if (jg == Cantera::npos || m_force_full_update) {
    // update transport properties only if a Jacobian is not being
    // evaluated, or if specifically requested
    // updateTransport(x, j0, j1);
  }

  int nsp = solution()->thermo()->nSpecies();

  if (jg == Cantera::npos) {
    double* Yleft = x + index(Offset::Y, jmin);
    m_kExcessLeft = std::distance(Yleft, std::max_element(Yleft, Yleft + nsp));

    double* Yright = x + index(Offset::Y, jmax);
    m_kExcessRight =
        std::distance(Yright, std::max_element(Yright, Yright + nsp));
  }

  // update the species diffusive mass fluxes whether or not a
  // Jacobian is being evaluated
  // updateDiffFluxes(x, j0, j1);
}

void AtmChemistry::evalContinuity(double* x, double* rsd, int* diag, double rdt,
                                  size_t jmin, size_t jmax) {
  std::cout << "I'm calling evalContinuity" << std::endl;
  // The left boundary has the same form for all cases.
  if (jmin == 0) {  // left boundary
    double dz = m_z[jmin + 1] - m_z[jmin];
    rsd[index(Offset::U, jmin)] =
        -(getRhoU(x, jmin + 1) - getRhoU(x, jmin)) / dz;
    diag[index(Offset::U, jmin)] = 0;  // Algebraic constraint
  }

  if (jmax == m_points - 1) {  // right boundary
    rsd[index(Offset::U, jmax)] = getRhoU(x, jmax) - getRhoU(x, jmax - 1);
    diag[index(Offset::U, jmax)] = 0;  // Algebraic constraint
  }

  // j0 and j1 are constrained to only interior points
  size_t j0 = std::max<size_t>(jmin, 1);
  size_t j1 = std::min(jmax, m_points - 2);

  // fixed mass flow rate
  for (size_t j = j0; j <= j1; j++) {
    rsd[index(Offset::U, j)] = getRhoU(x, j) - getRhoU(x, j - 1);
    diag[index(Offset::U, j)] = 0;  // Algebraic constraint
  }
}

void AtmChemistry::evalEnergy(double* x, double* rsd, int* diag, double rdt,
                              size_t jmin, size_t jmax) {
  std::cout << "I'm calling evalEnergy" << std::endl;
  if (jmin == 0) {  // left boundary
    rsd[index(Offset::T, jmin)] = getT(x, jmin);
  }

  if (jmax == m_points - 1) {  // right boundary
    rsd[index(Offset::T, jmax)] = getT(x, jmax);
  }

  auto thermo = solution()->thermo();
  int nsp = thermo->nSpecies();

  // j0 and j1 are constrained to only interior points
  size_t j0 = std::max<size_t>(jmin, 1);
  size_t j1 = std::min(jmax, m_points - 2);
  for (size_t j = j0; j <= j1; j++) {
    double sum = 0.0;
    for (size_t k = 0; k < nsp; k++) {
      double flxk = 0.5 * (m_flux(k, j - 1) + m_flux(k, j));
      double wt = thermo->molecularWeight(k);
      sum += m_wdot(k, j) * m_hk(k, j);
      sum += flxk * dhkdz(x, k, j) / wt;
    }

    rsd[index(Offset::T, j)] =
        -m_cp[j] * getRhoU(x, j) * dTdz(x, j) - divHeatFlux(x, j) - sum;
    rsd[index(Offset::T, j)] /= (m_rho[j] * m_cp[j]);
    rsd[index(Offset::T, j)] -= rdt * (getT(x, j) - getTprev(j));
    diag[index(Offset::T, j)] = 1;
  }
}

void AtmChemistry::evalSpecies(double* x, double* rsd, int* diag, double rdt,
                               size_t jmin, size_t jmax) {
  std::cout << "I'm calling evalSpecies" << std::endl;
  auto thermo = solution()->thermo();
  int nsp = thermo->nSpecies();

  if (jmin == 0) {  // left boundary
    double sum = 0.0;
    for (size_t k = 0; k < nsp; k++) {
      sum += getY(x, k, jmin);
      rsd[index(Offset::Y + k, jmin)] =
          -(m_flux(k, jmin) + getRhoU(x, jmin) * getY(x, k, jmin));
    }
    rsd[index(Offset::Y + m_kExcessLeft, jmin)] = 1.0 - sum;
  }

  if (jmax == m_points - 1) {  // right boundary
    double sum = 0.0;
    for (size_t k = 0; k < nsp; k++) {
      sum += getY(x, k, jmax);
      rsd[index(k + Offset::Y, jmax)] =
          m_flux(k, jmax - 1) + getRhoU(x, jmax) * getY(x, k, jmax);
    }
    rsd[index(Offset::Y + m_kExcessRight, jmax)] = 1.0 - sum;
    diag[index(Offset::Y + m_kExcessRight, jmax)] = 0;
  }

  // j0 and j1 are constrained to only interior points
  size_t j0 = std::max<size_t>(jmin, 1);
  size_t j1 = std::min(jmax, m_points - 2);
  for (size_t j = j0; j <= j1; j++) {
    for (size_t k = 0; k < nsp; k++) {
      double convec = getRhoU(x, j) * dYdz(x, k, j);
      double diffus =
          2.0 * (m_flux(k, j) - m_flux(k, j - 1)) / (z(j + 1) - z(j - 1));

      double wt = thermo->molecularWeight(k);
      rsd[index(Offset::Y + k, j)] =
          (wt * (m_wdot(k, j)) - convec - diffus) / m_rho[j] -
          rdt * (getY(x, k, j) - getYprev(k, j));
      diag[index(Offset::Y + k, j)] = 1;
    }
  }
}

void AtmChemistry::updateThermo(double const* x, size_t j0, size_t j1) {
  std::cout << "I'm calling updateThermo" << std::endl;
  auto thermo = solution()->thermo();

  for (size_t j = j0; j <= j1; j++) {
    thermo->setTemperature(getT(x, j));
    thermo->setPressure(getP(x, j));

    const double* yy = x + m_nv * j + Offset::Y;
    thermo->setMassFractions_NoNorm(yy);

    m_rho[j] = thermo->density();
    m_wtm[j] = thermo->meanMolecularWeight();
    m_cp[j] = thermo->cp_mass();

    thermo->getPartialMolarEnthalpies(&m_hk(0, j));
  }
}

void AtmChemistry::updateReactionRates(double const* x, size_t j0, size_t j1) {
  std::cout << "I'm calling updateReactionRates" << std::endl;
  auto kinetics = solution()->kinetics();

  for (size_t j = j0; j <= j1; j++) {
    kinetics->setActinicFluxLevel(j);
    kinetics->getNetProductionRates(&m_wdot(0, j));
  }
}

void AtmChemistry::updateTransport(double const* x, size_t j0, size_t j1) {
  std::cout << "I'm calling updateTransport" << std::endl;
  auto trans = solution()->transport();
  auto thermo = solution()->thermo();
  int nsp = thermo->nSpecies();

  // mixture averaged transport
  for (size_t j = j0; j < j1; j++) {
    setMidpoint(x, j);

    trans->getMixDiffCoeffs(&m_diff[j * nsp]);

    double rho = thermo->density();
    double wtm = thermo->meanMolecularWeight();
    for (size_t k = 0; k < nsp; k++) {
      double wt = thermo->molecularWeight(k);
      m_diff[k + j * nsp] *= wt * rho / wtm;
    }

    m_tcon[j] = trans->thermalConductivity();
  }
}

void AtmChemistry::updateDiffFluxes(const double* x, size_t j0, size_t j1) {
  std::cout << "I'm calling updateDiffFluxes" << std::endl;
  auto thermo = solution()->thermo();
  int nsp = thermo->nSpecies();

  // mixture averaged transport
  for (size_t j = j0; j < j1; j++) {
    double sum = 0.0;
    double dz = z(j + 1) - z(j);
    for (size_t k = 0; k < nsp; k++) {
      m_flux(k, j) =
          m_diff[k + nsp * j] * (getX(x, k, j) - getX(x, k, j + 1)) / dz;
      sum -= m_flux(k, j);
    }

    // correction flux to insure that \sum_k Y_k V_k = 0.
    for (size_t k = 0; k < nsp; k++) {
      m_flux(k, j) += sum * getY(x, k, j);
    }
  }
}
