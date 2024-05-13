// application
#include <application/application.hpp>

// Cantera
#include <cantera/numerics/funcs.h>

// C3M
#include "RadTran.hpp"
#include "actinic_flux.hpp"
#include "atm_chemistry_simulator.hpp"

AtmChemistrySimulator::AtmChemistrySimulator(
    std::vector<std::shared_ptr<Cantera::Domain1D>> domains)
    : OneDim(domains) {
  // resize the internal solution vector and the work array, and perform
  // domain-specific initialization of the solution vector.
  resize();
  for (size_t n = 0; n < nDomains(); n++) {
    domain(n)._getInitialSoln(m_state->data() + start(n));
  }
}

void AtmChemistrySimulator::initFromFile(const std::string& filename) {
  auto app = Application::GetInstance();
  auto stellar_input_file = app->FindResource(filename);
  auto stellar_input = ReadStellarRadiationInput(stellar_input_file, 1., 1.);

  // set wavelength
  std::vector<double> wavelength(10);
  std::vector<double> actinic_flux(10);

  for (int i = 0; i < 10; i++) {
    wavelength[i] = 20.0e-9 + i * 20.0e-9;
    actinic_flux[i] = 1.e25;
  }

  auto atm = find<AtmChemistry>("atm");
  m_actinic_flux = std::make_shared<ActinicFlux>();

  m_actinic_flux->setKinetics(atm->solution()->kinetics());
  m_actinic_flux->setAtmosphere(atm);
  // m_actinic_flux->setWavelength(stellar_input.first);
  m_actinic_flux->setWavelength(wavelength);
  // m_actinic_flux->setTOAFlux(stellar_input.second);
  m_actinic_flux->setTOAFlux(actinic_flux);
  m_actinic_flux->initialize();
  setTimeStepCallback(m_actinic_flux.get());

  atm->handleActinicFlux(m_actinic_flux);
}

void AtmChemistrySimulator::resize() {
  OneDim::resize();
  m_xnew.resize(size(), 0.);
}

int AtmChemistrySimulator::solve(double* x0, double* x1, int loglevel) {
  // x0 is the previous state
  // x1 is the residual
  eval(Cantera::npos, x0, x1, m_rdt, 1);

  auto p = left();
  while (p != nullptr) {
    double* X = x1 + p->loc();
    if (p->isConnector()) {
      p = p->right();
    } else {
      p = p->forwardSweep(X);
    }
  }

  p = right();
  while (p != nullptr) {
    double* X = x1 + p->loc();
    if (p->isConnector()) {
      p = p->left();
    } else {
      p = p->backwardSweep(X);
    }
  }

  bool good = true;
  for (size_t n = 0; n < size(); n++) {
    if (x0[n] + x1[n] < 0.) {
      good = false;
    }
  }

  if (good) {
    m_successiveSteps++;
    for (size_t n = 0; n < size(); n++) x1[n] += x0[n];
    if (m_successiveSteps > 3)
      return 100;
    else
      return 1;
  }

  m_successiveSteps = 0;
  return -1;
}

void AtmChemistrySimulator::setValue(std::shared_ptr<Cantera::Domain1D> pdom,
                                     size_t comp, size_t localPoint,
                                     double value) {
  size_t iloc = pdom->loc() + pdom->index(comp, localPoint);
  if (iloc > m_state->size()) {
    throw Cantera::CanteraError("AtmChemistrySimulator::setValue",
                                "Index out of bounds: {} > {}", iloc,
                                m_state->size());
  }
  (*m_state)[iloc] = value;
}

double AtmChemistrySimulator::value(std::shared_ptr<Cantera::Domain1D> pdom,
                                    size_t comp, size_t localPoint) const {
  size_t iloc = pdom->loc() + pdom->index(comp, localPoint);
  if (iloc > m_state->size()) {
    throw Cantera::CanteraError("AtmChemistrySimulator::value",
                                "Index out of bounds: {} > {}", iloc,
                                m_state->size());
  }
  return (*m_state)[iloc];
}

void AtmChemistrySimulator::setProfile(std::shared_ptr<Cantera::Domain1D> pdom,
                                       size_t comp,
                                       const std::vector<double>& pos,
                                       const std::vector<double>& values) {
  if (pos.front() != 0.0 || pos.back() != 1.0) {
    throw Cantera::CanteraError(
        "AtmChemistrySimulator::setProfile",
        "`pos` vector must span the range [0, 1]. Got a vector spanning "
        "[{}, {}] instead.",
        pos.front(), pos.back());
  }

  double z0 = pdom->zmin();
  double z1 = pdom->zmax();
  for (size_t n = 0; n < pdom->nPoints(); n++) {
    double zpt = pdom->z(n);
    double frac = (zpt - z0) / (z1 - z0);
    double v = Cantera::linearInterp(frac, pos, values);
    setValue(pdom, comp, n, v);
  }
}

void AtmChemistrySimulator::setFlatProfile(
    std::shared_ptr<Cantera::Domain1D> pdom, size_t comp, double v) {
  size_t np = pdom->nPoints();
  for (size_t n = 0; n < np; n++) {
    setValue(pdom, comp, n, v);
  }
}

void AtmChemistrySimulator::show() {
  for (size_t n = 0; n < nDomains(); n++) {
    if (domain(n).type() != "empty") {
      Cantera::writelog("\n\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> " +
                        domain(n).id() +
                        " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n");
      domain(n).show(m_state->data() + start(n));
    }
  }

  m_actinic_flux->show();
}

void AtmChemistrySimulator::update() {
  if (m_actinic_flux == nullptr) {
    throw Cantera::CanteraError("AtmChemistrySimulator::update",
                                "Actinic flux not initialized.");
  }

  for (size_t n = 0; n < nDomains(); n++) {
    domain(n).update(m_state->data() + start(n));
  }

  // update actinic flux
  m_actinic_flux->eval(0., m_state->data());
}

double AtmChemistrySimulator::timeStep(int nsteps, double dt, int loglevel) {
  update();

  return OneDim::timeStep(nsteps, dt, m_state->data(), m_xnew.data(), loglevel);
}
