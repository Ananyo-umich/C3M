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

  double h = 6.626e-34; //Planck's constant
  double c = 3e8;  //Speed of light
  double factor = 1/(h*c);
  Eigen::VectorXd wav = Eigen::Map<Eigen::VectorXd>(stellar_input.first.data(), stellar_input.first.size());;
  Eigen::VectorXd irr = Eigen::Map<Eigen::VectorXd>(stellar_input.second.data(), stellar_input.second.size());;
  Eigen::VectorXd actinicFlux = (wav.array()*irr.array()).matrix()*factor;
  // set wavelength
  std::vector<double> wavelength(wav.data(), wav.data() + wav.size());
  std::vector<double> actinic_flux(actinicFlux.data(), actinicFlux.data() + actinicFlux.size());


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
  
  //std::cout << "Starting C3M solve function" << std::endl;	
  // x0 is the previous state
  // x1 is the residual
  eval(Cantera::npos, x0, x1, m_rdt, 1);
  //std::cout << "C3M eval function complete!" << std::endl;
  auto p = left();
  while (p != nullptr) {
  //  std::cout << "flag 1 " << std::endl;
    double* X = x1 + p->loc();
    if (p->isConnector()) {
  //    std::cout << "flag 2" << std::endl;
      p = p->right();
    } else {
  //    std::cout << "flag 3" << std::endl;
      p = p->forwardSweep(X);
    }
  }
//  std::cout << "forward sweep complete" << std::endl;
  p = right();
  while (p != nullptr) {
    double* X = x1 + p->loc();
    if (p->isConnector()) {
      p = p->left();
    } else {
      p = p->backwardSweep(X);
    }
  }

//  std::cout << "backward sweep complete" << std::endl;
  bool good = true;
  for (size_t n = 0; n < size(); n++) {
    if (x0[n] + x1[n] < 0.) {
      good = false;
    }
  }

  if (good) {
    m_successiveSteps++;
    for (size_t n = 0; n < size(); n++) {
      x1[n] += x0[n];
    }
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

  if (m_actinic_flux != nullptr) m_actinic_flux->show();
}

double AtmChemistrySimulator::timeStep(int nsteps, double dt, int loglevel) {
  // update domain cached internal variable
  for (size_t n = 0; n < nDomains(); n++) {
    domain(n).update(m_state->data() + start(n));
  }
  std::cout << "Updated domain data" << std::endl;
  // update actinic flux
  if (m_actinic_flux != nullptr) m_actinic_flux->eval(0., m_state->data());
  if (m_actinic_flux != nullptr) m_actinic_flux->show();
  // if (m_actinic_flux != nullptr) m_actinic_flux->eval(0., m_state->data());

  // update boundary conditions
  eval(Cantera::npos, m_state->data(), m_xnew.data(), m_rdt, 1);
  std::cout << "Updated boundary conditions" << std::endl;
  
  //return AtmChemistrySimulator::solve(m_state->data(), m_xnew.data(), loglevel);
  return OneDim::timeStep(nsteps, dt, m_state->data(), m_xnew.data(), loglevel);
}
