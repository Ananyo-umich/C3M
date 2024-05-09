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

  auto iatm = domainIndex("atm");
  m_actinic_flux = std::make_shared<ActinicFlux>();

  m_actinic_flux->setKinetics(domain(iatm).solution()->kinetics());
  m_actinic_flux->setAtmosphere(
      std::dynamic_pointer_cast<AtmChemistry>(m_dom[iatm]));
  // m_actinic_flux->setWavelength(stellar_input.first);
  m_actinic_flux->setWavelength(wavelength);
  // m_actinic_flux->setTOAFlux(stellar_input.second);
  m_actinic_flux->setTOAFlux(actinic_flux);
  m_actinic_flux->initialize();
  setTimeStepCallback(m_actinic_flux.get());

  std::static_pointer_cast<AtmChemistry>(m_dom[iatm])
      ->handleActinicFlux(m_actinic_flux);
}

void AtmChemistrySimulator::resize() {
  std::cout << "total size = " << size() << std::endl;
  OneDim::resize();
  m_xnew.resize(size(), 0.);
}

void AtmChemistrySimulator::setValue(size_t dom, size_t comp, size_t localPoint,
                                     double value) {
  size_t iloc = domain(dom).loc() + domain(dom).index(comp, localPoint);
  AssertThrowMsg(iloc < m_state->size(), "AtmChemistrySimulator::setValue",
                 "Index out of bounds: {} > {}", iloc, m_state->size());
  (*m_state)[iloc] = value;
}

double AtmChemistrySimulator::value(size_t dom, size_t comp,
                                    size_t localPoint) const {
  size_t iloc = domain(dom).loc() + domain(dom).index(comp, localPoint);
  AssertThrowMsg(iloc < m_state->size(), "AtmChemistrySimulator::value",
                 "Index out of bounds: {} > {}", iloc, m_state->size());
  return (*m_state)[iloc];
}

void AtmChemistrySimulator::setProfile(size_t dom, size_t comp,
                                       const std::vector<double>& pos,
                                       const std::vector<double>& values) {
  if (pos.front() != 0.0 || pos.back() != 1.0) {
    throw Cantera::CanteraError(
        "AtmChemistrySimulator::setProfile",
        "`pos` vector must span the range [0, 1]. Got a vector spanning "
        "[{}, {}] instead.",
        pos.front(), pos.back());
  }
  Cantera::Domain1D& d = domain(dom);
  double z0 = d.zmin();
  double z1 = d.zmax();
  for (size_t n = 0; n < d.nPoints(); n++) {
    double zpt = d.z(n);
    double frac = (zpt - z0) / (z1 - z0);
    double v = Cantera::linearInterp(frac, pos, values);
    setValue(dom, comp, n, v);
  }
}

void AtmChemistrySimulator::setFlatProfile(size_t dom, size_t comp, double v) {
  size_t np = domain(dom).nPoints();
  for (size_t n = 0; n < np; n++) {
    setValue(dom, comp, n, v);
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
}

double AtmChemistrySimulator::timeStep(int nsteps, double dt, int loglevel) {
  std::cout << "nsteps = " << nsteps << std::endl;
  std::cout << "dt = " << dt << std::endl;
  std::cout << "loglevel = " << loglevel << std::endl;

  // update atmospheric thermo properties
  auto iatm = domainIndex("atm");
  auto atm = std::dynamic_pointer_cast<AtmChemistry>(m_dom[iatm]);
  atm->updateThermo(m_state->data(), atm->loc(),
                    atm->loc() + atm->nPoints() - 1);

  // update actinic flux
  m_time_step_callback->eval(dt, m_state->data());

  return OneDim::timeStep(nsteps, dt, m_state->data(), m_xnew.data(), loglevel);
}
