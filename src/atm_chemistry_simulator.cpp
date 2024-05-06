// Cantera
#include <cantera/numerics/funcs.h>

// C3M
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

void AtmChemistrySimulator::resize() {
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
