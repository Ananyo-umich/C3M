// cantera
#include <cantera/kinetics/Photolysis.h>

// C3M
#include "actinic_flux.hpp"

void Cantera::Kinetics::handleActinicFlux(
    std::shared_ptr<ActinicFlux> actinic_flux) {
  if (!actinic_flux->initialized()) {
    throw Cantera::CanteraError(
        "Kinetics::handleActinicFlux",
        "ActinicFlux not initialized. Call initialize().");
  }

  m_wavelength = actinic_flux->m_wavelength;
  m_actinicFlux = actinic_flux->m_actinicFlux;
  m_hasNewActinicFlux = true;
}

void ActinicFlux::setFromFile(const std::string &filename,
                              const std::string &format) {
  // Read wavelength and flux from file
  if (format == "kinetics") {
  } else if (format == "json") {
  } else {
  }
}

void ActinicFlux::initialize() {
  if (m_wavelength->empty()) {
    throw Cantera::CanteraError(
        "ActinicFlux::initialize",
        "Wavelength grid is empty. Call setWavelength().");
  }

  size_t nWaves = m_wavelength->size();
  size_t nPoints = m_atm.lock() ? m_atm.lock()->nPoints() : 1;

  if (m_toa_flux.size() != nWaves) {
    throw Cantera::CanteraError(
        "ActinicFlux::initialize",
        "TOA flux size does not match wavelength size.");
  }

  m_actinicFlux = std::make_shared<std::vector<double>>(nPoints * nWaves);
  m_delta_z.resize(nPoints);
  m_dtau.resize(nPoints, nWaves);

  // populate photolysis reaction indices
  m_photo_reactions.clear();
  auto kin = m_kin.lock();
  if (kin) {
    for (int i = 0; i < kin->nReactions(); i++) {
      auto rxn = kin->reaction(i);
      if (rxn->photolysis) m_photo_reactions.push_back(rxn);
    }
  }

  // populate layer thickness
  auto atm = m_atm.lock();
  if (atm) {
    m_delta_z[0] = atm->z(1) - atm->z(0);
    for (int j = 1; j < nPoints - 1; j++) {
      m_delta_z[j] = (atm->z(j + 1) - atm->z(j - 1)) / 2.;
    }
    m_delta_z[nPoints - 1] = atm->z(nPoints - 1) - atm->z(nPoints - 2);
  }

  // populate first TOA actinic flux
  for (size_t k = 0; k < nWaves; k++) {
    m_actinicFlux->at(k) = m_toa_flux[k];
  }

  m_initialized = true;
}

double ActinicFlux::eval(double dt, double *x) {
  if (!m_initialized) {
    throw Cantera::CanteraError(
        "ActinicFlux::eval", "ActinicFlux not initialized. Call initialize().");
  }

  auto atm = m_atm.lock();
  if (!atm) {
    throw Cantera::CanteraError(
        "ActinicFlux::eval",
        "Atmosphere not set in ActinicFlux. Call setAtmosphere().");
  }

  m_dtau.zero();

  size_t levels = atm->nPoints();

  for (size_t j = 0; j < levels; ++j) {
    double temp = atm->getT(x, j);
    double pres = atm->getP(x, j);

    // total number density
    double num_dens = pres / (temp * Cantera::Boltzmann);

    for (auto rxn : m_photo_reactions) {
      size_t n = atm->speciesIndex(rxn->reactants.begin()->first);

      // number density of the parent molecule
      double conc = atm->getX(x, n, j) * num_dens;

      // absorption cross section and optical thickness
      for (size_t k = 0; k < m_wavelength->size(); k++) {
        double wave = m_wavelength->at(k);
        auto &&cross =
            std::static_pointer_cast<Cantera::PhotolysisBase>(rxn->rate())
                ->getCrossSection(temp, wave);
        m_dtau(j, k) += conc * cross[0] * m_delta_z[j];
      }
    }
  }

  // calculate attenuation (simple)
  for (size_t k = 0; k < m_wavelength->size(); k++) {
    double tau = 0.;
    for (size_t j = 0; j < levels; j++) {
      tau += m_dtau(j, k);
      m_actinicFlux->at(j * m_wavelength->size() + k) =
          m_toa_flux[k] * exp(-tau);
    }
  }

  return 0.;
};
