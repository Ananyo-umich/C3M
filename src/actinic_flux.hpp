#ifndef SRC_ACTINIC_FLUX_HPP_
#define SRC_ACTINIC_FLUX_HPP_

// C/C++
#include <memory>
#include <string>

// cantera
#include <cantera/kinetics/Kinetics.h>
#include <cantera/kinetics/Reaction.h>
#include <cantera/numerics/Func1.h>
#include <cantera/oneD/Domain1D.h>

// C3M
#include "atm_chemistry.hpp"

class ActinicFlux : public Cantera::Func1 {
  friend class Cantera::Kinetics;

 public:
  ActinicFlux() = default;

  ActinicFlux(const std::vector<double>& wavelength,
              const std::vector<double>& flux,
              std::shared_ptr<Cantera::Kinetics> kin,
              std::shared_ptr<AtmChemistry> atm) {
    setKinetics(kin);
    setAtmosphere(atm);
    setWavelength(wavelength);
    setTOAFlux(flux);
    initialize();
  }

  //! set Cantera Kinetics object
  void setKinetics(std::shared_ptr<Cantera::Kinetics> kin) { m_kin = kin; }

  //! set Cantera Domain1D object
  void setAtmosphere(std::shared_ptr<AtmChemistry> atm) { m_atm = atm; }

  //! set TOA irradiance [w/(m^2.m)]
  void setTOAFlux(std::vector<double> const& flux) { m_toa_flux = flux; }

  //! set wavelength grid [m]
  void setWavelength(std::vector<double> const& wavelength) {
    m_wavelength = std::make_shared<std::vector<double>>(wavelength);
  }

  void setFromFile(const std::string& filename, const std::string& format);

  //! set cosine of stellar zenith angle
  void setStellarZenithAngle(double theta) { m_stelar_mu = cos(theta); }

  //! resize internal data structures based on wavelength and atmosphere
  //! dimensions
  void initialize();

  //! check if object is initialized
  bool initialized() const { return m_initialized; }

  std::string type() const override { return "actinic_flux"; }

  //! call back function after each successful time step
  //! @param dt time step [s]
  //! @param x global state vector
  double eval(double dt, double* x) override;

 protected:
  //! indicate if object is initialized
  bool m_initialized = false;

  //! cosine of stellar zenith angle
  double m_stelar_mu;

  //! Cantera Kinetics object
  std::weak_ptr<Cantera::Kinetics> m_kin;

  //! Cantera Domain1D object
  std::weak_ptr<AtmChemistry> m_atm;

  //! TOA irradiance [photons/(s.m^2.m)] (nWaves)
  std::vector<double> m_toa_flux;

  //! layer thickness [m] (nPoints)
  std::vector<double> m_delta_z;

  //! wavelength grid [m] (nWaves), shared with Kinetics
  std::shared_ptr<std::vector<double>> m_wavelength;

  //! actinic flux (nPoints x nWaves), shared with Kinetics
  std::shared_ptr<std::vector<double>> m_actinicFlux;

  //! photolysis reactions
  std::vector<std::shared_ptr<Cantera::Reaction>> m_photo_reactions;

  //! optical thickness [1] (nPoints x nWaves)
  Cantera::Array2D m_dtau;
};

#endif  // SRC_ATM_ACTINIC_FLUX_HPP_
