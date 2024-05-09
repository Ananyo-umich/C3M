#ifndef SRC_ATM_CHEMISTRY_HPP
#define SRC_ATM_CHEMISTRY_HPP

// C/C++
#include <memory>
#include <string>

// cantera
#include <cantera/base/Array.h>
#include <cantera/base/Solution.h>
#include <cantera/oneD/Domain1D.h>
#include <cantera/thermo/ThermoPhase.h>

class ActinicFlux;

//! Domain1D is a lightweight class that holds information about a 1D domain.
//! It does not contain the solution vector itself, but only information about
//! the domain, such as the number of grid points and the number of components.
//! Domain1D holds a pointer to the solution vector (m_state), which is owned by
//! OneDim. Domain1D has its own a Solution object (m_solution), which is a
//! manger class for Thermo, Kinetics, etc.

//! \brief AtmChemistry extends Domain1D to represent a column of atmosphere.
class AtmChemistry : public Cantera::Domain1D {
  //! Offsets of solution components in the 1D solution array.
  //! temperature, pressure, velocity, mass fractions
  enum Offset { T, P, U, Y };

 public:
  AtmChemistry(std::string const& id, std::shared_ptr<Cantera::Solution> sol,
               size_t npoints = 1);

  ~AtmChemistry() {}

  std::string domainType() const override;

  //! reset bad values in the solution vector
  void resetBadValues(double* xg) override;

  //! Change the grid size. Called after grid refinement.
  void resize(size_t components, size_t points) override;

  //! Set the state to be consistent with the solution at the midpoint
  //! between j and j + 1.
  void setMidpoint(const double* x, size_t j);

  //! Initial value of solution component @e n at grid point @e j.
  double initialValue(size_t n, size_t j) override { return 0.; }

  //! Name of the nth component.
  std::string componentName(size_t n) const override;

  //! Index of the species in the solution vector.
  size_t speciesIndex(const std::string& name) const {
    return componentIndex(name) - Offset::Y;
  }

  //! Pass actinic flux to kinetics object
  void handleActinicFlux(std::shared_ptr<ActinicFlux> actinic_flux) {
    solution()->kinetics()->handleActinicFlux(actinic_flux);
  }

  //! Returns true if the specified component is an active part of the solver
  //! state
  virtual bool componentActive(size_t n) const { return m_do_species[n]; }

  /**
   * Evaluate the residual functions for atmospheric chemistry
   * If jGlobal == npos, the residual function is evaluated at all grid points.
   * Otherwise, the residual function is only evaluated at grid points j-1, j,
   * and j+1. This option is used to efficiently evaluate the Jacobian
   * numerically.
   *
   * These residuals at all the boundary grid points are evaluated using a
   * default boundary condition that may be modified by a boundary object that
   * is attached to the domain. The boundary object connected will modify these
   * equations by subtracting the boundary object's values for V, T, mdot, etc.
   * As a result, these residual equations will force the solution variables to
   * the values of the connected boundary object.
   *
   *  @param jGlobal  Global grid point at which to update the residual
   *  @param[in] xGlobal  Global state vector
   *  @param[out] rsdGlobal  Global residual vector
   *  @param[out] diagGlobal  Global boolean mask indicating whether each
   * solution component has a time derivative (1) or not (0).
   *  @param[in] rdt Reciprocal of the timestep (`rdt=0` implies steady-state.)
   */
  void eval(size_t jGlobal, double* xGlobal, double* rsdGlobal,
            integer* diagGlobal, double rdt) override;

  //! @return temperature at grid point `j`
  double getT(double const* x, size_t j) const {
    return x[index(Offset::T, j)];
  }

  //! @return pressure at grid point `j`
  double getP(double const* x, size_t j) const {
    return x[index(Offset::P, j)];
  }

  //! @return mass fraction of species `k` at grid point `j`
  double getY(const double* x, size_t k, size_t j) const {
    return x[index(Offset::Y + k, j)];
  }

  //! @return mole fraction of species `k` at grid point `j`
  double getX(const double* x, size_t k, size_t j) const {
    double wt = solution()->thermo()->molecularWeight(k);
    return m_wtm[j] * getY(x, k, j) / wt;
  }

  /**
   * Update the thermodynamic properties from point j0 to point j1
   * (inclusive), based on solution x.
   *
   * The gas state is set to be consistent with the solution at the
   * points from j0 to j1.
   *
   * Properties that are computed and cached are:
   * * #m_rho (density)
   * * #m_wtm (mean molecular weight)
   * * #m_cp (specific heat capacity)
   * * #m_hk (species specific enthalpies)
   */
  virtual void updateThermo(double const* x, size_t j0, size_t j1);

 protected:
  //---------------------------------------------------------
  //             physical functions
  //---------------------------------------------------------

  virtual void updateProperties(size_t jg, double* x, size_t jmin, size_t jmax);
  virtual void evalContinuity(double* x, double* rsd, int* diag, double rdt,
                              size_t jmin, size_t jmax);
  virtual void evalEnergy(double* x, double* rsd, int* diag, double rdt,
                          size_t jmin, size_t jmax);

  //! update #m_wdot (species production rates)
  virtual void updateReactionRates(double const* x, size_t j0, size_t j1);

  //! Update the transport properties at grid points in the range from `j0`
  //! to `j1`, based on solution `x`.
  virtual void updateTransport(double const* x, size_t j0, size_t j1);

  //! Update the diffusive mass fluxes.
  virtual void updateDiffFluxes(const double* x, size_t j0, size_t j1);

  //! @return previous temperature at grid point `j`
  double getTprev(size_t j) const { return prevSoln(Offset::T, j); }

  //! @return The fixed temperature value at point j.
  double getTfixed(size_t j) const { return m_fixedtemp[j]; }

  //! @return velocity at grid point `j`
  double getU(const double* x, size_t j) const {
    return x[index(Offset::U, j)];
  }

  //! @return density at grid point `j`
  double getRho(size_t j) const { return m_rho[j]; }

  //! @return density multiplied by velocity at grid point `j`
  double getRhoU(double const* x, size_t j) const {
    return m_rho[j] * x[index(Offset::U, j)];
  }

  //! @return previous mass fraction of species `k` at grid point `j`
  double getYprev(size_t k, size_t j) const {
    return prevSoln(Offset::Y + k, j);
  }

  //! @return upwind mass fraction derivative at grid point `j` of species `k`
  double dYdz(const double* x, size_t k, size_t j) const {
    size_t jloc = (getU(x, j) > 0.0 ? j : j + 1);
    double dz = m_z[jloc] - m_z[jloc - 1];
    return (getY(x, k, jloc) - getY(x, k, jloc - 1)) / dz;
  }

  //! @return upwind temperature derivative at grid point `j`
  double dTdz(const double* x, size_t j) const {
    size_t jloc = (getU(x, j) > 0.0 ? j : j + 1);
    double dz = m_z[jloc] - m_z[jloc - 1];
    return (getT(x, jloc) - getT(x, jloc - 1)) / dz;
  }

  //! @return upwind enthalpy derivative at grid point `j` of species `k`
  double dhkdz(const double* x, size_t k, size_t j) const {
    if (getU(x, j) > 0.0) {
      double dz = m_z[j] - m_z[j - 1];
      return (m_hk(k, j) - m_hk(k, j - 1)) / dz;
    } else {
      double dz = m_z[j + 1] - m_z[j];
      return (m_hk(k, j + 1) - m_hk(k, j)) / dz;
    }
  }

  //! @return conductive heat divergence at grid point `j`
  double divHeatFlux(const double* x, size_t j) const {
    double c1 = m_tcon[j - 1] * (getT(x, j) - getT(x, j - 1));
    double c2 = m_tcon[j] * (getT(x, j + 1) - getT(x, j));
    return -2.0 * (c2 / (z(j + 1) - z(j)) - c1 / (z(j) - z(j - 1))) /
           (z(j + 1) - z(j - 1));
  }

  /**
   * Evaluate the species equations' residuals.
   *
   * The function calculates the species equations as
   * @f[
   *    \rho u \frac{dY_k}{dz} + \frac{dj_k}{dz} = W_k\dot{\omega}_k
   * @f]
   *
   * The species equations include terms for temporal and spatial variations
   * of species mass fractions (@f$ Y_k @f$). The default boundary condition
   * is zero flux for species at the left and right boundary.
   *
   * For argument explanation, see evalContinuity().
   */
  virtual void evalSpecies(double* x, double* rsd, int* diag, double rdt,
                           size_t jmin, size_t jmax);

  //---------------------------------------------------------
  //             member data
  //---------------------------------------------------------
  //! active component indicator (nSpecies)
  std::vector<bool> m_do_species;

  //! cached density (nPoints)
  std::vector<double> m_rho;

  //! cached mean molecular weight (nPoints)
  std::vector<double> m_wtm;

  //! cached specific heat capacity (nPoints)
  std::vector<double> m_cp;

  //! cached thermal conductivity (nPoints)
  std::vector<double> m_tcon;

  //! cached species diffusion coefficients (nSpecies x nPoints)
  std::vector<double> m_diff;

  // fixed T
  std::vector<double> m_fixedtemp;

  //! cached species specific enthalpies (nSpecies x nPoints)
  Cantera::Array2D m_hk;

  //! cached species production rates (nSpecies x nPoints)
  Cantera::Array2D m_wdot;

  //! cached species diffusive mass fluxes (nSpecies x nPoints)
  Cantera::Array2D m_flux;

  //! Index of species with a large mass fraction at each boundary, for which
  //! the mass fraction may be calculated as 1 minus the sum of the other mass
  //! fractions
  size_t m_kExcessLeft = 0;
  size_t m_kExcessRight = 0;

  //! energy equation indicator
  bool m_do_energy;

 private:
  std::vector<double> m_ymid;
};

#endif  // SRC_ATM_CHEMISTRY_HPP
