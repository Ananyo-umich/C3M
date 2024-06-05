#ifndef SRC_ATM_CHEMISTRY_HPP
#define SRC_ATM_CHEMISTRY_HPP

// C/C++
#include <memory>
#include <string>

// Eigen
#include <Eigen/Core>
#include <Eigen/Sparse>

// cantera
#include <cantera/base/Array.h>
#include <cantera/base/Solution.h>
#include <cantera/kinetics/Kinetics.h>
#include <cantera/oneD/Domain1D.h>
#include <cantera/thermo/ThermoPhase.h>

class ActinicFlux;

//! Domain1D is a lightweight class that holds information about a 1D domain.
//! It does not contain the solution vector itself, but only information about
//! the domain, such as the number of grid points and the number of components
//! (strides).
//!
//! Domain1D has its own a Solution object, accessed by solution(), which is a
//! manger class for Thermo, Kinetics, Transport, etc.
//!
//! There are two kinds of domains, Bulk and Connector. Bulk is a domain that
//! contains interior points, while Connector is a domain that contains boundary
//! points. Both Bulk and Connector are derived from Domain1D and have their own
//! implementations of physical functions.
//!
//! Bulk and Connector arranged as the following:
//!
//!             .----.----.----.---- Bulk 1 managed interior points
//!             |    |    |    |
//!        .----+----+----+----+-----.-------------------------.
//!        |                         |                         |
//! (Connector 0) -- (Bulk 1) -- (Connector 1) -- (Bulk 2) -- (Connector 2)
//!        |                         |                         |
//!        * -- + -- + -- + -- + --  * -- + -- + -- + -- + --  *
//!        ^                         ^                         ^
//!        |-- boundary point        |-- boundary point        |-- boundary
//!
//!
//! Each Bulk object has a left and right boundary point (*). The boundary point
//! is a point where the solution is managed by a Connector.
//! The Bulk object manages the interior points (+) of the Domain1D object.

//! \brief AtmChemistry extends Domain1D to represent a column of atmosphere.
class AtmChemistry : public Cantera::Domain1D {
  //! Offsets of solution components in the 1D solution array.
  //! temperature, pressure, velocity, mass fractions
  enum Offset { T, P, U };

 public:
  AtmChemistry(std::string const& id, std::shared_ptr<Cantera::Solution> sol,
               size_t npoints = 1, size_t stride = Cantera::npos);

  ~AtmChemistry() {}

  //---------------------------------------------------------
  //             overriden general functions
  //---------------------------------------------------------
  //! concrete type of the domain
  std::string domainType() const override;

  //! reset bad values in the global solution vector
  //! @param xg global mole fraction solution vector
  void resetBadValues(double* xg) override;

  //! @return the size of the solution vector
  size_t size() const override;

  //! change the grid size. Called after grid refinement.
  void resize(size_t stride, size_t points) override;

  //! called to set up initial grid, and after grid refinement
  void setupGrid(size_t n, const double* z) override;

  //! Initial value of solution component @e n at grid point @e j.
  double initialValue(size_t n, size_t j) override;

  //! Name of the nth component.
  std::string componentName(size_t n) const override;

  //! Update the following thermo properties at levels
  //! (1) mean molecular weight, m_wtm
  //! (2) molar enthalpy, m_hk
  //! Update the following thermo properties at half-levels
  //! (1) temperature, m_Tmid
  //! (2) pressure, m_Pmid
  //! (3) composition, m_Xmid
  //! (4) molecular weight anomaly, m_dwtm
  //!
  //! @param x mole fraction vector
  void update(double const* x) override;

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
   *  @param[in] yGlobal  Global state vector
   *  @param[out] rsdGlobal  Global residual vector
   *  @param[out] diagGlobal  Global boolean mask indicating whether each
   * solution component has a time derivative (1) or not (0).
   *  @param[in] rdt Reciprocal of the timestep (`rdt=0` implies steady-state.)
   */
  void eval(size_t jGlobal, double* yGlobal, double* rsdGlobal,
            integer* diagGlobal, double rdt) override;

 public:
  //---------------------------------------------------------
  //             special functions
  //---------------------------------------------------------
  //! Returns true if the specified component is an active part of the solver
  //! state
  virtual bool componentActive(size_t n) const { return m_do_species[n]; }

  void setGraivty(double g) { m_grav = g; }

  //! @return temperature at grid point `j`
  double getT(size_t j) const {
    auto tmp = m_hydro.lock();
    if (!tmp) {
      throw Cantera::CanteraError("AtmChemistry::getT",
                                  "Hydro not initialized. Call setHydro().");
    }
    return tmp->at(index_hydro(Offset::T, j));
  }

  //! @return pressure at grid point `j`
  double getP(size_t j) const {
    auto tmp = m_hydro.lock();
    if (!tmp) {
      throw Cantera::CanteraError("AtmChemistry::getP",
                                  "Hydro not initialized. Call setHydro().");
    }
    return tmp->at(index_hydro(Offset::P, j));
  }

  //! @return velocity at grid point `j`
  double getU(size_t j) const {
    auto tmp = m_hydro.lock();
    if (!tmp) {
      throw Cantera::CanteraError("AtmChemistry::getP",
                                  "Hydro not initialized. Call setHydro().");
    }
    return tmp->at(index_hydro(Offset::U, j));
  }

  //! @return mass fraction of species `k` at grid point `j`
  //! @param x mole fraction array
  double getY(const double* x, size_t k, size_t j) const {
    double wt = solution()->thermo()->molecularWeight(k);
    return getX(x, k, j) * wt / m_wtm[j];
  }

  //! @return mole fraction of species `k` at grid point `j`
  //! @param x mole fraction array
  double getX(const double* x, size_t k, size_t j) const {
    return x[index(k, j)];
  }

  //! @return layer thickness
  double getH(size_t j) const { return m_delta_z[j]; }

 public:
  //---------------------------------------------------------
  //             hydro coupler
  //---------------------------------------------------------
  //! set const ptr to hydro dara array with stride
  void setHydro(std::shared_ptr<std::vector<double>> hydro, size_t stride) {
    m_hydro = hydro;
    m_hydro_stride = stride;
  }

  //---------------------------------------------------------
  //             kinetics coupler (actinic flux)
  //---------------------------------------------------------
  //! Pass actinic flux to kinetics object
  void handleActinicFlux(std::shared_ptr<ActinicFlux> actinic_flux) {
    solution()->kinetics()->handleActinicFlux(actinic_flux);
  }
  
  //! @set the Eddy diffusion coefficient 'Kzz' at grid point `j`
  //! @param Kzz value of Eddy diffusion coefficient
  virtual void setEddyDiffusionCoeff(double Kzz, size_t j);

  //! @return Eddy diffusion coefficient at grid point `j`
  virtual double getEddyDiffusionCoeff(size_t j);

  //! @set matrix of binary diffusion coefficients at grid point `j`
  //! @param BinDiff Binary diffusion matrix
  virtual void setBinaryDiffusionCoeff(Eigen::MatrixXd BinDiff, size_t j);

  virtual Eigen::MatrixXd getBinaryDiffusionCoeff(size_t j) const {
    size_t nsp = nSpecies();
    Eigen::MatrixXd Kbinary(nsp, nsp);
    Kbinary.setZero();
    return Kbinary;
  }

  virtual void setGravity(double grav);

 protected:
  //---------------------------------------------------------
  //             physical functions
  //---------------------------------------------------------
  //! update #m_wdot (species production rates)
  virtual void updateReaction(double rdt, double const* x, size_t j0,
                              size_t j1);

  //! Update the convection properties at grid points in the range from `j0`
  //! to `j1`, based on solution `x`.
  virtual void updateConvection(double const* x, size_t j0, size_t j1);

  //! Update the diffusive mass fluxes.
  virtual void updateDiffusion(const double* x, size_t j0, size_t j1);


 protected:
  //---------------------------------------------------------
  //             helper functions
  //---------------------------------------------------------
  //! @return number of variables in a layer
  size_t stride() const { return m_nv; }

  size_t index_hydro(size_t n, size_t j) const {
    return m_hydro_stride * j + n;
  }

  Eigen::VectorXd dXdz(const double* x, size_t j) const {
    size_t nsp = nSpecies();
    Eigen::VectorXd dXdz(nsp);
    for (size_t k = 0; k < nsp; ++k) {
      dXdz(k) = (getX(x, k, j) - getX(x, k, j - 1)) / (m_z[j] - m_z[j - 1]);
    }
    return dXdz;
  }

  //! TODO(cli) move to transport

 protected:
  //---------------------------------------------------------
  //             member data
  //---------------------------------------------------------
  //! cached gravity
  double m_grav = 0.0;

  //! active species indicator (nSpecies)
  std::vector<bool> m_do_species;

  //! layer thickness [m] (nPoints)
  std::vector<double> m_delta_z;

  //! cached mean molecular weight [kg/mol] (nPoints)
  std::vector<double> m_wtm;

  //! cached species molar enthalpies [J/mol] (nSpecies x nPoints)
  Cantera::Array2D m_hk;

  //! cached species diffusion flux (nSpecies x nPoints)
  Cantera::Array2D m_diff_flux;

  //! cached species diffusion tendency (nSpecies x nPoints)
  Cantera::Array2D m_diff_t;

  //! cached eddy diffusion coefficients (nPoints)
  std::vector<double> m_Keddy;

  //! cached binary diffusion coefficients (nPoints)
  std::vector<Eigen::MatrixXd> m_Kbinary;

  //! cached species production rates (nSpecies x nPoints)
  Cantera::Array2D m_wdot;

  //! cached species convection tendency (nSpecies x nPoints)
  Cantera::Array2D m_conv_t;

 private:
  //---------------------------------------------------------
  //             half grid member data
  //---------------------------------------------------------

  //! cached species mole fraction at half grid (nPoints)
  std::vector<double> m_Xmid;

  //! cached temperature at half grid (nPoints)
  std::vector<double> m_Tmid;

  //! cached pressure at half grid (nPoints)
  std::vector<double> m_Pmid;

  //! cached molecular weight difference at half grid (nPoints)
  std::vector<Eigen::VectorXd> m_dwtm;

  //! cached pointer to hydro
  std::weak_ptr<std::vector<double> const> m_hydro;

  //! stride between two levels in hydro
  size_t m_hydro_stride = 0;
};

#endif  // SRC_ATM_CHEMISTRY_HPP
