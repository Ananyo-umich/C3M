#ifndef SRC_ATM_CHEMISTRY_SIMULATOR_HPP_
#define SRC_ATM_CHEMISTRY_SIMULATOR_HPP_

// C/C++
#include <memory>
#include <string>
#include <vector>

// cantera
#include <cantera/oneD/Sim1D.h>

class ActinicFlux;

//! AtmChemistrySimulator extends OneDim to enable simulation of atmospheric
//! chemistry over a 1D domain (column).

class AtmChemistrySimulator : public Cantera::OneDim {
 public:
  //! Default constructor.
  AtmChemistrySimulator() {}

  //! Standard constructor.
  AtmChemistrySimulator(
      std::vector<std::shared_ptr<Cantera::Domain1D>> domains);

  //! Destructor.
  ~AtmChemistrySimulator() {}

  //---------------------------------------------------------
  //             overriden general functions
  //---------------------------------------------------------
  //! Resize the solution vector and the work array.
  void resize() override;

  /**
   * Solve F(x) = 0, where F(x) is the multi-domain residual function.
   * @param x0         Starting estimate of solution.
   * @param x1         Final solution satisfying F(x1) = 0.
   * @param loglevel   Controls amount of diagnostic output.
   */
  int solve(double* x0, double* x1, int loglevel) override;

 public:
  //---------------------------------------------------------
  //             special functions
  //---------------------------------------------------------
  void initFromFile(const std::string& filename);

  /**
   * Set a single value in the solution vector.
   * @param dom domain number, beginning with 0 for the leftmost domain.
   * @param comp component number
   * @param localPoint grid point within the domain, beginning with 0 for
   *     the leftmost grid point in the domain.
   * @param value the value.
   */
  void setValue(std::shared_ptr<Cantera::Domain1D> pdom, size_t comp,
                size_t localPoint, double value);

  /**
   * Get one entry in the solution vector.
   * @param dom domain number, beginning with 0 for the leftmost domain.
   * @param comp component number
   * @param localPoint grid point within the domain, beginning with 0 for
   *     the leftmost grid point in the domain.
   */
  double value(std::shared_ptr<Cantera::Domain1D> pdom, size_t comp,
               size_t localPoint) const;

  /**
   * Specify a profile for one component of one domain.
   * @param dom domain number, beginning with 0 for the leftmost domain.
   * @param comp component number
   * @param pos A vector of relative positions, beginning with 0.0 at the
   *     left of the domain, and ending with 1.0 at the right of the domain.
   * @param values A vector of values corresponding to the relative position
   *     locations.
   *
   * Note that the vector pos and values can have lengths different than the
   * number of grid points, but their lengths must be equal. The values at
   * the grid points will be linearly interpolated based on the (pos,
   * values) specification.
   */
  void setProfile(std::shared_ptr<Cantera::Domain1D> pdom, size_t comp,
                  const std::vector<double>& pos,
                  const std::vector<double>& values);

  //! Set component 'comp' of domain 'dom' to value 'v' at all points.
  void setFlatProfile(std::shared_ptr<Cantera::Domain1D> pdom, size_t comp,
                      double v);

  //! Show logging information on current solution for all domains.
  void show();

  /**
   * Take time steps using Backward Euler.
   *
   * @param nsteps number of steps
   * @param dt initial step size
   * @param loglevel controls amount of printed diagnostics [0-8]
   * @returns size of last timestep taken
   */
  double timeStep(int nsteps, double dt, int loglevel);

  //! downcast function
  //! @return a shared pointer to the domain of type T
  template <typename T = Cantera::Domain1D>
  std::shared_ptr<T> find(std::string const& name) const {
    size_t index = domainIndex(name);
    return std::dynamic_pointer_cast<T>(m_dom[index]);
  }

  //! calls the update function of each domain
  void update();

 protected:
  //! a work array used to hold the residual or the new solution
  std::vector<double> m_xnew;

  //! actinic flux handler
  std::shared_ptr<ActinicFlux> m_actinic_flux;

  //! number of successive steps
  size_t m_successiveSteps = 0;
};

#endif  // SRC_ATM_CHEMISTRY_SIMULATOR_HPP_:
