
// Author: Ananyo Bhattacharya
// Affiliation: University of Michigan
// Email: ananyo@umich.edu

// C/C++ headers
#include <cantera/base/Solution.h>
#include <cantera/kinetics.h>
#include <cantera/kinetics/Reaction.h>
#include <cantera/numerics/eigen_dense.h>
#include <cantera/numerics/eigen_sparse.h>
#include <cantera/thermo.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

void DiffusiveFlux(Eigen::MatrixXd Atm, Eigen::MatrixXd Xf,
                   Eigen::MatrixXd flux, Eigen::VectorXd alpha, double grav,
                   std::string PlanetName, std::string FileName) {
  // Extracting input file details
  auto sol = Cantera::newSolution(FileName);
  auto gas = sol->thermo();
  Cantera::ThermoPhase *gasThermo = gas.get();

  // Unpack matrix into atmospheric data: temperature, pressure, number density
  int iTemp = 0;
  int iPress = 1;
  int iKzz = 2;
  int iAlt = 3;
  int iNd = 4;

  // Unpack matrix into atmospheric compostion:
  int nSize = Xf.row(0).size();
  int nsp = Xf.col(0).size();
  double Nt;
  Eigen::VectorXd mWt(nsp);  // Molecular weight

  for (int insp = 0; insp < nsp; insp++) {
    mWt(insp) = gas->molecularWeight(insp);  // Kg/kml
  }

  // Initialize diffusion variables
  double K_next, K_prev, dz_prev, dz_next;
  Eigen::VectorXd I = Eigen::VectorXd::Ones(nsp);
  // Loop over all the rows to get diffusive flux at mid point of grid values
  int j = 1;
  while ((j > 0) & (j < nSize)) {
    double T_next = (Atm(iTemp, j) + Atm(iTemp, j + 1)) / 2;
    double T_prev = (Atm(iTemp, j) + Atm(iTemp, j - 1)) / 2;
    double T_this = Atm(iTemp, j);

    gas->setMoleFractions(&Xf.col(j - 1)[0]);
    double mm_prev = gas->meanMolecularWeight();  // Mean molcular weight

    gas->setMoleFractions(&Xf.col(j + 1)[0]);
    double mm_next = gas->meanMolecularWeight();  // Mean molcular weight

    // Eddy diffusion coefficient
    K_next = (Atm(iKzz, j + 1) + Atm(iKzz, j)) / 2;
    K_prev = (Atm(iKzz, j - 1) + Atm(iKzz, j)) / 2;
    // Altitude grid step size
    dz_prev = (Atm(iAlt, j - 1) - Atm(iAlt, j));
    dz_next = (Atm(iAlt, j) - Atm(iAlt, j + 1));
    // Molecular diffusion coefficient
    Eigen::VectorXd D_this = Eigen::VectorXd::Zero(
        nsp);  // handleCustomMolecularDiffusion(PlanetName, gasThermo,
               // Atm(iPress,j), Atm(iTemp,j), mWt);
    Eigen::VectorXd D_next = Eigen::VectorXd::Zero(
        nsp);  // handleCustomMolecularDiffusion(PlanetName, gasThermo,
               // Atm(iPress,j+1), Atm(iTemp,j+1), mWt);
    Eigen::VectorXd D_prev = Eigen::VectorXd::Zero(
        nsp);  // handleCustomMolecularDiffusion(PlanetName, gasThermo,
               // Atm(iPress,j-1), Atm(iTemp,j-1), mWt);
    Eigen::VectorXd d_next =
        ((mm_next * I - mWt).array() * D_next.array()).matrix() *
        (1 / (2 * dz_next * dz_next)) *
        (grav * 1e-3 * dz_next / (6.022E23 * 1.38E-23 * T_next));
    Eigen::VectorXd d_prev =
        ((mm_prev * I - mWt).array() * D_prev.array()).matrix() *
        (1 / (2 * dz_prev * dz_prev)) *
        (grav * 1e-3 * dz_prev / (6.022E23 * 1.38E-23 * T_prev));
    // Thermal diffusion
    Eigen::VectorXd dTdz_prev =
        alpha * (T_this - T_prev) * 2 / ((T_this + T_prev) * dz_prev);
    Eigen::VectorXd dTdz_next =
        alpha * (T_next - T_this) * 2 / ((T_this + T_next) * dz_next);
    // Total contribution due to molecular diffusion and eddy diffusion
    d_next = d_next + ((dTdz_next.array() * D_next.array()).matrix() *
                       (1 / (2 * dz_next * dz_next)));
    d_prev = d_prev + ((dTdz_prev.array() * D_prev.array()).matrix() *
                       (1 / (2 * dz_prev * dz_prev)));
    Eigen::VectorXd k_next = (((K_next * I) + D_next) / (dz_next * dz_next));
    Eigen::VectorXd k_prev = (((K_prev * I) + D_prev) / (dz_prev * dz_prev));

    double N_next = Atm(iNd, j + 1);
    double N_this = Atm(iNd, j);
    double N_prev = Atm(iNd, j - 1);
    double N_n = (N_this + N_next) / 2;
    double N_p = (N_this + N_prev) / 2;

    // Diffusive flux at grid point (other than boundary points)
    // Contribution of U_i
    Eigen::VectorXd B = -1 * ((k_next * N_n) + (d_next * N_n) + (k_prev * N_p) -
                              (d_prev * N_p));
    // Contribution of U_i-1
    Eigen::VectorXd A = ((k_prev * N_p) + (d_prev * N_p));
    // Contribution of U_i+1
    Eigen::VectorXd C = ((k_next * N_n) - (d_next * N_n));

    // Total flux at grid point
    flux.col(j) =
        (A.array() * Xf.col(j - 1).array() + B.array() * Xf.col(j).array() +
         C.array() * Xf.col(j + 1).array())
            .matrix();

    // Update the index
    j++;
  }
}
