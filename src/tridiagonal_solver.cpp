// C/C++ headers
#include <iostream>
#include <vector>

// Eigen
#include <Eigen/Sparse>

// Cantera
#include <cantera/oneD/Domain1D.h>

Cantera::Domain1D *Cantera::Domain1D::forwardSweep(double const *residual) {
  auto &a = m_B;
  auto &b = m_A;
  auto &c = m_C;
  auto &delta = m_D;

  size_t il = isConnector() ? 0 : 1;
  size_t iu = m_points - 2;
  size_t nv = delta[0].size();

  Eigen::VectorXd rhs(nv);
  for (int n = 0; n < nv; ++n) rhs(n) = residual[n];

  Eigen::SparseMatrix<double> identity(nv, nv);
  for (int k = 0; k < nv; ++k) identity.insert(k, k) = 1.;

  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

  if (m_left) {
    m_left->modifyRight(a[il - 1], delta[il - 1]);
    // a[il] = (a[il] - b[il] * a[il - 1]).inverse().eval();
    solver.compute(a[il] - b[il] * a[il - 1]);
    a[il] = solver.solve(identity);
    delta[il] = a[il] * (rhs - b[il] * delta[il - 1]);
    a[il] = a[il] * c[il];
  } else {
    // a[il] = a[il].inverse().eval();
    solver.compute(a[il]);
    a[il] = solver.solve(identity);
    delta[il] = a[il] * rhs;
    a[il] = a[il] * c[il];
  }

  for (int i = il + 1; i <= iu; ++i) {
    for (int n = 0; n < nv; ++n) rhs(n) = residual[n + i * m_points];
    // a[i] = (a[i] - b[i] * a[i - 1]).inverse().eval();
    solver.compute(a[i] - b[i] * a[i - 1]);
    a[i] = solver.solve(identity);
    delta[i] = a[i] * (rhs - b[i] * delta[i - 1]);
    a[i] = a[i] * c[i];
  }

  return m_right;
}

Cantera::Domain1D *Cantera::Domain1D::backwardSweep(double *result) {
  auto &a = m_B;
  auto &delta = m_D;

  size_t il = 1;
  size_t iu = m_points - 2;
  size_t nv = delta[0].size();

  if (m_right) {
    m_right->modifyLeft(delta[iu + 1]);
    delta[iu] -= a[iu] * delta[iu + 1];
  }

  for (int i = iu - 1; i >= il; --i) {
    delta[i] -= a[i] * delta[i + 1];
  }

  for (int i = il; i <= iu; ++i) {
    for (int n = 0; n < nv; ++n) result[n + i * nv] = delta[i](n);
  }

  return m_left;
}
