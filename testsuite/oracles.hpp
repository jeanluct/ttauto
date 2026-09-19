// <LICENSE
//   ttauto: a C++ library for building train track automata
//
//   https://github.com/jeanluct/ttauto
//
//   Copyright (C) 2010-2026  Jean-Luc Thiffeault   <jeanluc@math.wisc.edu>
//                            Erwan Lanneau <erwan.lanneau@ujf-grenoble.fr>
//
//   This file is part of ttauto.
//
//   ttauto is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   ttauto is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with ttauto.  If not, see <http://www.gnu.org/licenses/>.
// LICENSE>

// Oracles for the testsuite: slow, independent implementations kept to
// cross-check the library.

#ifndef TTAUTO_TESTSUITE_ORACLES_HPP
#define TTAUTO_TESTSUITE_ORACLES_HPP

#include <cmath>
#include <utility>
#include <vector>
#include <jlt/mathmatrix.hpp>
#include <jlt/polynomial.hpp>
#include "traintracks/traintrack.hpp"

namespace oracles {

// One-fold transition matrix TM(new, old) obtained the original way: put
// unit weight on each old edge in turn, fold a copy, read the new weights.
// Costs n folds; the library now reads the matrix off the fold record.
inline jlt::mathmatrix<int>
unit_weight_transition_matrix(const traintracks::traintrack& tt0, const int f)
{
  const int n = tt0.edges();
  jlt::mathmatrix<int> TM(n,n,0);
  for (int i = 0; i < n; ++i)
    {
      traintracks::traintrack tt(tt0);
      traintracks::traintrack::dblVec wv(n,0.0);
      wv[i] = 1;
      tt.weights(wv.begin());
      tt.fold(f);
      wv = tt.weights();
      for (int j = 0; j < n; ++j) TM(j,i) = (int)wv[j];
    }
  return TM;
}

// Determinant of an integer matrix by fraction-free (Bareiss) elimination.
// Exact in long long for the small matrices of the testsuite.
inline long long bareiss_determinant(std::vector<std::vector<long long> > A)
{
  const int n = A.size();
  if (n == 0) return 1;
  long long sign = 1, prev = 1;
  for (int k = 0; k < n-1; ++k)
    {
      if (A[k][k] == 0)
        {
          int r = k+1;
          while (r < n && A[r][k] == 0) ++r;
          if (r == n) return 0;
          std::swap(A[k],A[r]);
          sign = -sign;
        }
      for (int i = k+1; i < n; ++i)
        for (int j = k+1; j < n; ++j)
          A[i][j] = (A[i][j]*A[k][k] - A[i][k]*A[k][j]) / prev;
      prev = A[k][k];
    }
  return sign*A[n-1][n-1];
}

// Characteristic polynomial det(x I - M) from the sums of principal
// minors: the coefficient of x^(n-k) is (-1)^k times the sum of the k x k
// principal minors.  Independent of the trace-based recursion in
// jlt::mathmatrix::charpoly().  Exponential in n; fine for n <= 12.
inline jlt::polynomial<int> charpoly_by_minors(const jlt::mathmatrix<int>& M)
{
  const int n = M.rows();
  jlt::polynomial<int> p;
  for (int k = 0; k <= n; ++k) p[k] = 0;
  p[n] = 1;
  for (unsigned long mask = 1; mask < (1ul << n); ++mask)
    {
      std::vector<int> idx;
      for (int i = 0; i < n; ++i) if (mask & (1ul << i)) idx.push_back(i);
      const int k = idx.size();
      std::vector<std::vector<long long> > A(k,std::vector<long long>(k));
      for (int i = 0; i < k; ++i)
        for (int j = 0; j < k; ++j) A[i][j] = M(idx[i],idx[j]);
      const long long d = bareiss_determinant(A);
      p[n-k] += (k % 2 == 0 ? 1 : -1) * (int)d;
    }
  return p;
}

// Spectral radius of a nonnegative primitive matrix by power iteration on
// the 1-norm.  Independent of the characteristic polynomial and of the
// Newton solver in ttauto::findroot().
inline double power_iteration_radius(const jlt::mathmatrix<int>& M,
                                     const int itmax = 20000,
                                     const double tol = 1e-13)
{
  const int n = M.rows();
  std::vector<double> v(n,1.0), w(n);
  double lambda = 0;
  for (int it = 0; it < itmax; ++it)
    {
      double norm = 0;
      for (int i = 0; i < n; ++i)
        {
          w[i] = 0;
          for (int j = 0; j < n; ++j) w[i] += M(i,j)*v[j];
          norm += w[i];
        }
      const double lambda_new = norm;   // since sum(v) == 1
      for (int i = 0; i < n; ++i) v[i] = w[i]/norm;
      if (std::abs(lambda_new - lambda) < tol) return lambda_new;
      lambda = lambda_new;
    }
  return lambda;
}

} // namespace oracles

#endif // TTAUTO_TESTSUITE_ORACLES_HPP
