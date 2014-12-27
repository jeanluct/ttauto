// <LICENSE
//   ttauto: a C++ library for building train track automata
//
//   https://github.com/jeanluct/ttauto
//
//   Copyright (C) 2010-2014  Jean-Luc Thiffeault   <jeanluc@math.wisc.edu>
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

#ifndef TRAINTRACKS_UTIL_HPP
#define TRAINTRACKS_UTIL_HPP

#include <jlt/vector.hpp>
#include <jlt/mathmatrix.hpp>
#include <jlt/polynomial.hpp>
#include <jlt/exceptions.hpp>

namespace ttauto {

// If two vectors are cyclically equivalent, return a vector p0v of
// offsets between them such that v1[v] == v2[(v+p0v[i]) % size()].
// Return empty p0v if they are not cyclically equivalent.
template<class Vec>
jlt::vector<int> cyclic_shift(const Vec& v1, const Vec& v2)
{
  if (v1.size() != v2.size()) return jlt::vector<int>();

  jlt::vector<int> p0v;

  const int k = v1.size();

  // Shift the starting point of the second path to test for 'cyclic
  // equality'.
  // Try all possibilities, and return a list of matching offsets.
  for (int p0 = 0; p0 < k; ++p0)
    {
      int p = 0;
      for (; p < k; ++p)
	{
	  int p1 = (p+p0) % k;
	  if (!(v1[p] == v2[p1])) break;
	}
      // We made it through without breaking out, so must tbe equal.
      if (p == k) p0v.push_back(p0);
    }

  return p0v;
}


template<class Vec>
inline bool cyclic_equal(const Vec& v1, const Vec& v2)
{
  return !(cyclic_shift(v1,v2).empty());
}


template<class Mat>
inline bool pattern_equal(const Mat& A, const Mat& B)
{
  typename Mat::const_iterator i;
  typename Mat::const_iterator j = B.begin();

  for (i = A.begin(); i != A.end(); ++i, ++j)
    {
      if ((*i != 0 && *j == 0) || (*i == 0 && *j != 0)) return false;
    }

  return true;
}


template<class TrTr>
jlt::mathmatrix<int> fold_transition_matrix(const TrTr& tt0, const int f)
{
  // Chaining two folds amounts to left-multiplication of
  // transition matrices: TM(12)=TM(2)*TM(1).

  TrTr tt(tt0);
  const int n = tt0.edges();
  jlt::mathmatrix<int> TM(n,n);

  for (int i = 0; i < n; ++i)
    {
      tt = tt0;
      // Set all weights but one to zero.
      typename TrTr::dblVec wv(n);
      wv[i] = 1;
      tt.weights(wv.begin());

      // Compute new weights.
      tt.fold(f);
      wv = tt.weights();

      // Copy to ith row of matrix.
      for (int j = 0; j < n; ++j) TM(j,i) = (int)wv[j];
    }

  // Check that the transition matrix is permutation+1 or identity.
  if (TM != jlt::identity_matrix<int>(n))
    {
      bool notfoundcol2 = false, notfoundrow2 = false, bad = false;
      for (int i = 0; i < n; ++i)
	{
	  int colsum = 0, rowsum = 0;
	  for (int j = 0; j < n; ++j)
	    {
	      if (!(TM(i,j) == 0 || TM(i,j) == 1))
		{
		  std::cerr << "Matrix should contains only ones or zeros ";
		  std::cerr << " in ttauto::fold_transition_matrix.\n";
		  std::exit(1);
		}
	      colsum += TM(i,j);
	      rowsum += TM(j,i);
	    }
	  if (colsum == 2)
	    notfoundcol2 = false;
	  else if (colsum != 1)
	    bad = true;

	  if (rowsum == 2)
	    notfoundrow2 = false;
	  else if (rowsum != 1)
	    bad = true;
	}
      if (bad || notfoundrow2 || notfoundcol2)
	{
	  std::cerr << "Matrix should be permutation+1 or identity";
	  std::cerr << " in ttauto::fold_transition_matrix:\n";
	  TM.printMatrixForm(std::cerr);
	  std::exit(1);
	}
    }

  return TM;
}


// Small function to mod and make the result nonnegative.
template<class Int, class Int2>
inline Int2 mod(const Int m, const Int2 n)
{
  // Be careful with unsigned types for Int!
  Int mm = m % n;
  // mm has a value between -(n-1) and (n-1).  Make it nonnegative.
  return (Int2)(mm >= 0 ? mm : mm + n);
}


// Find the largest real root of a polynomial via Newton iteration.
//
// The polynomial comes from a Perron-Frobenius matrix, so the largest
// real root is guaranteed to be real and unique.
inline double findroot(const jlt::polynomial<int>& p,
		       const double x0, const double tol)
{
  double px(p(x0)), x(x0);

  int i = 0;
  const int itmax = 100;

  while (std::abs(px) > tol && i++ < itmax)
    {
      x = x - px / p.derivative_at(x);
      px = p(x);
    }

  if (i == itmax)
    throw
      jlt::failed_to_converge<double>
      ("Failed to converge to specified accuracy.\n",std::abs(px));

  return x;
}

} // namespace ttauto

#endif // TRAINTRACKS_UTIL_HPP
