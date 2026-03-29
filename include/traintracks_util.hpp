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
#include "freeword.hpp"
#include "ttmap_labels.hpp"

namespace traintracks {

struct permplus1_decode
{
  bool is_perm;
  int row2;
  int col2;
  jlt::mathvector<int> perm;
  jlt::mathvector<int> perm_inv;
};


inline permplus1_decode decode_transposed_permplus1(const jlt::mathmatrix<int>& TM)
{
  const int n = TM.dim();
  permplus1_decode d{true,-1,-1,jlt::mathvector<int>(n,-1),jlt::mathvector<int>(n,-1)};

  for (int i = 0; i < n; ++i)
    {
      int colsum = 0, rowsum = 0;
      for (int j = 0; j < n; ++j)
        {
          if (!(TM(i,j) == 0 || TM(i,j) == 1))
            {
              std::cerr << "Matrix should contain only ones or zeros in ";
              std::cerr << "traintracks::decode_transposed_permplus1.\n";
              std::exit(1);
            }
          colsum += TM(i,j);
          rowsum += TM(j,i);
        }

      if (colsum == 2)
        {
          if (d.row2 != -1)
            {
              std::cerr << "More than one +1 row in ";
              std::cerr << "traintracks::decode_transposed_permplus1.\n";
              std::exit(1);
            }
          d.row2 = i;
        }
      else if (colsum != 1)
        {
          std::cerr << "Bad row sum in traintracks::decode_transposed_permplus1.\n";
          std::exit(1);
        }

      if (rowsum == 2)
        {
          if (d.col2 != -1)
            {
              std::cerr << "More than one +1 col in ";
              std::cerr << "traintracks::decode_transposed_permplus1.\n";
              std::exit(1);
            }
          d.col2 = i;
        }
      else if (rowsum != 1)
        {
          std::cerr << "Bad col sum in traintracks::decode_transposed_permplus1.\n";
          std::exit(1);
        }
    }

  d.is_perm = (d.row2 == -1 && d.col2 == -1);
  if ((d.row2 == -1) != (d.col2 == -1))
    {
      std::cerr << "Inconsistent +1 row/col in ";
      std::cerr << "traintracks::decode_transposed_permplus1.\n";
      std::exit(1);
    }

  jlt::mathmatrix<int> P(TM);
  if (!d.is_perm) P(d.row2,d.col2) = 0;

  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      if (P(i,j) == 1)
        {
          d.perm[i] = j;
          d.perm_inv[j] = i;
        }

  for (int i = 0; i < n; ++i)
    if (d.perm[i] < 0 || d.perm_inv[i] < 0)
      {
        std::cerr << "Not a permutation(+1) matrix in ";
        std::cerr << "traintracks::decode_transposed_permplus1.\n";
        std::exit(1);
      }

  return d;
}


template<class TrTr>
inline permplus1_decode decode_fold_map_structure(const TrTr& tt0, const int f)
{
  // Main-edge map (AM) is represented in transposed transition convention,
  // so decode the transposed fold transition matrix.
  return decode_transposed_permplus1(fold_transition_matrix(tt0,f).transpose());
}

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
  // Conventions:
  // - We represent one fold by TM(f).
  // - Applying f1 then f2 corresponds to left-multiplication:
  //     TM(f2 followed by f1) = TM(f2) * TM(f1).

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
      bool notfoundcol2 = true, notfoundrow2 = true, bad = false;
      for (int i = 0; i < n; ++i)
	{
	  int colsum = 0, rowsum = 0;
	  for (int j = 0; j < n; ++j)
	    {
	      if (!(TM(i,j) == 0 || TM(i,j) == 1))
		{
		  std::cerr << "Matrix should contains only ones or zeros ";
		  std::cerr << " in traintracks::fold_transition_matrix.\n";
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
	  std::cerr << " in traintracks::fold_transition_matrix:\n";
	  TM.printMatrixForm(std::cerr);
	  std::exit(1);
	}
    }

  return TM;
}


template<class TrTr>
free_auto<int> fold_traintrack_map(const TrTr& tt0, const int f)
{
  // Conventions:
  // - We represent one fold by AM(f).
  // - Applying f1 then f2 corresponds to right-composition in free_auto:
  //     AM(f2 followed by f1) = AM(f1) * AM(f2).
  // This matches free_auto<T>::operator*= semantics.

  const int ninf = tt0.total_prongs();
  const int n = tt0.edges();
  const ttmap_labeler labels(n,ninf);
  free_auto<int> AM(labels.num_generators());

  // Need to find two main edges and one infinitesimal edge.
  // main edge a: folding from
  // main edge b: folding onto
  // infinitesimal edge: in between

  // The cusp should give us the in-between edge.

  permplus1_decode dec = decode_fold_map_structure(tt0,f);

  // Copy permutation over to train track map (main generators only).
  for (int i = 0; i < n; ++i)
    AM[labels.main_gen(i)] = {labels.main_gen(dec.perm[i])};

  if (dec.is_perm) return AM;

  // Now deal with the folded edges.
  int e1 = labels.main_gen(dec.row2);
  int e21 = labels.main_gen(dec.perm[dec.row2]);
  int e22 = labels.main_gen(dec.col2);

  // One main generator flips orientation under a non-permutation fold.
  AM[labels.main_gen(dec.perm_inv[dec.col2])] = {-e22};

  // Infinitesimal edge chosen from fold cusp geometry.
  int infinitesimal = tt0.fold_infinitesimal_generator(f,n);

  // Fold direction determines ordering around the inserted infinitesimal edge.
  if (f % 2 == 0)
    AM[e1] = {e21,infinitesimal,e22}; // fold counterclockwise
  else
    AM[e1] = {e22,infinitesimal,e21}; // fold clockwise

  /* Check automorphism (see transition matrix) */

  return AM;
}


template<class TrTr>
jlt::mathmatrix<int> transition_matrix_from_map(const TrTr& tt,
						const free_auto<int>& AM)
{
  // Returns the main-edge transition matrix in non-transposed form.
  // Infinitesimal generators are ignored for this projection.
  const int n = tt.edges();
  const ttmap_labeler labels(n,tt.total_prongs());
  jlt::mathmatrix<int> TM(n,n,0);

  if ((int)AM.numgens() < n)
    {
      std::cerr << "Automorphism has too few generators in ";
      std::cerr << "traintracks::transition_matrix_from_map.\n";
      std::exit(1);
    }

  for (int src = 0; src < n; ++src)
    {
      const int g = src + 1;
      for (auto img : AM.get_action(g))
	{
	  if (!labels.is_valid_generator(img))
	    {
	      std::cerr << "Bad generator in ";
	      std::cerr << "traintracks::transition_matrix_from_map.\n";
	      std::exit(1);
	    }
	  if (!labels.is_main_generator(img)) continue;
	  int col = labels.main_generator_index(img);
	  ++TM(col,src);
	}
    }

  return TM;
}


template<class TrTr>
jlt::mathmatrix<int> transition_matrix_from_map_transposed(const TrTr& tt,
							   const free_auto<int>& AM)
{
  jlt::mathmatrix<int> TM = transition_matrix_from_map(tt,AM);
  TM.transpose();
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

} // namespace traintracks

#endif // TRAINTRACKS_UTIL_HPP
