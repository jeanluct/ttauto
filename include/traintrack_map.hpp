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

#ifndef TRAINTRACK_MAP_HPP
#define TRAINTRACK_MAP_HPP

#include <cstdlib>
#include <iostream>
#include <jlt/mathmatrix.hpp>
#include <jlt/vector.hpp>
#include "freeauto.hpp"
#include "ttmap_labels.hpp"

namespace traintracks {

template<class TrTr>
jlt::mathmatrix<int> fold_transition_matrix(const TrTr& tt0, const int f);

struct permplus1_decode
{
  // True when the matrix is a pure permutation (no extra +1 entry).
  bool is_perm;
  // Row/column indices of the extra +1 when is_perm is false.
  int row2;
  int col2;
  // Row permutation and its inverse, with optional +1 entry removed.
  jlt::mathvector<int> perm;
  jlt::mathvector<int> perm_inv;
};


inline permplus1_decode decode_permplus1(const jlt::mathmatrix<int>& TM)
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
              std::cerr << "traintracks::decode_permplus1.\n";
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
              std::cerr << "traintracks::decode_permplus1.\n";
              std::exit(1);
            }
          d.row2 = i;
        }
      else if (colsum != 1)
        {
          std::cerr << "Bad row sum in traintracks::decode_permplus1.\n";
          std::exit(1);
        }

      if (rowsum == 2)
        {
          if (d.col2 != -1)
            {
              std::cerr << "More than one +1 col in ";
              std::cerr << "traintracks::decode_permplus1.\n";
              std::exit(1);
            }
          d.col2 = i;
        }
      else if (rowsum != 1)
        {
          std::cerr << "Bad col sum in traintracks::decode_permplus1.\n";
          std::exit(1);
        }
    }

  d.is_perm = (d.row2 == -1 && d.col2 == -1);
  if ((d.row2 == -1) != (d.col2 == -1))
    {
      std::cerr << "Inconsistent +1 row/col in ";
      std::cerr << "traintracks::decode_permplus1.\n";
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
        std::cerr << "traintracks::decode_permplus1.\n";
        std::exit(1);
      }

  return d;
}


// Decode fold structure once and convert to map-indexing conventions.
template<class TrTr>
inline permplus1_decode decode_fold_map_structure(const TrTr& tt0, const int f)
{
  // fold_traintrack_map uses generator-action composition conventions where
  // permutation data is read in the dual (column-action) indexing. Decode the
  // non-transposed transition matrix, then convert indices/permutations once.
  permplus1_decode d0 = decode_permplus1(fold_transition_matrix(tt0,f));
  const int n = d0.perm.size();
  permplus1_decode d{d0.is_perm,d0.col2,d0.row2,
                     jlt::mathvector<int>(n,-1),jlt::mathvector<int>(n,-1)};
  for (int i = 0; i < n; ++i)
    {
      d.perm[i] = d0.perm_inv[i];
      d.perm_inv[i] = d0.perm[i];
    }
  return d;
}


// Validate that the map and transition matrix agree on main edges.
template<class TrTr>
inline void check_fold_map_main_transition(const TrTr& tt0, const int f,
					   const freeauto<int>& AM)
{
  const int n = tt0.edges();
  const ttmap_labeler labels(n,tt0.total_prongs());
  jlt::mathmatrix<int> TMfromAM(n,n,0);

  for (int src = 0; src < n; ++src)
    {
      const int g = src + 1;
      for (auto img : AM.get_action(g))
	{
	  if (!labels.is_valid_generator(img))
	    {
	      std::cerr << "Bad generator in traintracks::check_fold_map_main_transition.\n";
	      std::exit(1);
	    }
	  if (!labels.is_main_generator(img)) continue;
	  int col = labels.main_generator_index(img);
	  ++TMfromAM(col,src);
	}
    }
  jlt::mathmatrix<int> TM = fold_transition_matrix(tt0,f);

  if (TMfromAM != TM)
    {
      std::cerr << "Map/transition mismatch in traintracks::check_fold_map_main_transition.\n";
      std::exit(1);
    }
}


// Build one-fold train-track map on main and infinitesimal generators.
template<class TrTr>
freeauto<int> fold_traintrack_map(const TrTr& tt0, const int f)
{
  // Conventions:
  // - We represent one fold by AM(f).
  // - Applying f1 then f2 corresponds to right-composition in freeauto:
  //     AM(f2 followed by f1) = AM(f1) * AM(f2).
  // This matches freeauto<T>::operator*= semantics.

  const int ninf = tt0.total_prongs();
  const int n = tt0.edges();
  const ttmap_labeler labels(n,ninf);
  freeauto<int> AM(labels.num_generators());

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

  // Main-edge transition consistency check.
  check_fold_map_main_transition(tt0,f,AM);

  return AM;
}


template<class TrTr>
jlt::mathmatrix<int> transition_matrix_from_map(const TrTr& tt,
					 const freeauto<int>& AM)
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


} // namespace traintracks

#endif // TRAINTRACK_MAP_HPP
