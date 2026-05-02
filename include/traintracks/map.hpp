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

#ifndef TRAINTRACKS_MAP_HPP
#define TRAINTRACKS_MAP_HPP

#include <cstdlib>
#include <iostream>
#include <jlt/freeauto.hpp>
#include <jlt/freeword.hpp>
#include <jlt/mathmatrix.hpp>
#include <jlt/vector.hpp>
#include "traintracks/mathmatrix_permplus1.hpp"
#include "traintracks/map_labels.hpp"

namespace traintracks {

template<class TrTr>
mathmatrix_permplus1 fold_transition_matrix(const TrTr& tt0, const int f)
{
  // Conventions:
  // - We represent one fold by TM(f).
  // - Applying f1, then f2, composes by left-multiplication:
  //     TM_total = TM(f2) * TM(f1).

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

  // Validate fold matrix shape once through the canonical sparse decoder.
  // This enforces the same permutation/permutation+1 constraints used
  // everywhere else in the map pipeline.
  return mathmatrix_permplus1(TM);
}


// Validate that the map and transition matrix agree on main edges.
template<class TrTr>
inline void check_fold_map_main_transition(const TrTr& tt0, const int f,
					   const jlt::freeauto<int>& AM)
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
  jlt::mathmatrix<int> TM = fold_transition_matrix(tt0,f).full();

  if (TMfromAM != TM)
    {
      std::cerr << "Map/transition mismatch in traintracks::check_fold_map_main_transition.\n";
      std::exit(1);
    }
}


// Build one-fold train-track map on main and infinitesimal generators.
template<class TrTr>
jlt::freeauto<int> fold_traintrack_map(const TrTr& tt0, const int f)
{
  // Conventions:
  // - We represent one fold by AM(f).
  // - Applying f1 then f2 corresponds to right-composition in freeauto:
  //     AM(f2 followed by f1) = AM(f1) * AM(f2).
  // This matches freeauto<T>::operator*= semantics.

  const int ninf = tt0.total_prongs();
  const int n = tt0.edges();
  const ttmap_labeler labels(n,ninf);
  jlt::freeauto<int> AM(labels.num_generators());

  // Need to find two main edges and one infinitesimal edge.
  // main edge a: folding from
  // main edge b: folding onto
  // infinitesimal edge: in between

  // The cusp should give us the in-between edge.

  mathmatrix_permplus1 dec = fold_transition_matrix(tt0,f);

  // Copy permutation over to train track map (main generators only).
  for (int i = 0; i < n; ++i)
    AM[labels.main_gen(i)] = jlt::freeword<int>({labels.main_gen(dec.column_perm()[i])});

  if (dec.is_perm()) return AM;

  // Now deal with the folded edges.
  int e1 = labels.main_gen(dec.plus1_col());
  int e21 = labels.main_gen(dec.column_perm()[dec.plus1_col()]);
  int e22 = labels.main_gen(dec.plus1_row());

  // One main generator flips orientation under a non-permutation fold.
  AM[labels.main_gen(dec.row_perm()[dec.plus1_row()])] = {-e22};

  // Infinitesimal edge chosen from fold cusp geometry.
  int infinitesimal = tt0.fold_infinitesimal_generator(f,n);

  // Fold direction determines ordering around the inserted infinitesimal edge.
  if (f % 2 == 0)
    AM[e1] = {e21,infinitesimal,e22}; // fold clockwise
  else
    AM[e1] = {e22,infinitesimal,e21}; // fold counterclockwise

  // Main-edge transition consistency check (debug mode only).
  if (TrTr::debug) check_fold_map_main_transition(tt0,f,AM);

  return AM;
}


template<class TrTr>
jlt::mathmatrix<int> transition_matrix_from_map(const TrTr& tt,
					 const jlt::freeauto<int>& AM)
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

#endif // TRAINTRACKS_MAP_HPP
