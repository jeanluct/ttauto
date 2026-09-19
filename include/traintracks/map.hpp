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
#include "traintracks/fold_map.hpp"

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


// Build the one-fold train-track map on main and side (infinitesimal or
// peripheral) generators, in the canonical numberings of tt0 and of the
// folded track.  See fold_map.hpp for the conventions; the fold is
// performed on a copy of tt0.  Composition convention, as for
// freeauto<T>::operator*=: applying f1 then f2 is AM(f1) * AM(f2).
template<class TrTr>
jlt::freeauto<int> fold_traintrack_map(const TrTr& tt0, const int f)
{
  TrTr tt(tt0);
  fold_map_data fm;
  if (!tt.fold_with_map(f,fm))
    return identity_traintrack_map(tt0.edges(),tt0.total_prongs());

  // Main-edge transition consistency check (debug mode only).
  if (TrTr::debug) check_fold_map_main_transition(tt0,f,fm.to_freeauto());

  return fm.to_freeauto();
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
