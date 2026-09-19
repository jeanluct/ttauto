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

#ifndef TRAINTRACKS_FOLD_MAP_HPP
#define TRAINTRACKS_FOLD_MAP_HPP

#include <iosfwd>
#include <vector>
#include <jlt/freeauto.hpp>
#include "traintracks/coding.hpp"

namespace traintracks {

// Everything one fold does to the edges and prongs of a train track,
// expressed in the canonical numberings (ttnumbering) of the track before
// and after the fold.  Filled by traintrack::fold_with_map().
//
// Conventions (see devel/iss002/issue2_gates.tex, Section 4):
//
// - The fold takes the cusp between two edges at one prong of the cusp
//   multigon and slides the end of edge `moved` along edge `onto` to the
//   target multigon, where `onto` is attached at prong target_from; the
//   moved end then runs along the side of the target multigon to the
//   adjacent prong target_to and re-attaches there.
// - Every edge other than `moved` maps to itself, relabelled:
//   edge_image[e] is the signed new letter of old edge e, the sign
//   recording whether the canonical orientations before and after agree.
// - The old edge `moved`, in its old canonical orientation, maps to the
//   three-letter path moved_word = {a, S, b} where S = side is the side of
//   the target multigon traversed between target_from and target_to (in
//   the new numbering, signed), and {a, b} are the new letters of `moved`
//   and `onto` in the order in which the path runs through them.
//   moved_first tells whether the moved edge's new copy comes first.
// - Side q of the old track becomes side prong_image[q] of the new one.
//
// Letters follow ttnumbering: main edge e is +-(e+1), side q is
// +-(nmain+1+q).  The struct is only valid after a successful fold.
struct fold_map_data
{
  int nmain = 0;           // number of main edges
  int nsides = 0;          // number of sides (= total prongs)

  int moved = -1;          // old edge number of the folded (moved) edge
  int onto = -1;           // old edge number of the edge it is folded onto
  int dir = 0;             // fold direction, +1 or -1, as in traintrack::fold
  int cusp_prong = -1;     // old prong number of the cusp
  int target_from = -1;    // new prong number where `onto` meets the target
  int target_to = -1;      // new prong number where `moved` re-attaches
  int side = 0;            // signed new side letter traversed by the path
  bool moved_first = false;
  std::vector<int> moved_word;    // three signed new letters

  std::vector<int> edge_image;    // old edge -> signed new letter
  std::vector<int> prong_image;   // old prong -> new prong number

  ttnumbering before;
  ttnumbering after;

  // The one-fold train-track map as a free-group automorphism on
  // nmain + nsides generators, images in the new numbering.
  jlt::freeauto<int> to_freeauto() const;

  std::ostream& print(std::ostream& strm) const;
};

// Identity map on nmain + nsides generators (what a fold that cannot be
// performed contributes).
jlt::freeauto<int> identity_traintrack_map(const int nmain, const int nsides);

} // namespace traintracks

#endif // TRAINTRACKS_FOLD_MAP_HPP
