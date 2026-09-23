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

#ifndef TRAINTRACKS_EMBEDDING_HPP
#define TRAINTRACKS_EMBEDDING_HPP

#include <iosfwd>
#include <vector>
#include "traintracks/coding.hpp"
#include "traintracks/fold_map.hpp"

namespace traintracks {

// The proper embedding of a train track in the punctured disc, as far as it
// is combinatorial: the walk around the boundary region, and the resulting
// order of the punctures.
//
// The multigons and main edges of a train track form a tree
// (traintrack::edges() == multigons()-1), so the cyclic order of prongs
// within a multigon and of slots within a prong -- which is what the coding
// records -- is a complete rotation system, and the track has one planar
// embedding up to isotopy and reflection.  Every main edge is a bridge, so
// the boundary region runs along both of its sides; every side (the loop of
// a punctured monogon, or an infinitesimal edge of an unpunctured polygon)
// has a complementary region on one side and the boundary region on the
// other, so the boundary region runs along it once.  Hence
//
//   walk.size() == 2*nedges() + nprongs()
//
// and the boundary region is the unique face of that length.  Each puncture
// is met exactly once, when the walk traverses its loop, and that is what
// puts the punctures in a cyclic order in the disc.
//
// Making the order linear -- punctures 1..n along the real line, as
// "properly embedded" requires (Lanneau and Thiffeault, Prop. 4) -- needs
// one more datum, namely where the cycle is cut.  The coding does not fix
// it: two tracks with the same coding can sit in the disc so that different
// punctures are leftmost.  The cut is a point of the boundary region, not a
// puncture, and it is named here by the dart of the walk that follows it.
// The canonical choice is the loop of the root monogon of the coding, which
// makes that monogon puncture 1; following a folding path instead carries
// the cut across each fold, since the punctures stay where they are while
// the track is folded and only the track changes.
struct tt_embedding
{
  // Signed letters (ttnumbering encoding) traversed by the boundary region,
  // in order, starting at the loop of the puncture at position 1.
  std::vector<int> walk;

  // Prong numbers of the punctured multigons, in walk order: entry i is the
  // puncture at position i+1.
  std::vector<int> puncture_order;

  // Prong number -> position 1..n, and 0 for a prong that is not the one
  // recording a puncture.
  std::vector<int> position_of;

  // The dart the walk starts at, that is the one just after the cut.
  int cut_dart = 0;

  // Corners of the boundary region between two consecutive main letters,
  // that is the exterior cusps, as indices into walk: the cusp lies between
  // walk[i] and walk[i+1].
  std::vector<int> cusp_after;

  int npunctures() const { return (int)puncture_order.size(); }

  std::ostream& print(std::ostream& strm) const;
};

// All directions at prong q in cyclic order, including the sides of
// unpunctured multigons: the incoming side reversed, the main letters in
// slot order, then the outgoing side.  ttnumbering::directions_at_prong
// omits the sides of unpunctured multigons, since a gate is about tangent
// directions and an infinitesimal edge is not one; walking the faces of the
// ribbon graph needs them.
std::vector<int> all_directions_at_prong(const ttnumbering& num, const int q);

// The boundary walk of a normalised track, starting at cut_dart.  Pass 0
// for the canonical choice, the loop of the root monogon (prong 0), which
// puts that puncture at position 1.  cut_dart must lie on the boundary
// region: a main letter of either sign, or a positive side letter.  Exits
// if the track is not the tree of multigons this assumes, or if a punctured
// multigon is not a monogon, which braid extraction cannot yet handle.
tt_embedding outer_embedding(const ttnumbering& num, const int cut_dart = 0);

// The dart that cut_dart of the track before a fold becomes after it.
int transported_cut_dart(const ttnumbering& before, const fold_map_data& fm,
                         const int cut_dart);

} // namespace traintracks

#endif // TRAINTRACKS_EMBEDDING_HPP
