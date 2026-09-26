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

#ifndef TRAINTRACKS_COLLAPSED_LAYOUT_HPP
#define TRAINTRACKS_COLLAPSED_LAYOUT_HPP

#include <utility>
#include <vector>
#include "traintracks/coding.hpp"
#include "traintracks/embedding.hpp"

namespace traintracks {

// A drawing of a train track in its collapsed representation: every
// multigon shrunk to a point and every main edge drawn as a chain of
// cubic Bezier arcs, joined smoothly.  This is what examples/ttplot turns into TikZ; it lives here so
// that its planarity can be tested.
//
// The drawing is a proper embedding: the punctures lie on the real axis
// at positions 1..n, in the order outer_embedding's boundary walk meets
// them, and the track lies in the upper half plane.  A multigon's prongs
// are equally spaced by angle about it, and every edge at a prong leaves
// along that prong's direction, so edges sharing a prong are tangent
// there.  Arcs should never cross; where they meet they are tangent.
// That holds for every automaton vertex for n = 3..7, which the testsuite
// checks, but not yet for all of n = 8.

struct vec2 { double x; double y; };

// A cubic Bezier arc from p0 to p3 with control points p1 and p2.
struct cubic { vec2 p0, p1, p2, p3; };

vec2 cubic_point(const cubic& c, const double t);

struct collapsed_layout
{
  std::vector<vec2> mg;                       // by multigon index
  std::vector<std::pair<int,int> > edge_mg;   // by edge: its two multigons
  std::vector<std::pair<int,int> > edge_pr;   // by edge: its two prongs
  std::vector<vec2> prong_dir;                // by prong: outgoing tangent
  std::vector<int> puncture_pos;              // by multigon, 1..n, 0 if none
  std::vector<std::vector<cubic> > arc;       // by edge: pieces, tail to head
};

// Requires every punctured multigon to be a monogon, as outer_embedding
// does.
collapsed_layout make_collapsed_layout(const ttnumbering& num,
                                       const tt_embedding& emb);

} // namespace traintracks

#endif // TRAINTRACKS_COLLAPSED_LAYOUT_HPP
