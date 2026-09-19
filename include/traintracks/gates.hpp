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

#ifndef TRAINTRACKS_GATES_HPP
#define TRAINTRACKS_GATES_HPP

#include <iosfwd>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <jlt/freeauto.hpp>
#include "traintracks/coding.hpp"

namespace traintracks {

// Bestvina-Handel gate test for a train-track map, without composing
// words.  See devel/iss002/issue2_gates.tex, Sections 2, 3 and 6.
//
// Dictionary.  A Bestvina-Handel vertex is an unpunctured multigon taken
// as a whole, or one prong of a punctured multigon.  The directions at a
// vertex are the signed main letters whose tail is there and, at a
// punctured prong, the two peripheral directions (the incoming side
// reversed and the outgoing side).  Sides of unpunctured multigons are
// infinitesimal edges, not directions.  The derivative D sends a letter
// to the first letter of its image; two directions at the same vertex are
// in the same gate when some power of D identifies them.  An image word
// x S y takes the turn (xbar, y) at the multigon of S if S is
// infinitesimal, and the two turns (xbar, S), (Sbar, y) if S is
// peripheral.  Gates joined by a realised turn are joined by an
// infinitesimal edge of the Bestvina-Handel train track; the map is
// pseudo-Anosov only if at every vertex these edges connect all the gates
// (Bestvina-Handel Prop. 3.3.2 and Section 3.4).
//
// Composition (Section 6 of the note): D(g o h) = D(g) o D(h) and
// T(g o h) = T(g) u D(g)(T(h)), so a path is analysed by pushing one
// fold_derivative per branch into a gate_accumulator.

// Letters are signed integers in 1..ngen (ngen = main edges + sides).
typedef std::pair<int,int> turn;   // unordered pair, stored with first < second

// The derivative and realised turns of one train-track map whose image
// words live in the numbering `target`.
struct fold_derivative
{
  int ngen = 0;
  std::vector<int> D;      // indexed by letter + ngen, size 2*ngen + 1
  std::set<turn> turns;    // realised turns, in the target's letters

  fold_derivative() {}
  // From any train-track map (one-step or composed); words are read in the
  // numbering of the target track.
  fold_derivative(const jlt::freeauto<int>& AM, const ttnumbering& target);

  int apply(const int letter) const { return D[letter + ngen]; }
};

// Result of the gate test at one vertex.
struct gate_vertex_report
{
  std::string name;                       // human-readable vertex description
  int multigon = -1;                      // multigon index
  int prong = -1;                         // prong number, or -1 for a whole unpunctured multigon
  std::vector<int> directions;            // in cyclic order
  std::vector<std::vector<int> > gates;   // partition of directions, in cyclic order
  std::vector<turn> joins;                // realised turns between distinct gates
  int components = 0;                     // connected components of the gate graph
  bool connected = false;
  bool shape_ok = false;                  // Bestvina-Handel Props. 3.3.3-3.3.4 pattern
};

struct gate_analysis
{
  bool connected = false;        // every vertex connected: passes the test
  bool shapes_ok = false;        // every vertex has an allowed infinitesimal-edge pattern
  std::vector<gate_vertex_report> vertices;

  std::ostream& print(std::ostream& strm) const;
};

// Accumulates derivative data along a folding path.  push_back() the
// fold_derivative of each branch in path order; analyse() at the end with
// the numbering of the initial (= final) track of a closed path.
class gate_accumulator
{
public:
  gate_accumulator() {}
  explicit gate_accumulator(const int ngen);

  int ngen() const { return ng; }
  int length() const { return len; }

  void push_back(const fold_derivative& d);

  const std::vector<int>& derivative() const { return Dacc; }
  const std::set<turn>& realised_turns() const { return Tacc; }

  // Gate test for the composed map, read as a self-map of the track with
  // numbering N (the path must be closed).
  gate_analysis analyse(const ttnumbering& N) const;

private:
  int ng = 0;
  int len = 0;
  std::vector<int> Dacc;   // identity initially
  std::set<turn> Tacc;
};

// Word-based version for tests and diagnostics: gate test of a composed
// map AM on the track with numbering N (AM must be a self-map of that
// track).  Must agree with the accumulator.
gate_analysis analyse_gates(const ttnumbering& N, const jlt::freeauto<int>& AM);

// Bestvina-Handel vertex of a direction: the prong number for a punctured
// multigon, or -(multigon index + 1) for an unpunctured multigon.
int bh_vertex_of(const ttnumbering& N, const int letter);

} // namespace traintracks

#endif // TRAINTRACKS_GATES_HPP
