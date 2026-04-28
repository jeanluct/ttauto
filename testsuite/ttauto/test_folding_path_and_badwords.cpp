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

#include <cassert>
#include <unordered_map>
#include "traintracks/traintrack.hpp"
#include "ttauto/badwords.hpp"
#include "ttauto/folding_path.hpp"
#include "ttauto/ttfoldgraph.hpp"


int main()
{
  using traintracks::traintrack;
  using ttauto::folding_path;

  typedef traintrack TrTr;

  const int n = 6;
  traintrack tt(n,3);
  ttauto::ttfoldgraph<TrTr> ttg(tt);

  // Graph should be non-empty for this deterministic fixture.
  assert(ttg.vertices() > 0);

  folding_path<TrTr> p1(ttg,1);
  p1.push_back(0);
  p1.push_back(0);
  p1.push_back(1);
  p1.push_back(1);

  folding_path<TrTr> p2(p1);
  p2.cycle_path(-1);

  // Closed-path cyclic equality must hold after a one-step cycle.
  assert(p1 == p2);

  typedef folding_path<TrTr>::hash path_hash;
  typedef std::unordered_map<folding_path<TrTr>,int,path_hash> pathlist;
  pathlist pl;
  pl[p1] = 1;
  pl[p2] = 2;

  // Hash/equality policy should deduplicate cyclically equivalent closed paths.
  assert(pl.size() == 1);

  const int maxplen = 4;

  // badwords() should return a matrix indexed by (vertex, path-length).
  auto pbad = ttauto::badwords(ttg,maxplen);
  assert(pbad.dim1() == ttg.vertices());
  assert(pbad.dim2() == maxplen);

  return 0;
}
