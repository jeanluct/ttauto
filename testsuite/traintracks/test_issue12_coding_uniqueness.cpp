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
#include <sstream>
#include <string>
#include "traintracks/traintrack.hpp"
#include "ttauto/ttfoldgraph.hpp"


static traintracks::traintrack::intVec parse_compact_coding(const char* text)
{
  traintracks::traintrack::intVec out;
  std::istringstream in(text);
  std::string blk;
  while (in >> blk)
    {
      assert(blk.size() == 5);
      for (std::string::size_type i = 0; i < blk.size(); ++i)
        {
          assert(blk[i] >= '1' && blk[i] <= '9');
        }

      const int prong = (blk[0] - '0') - 1;
      const int nprongs = (blk[1] - '0');
      const int label = (blk[2] - '0') - 1;
      const int edge = (blk[3] - '0') - 1;
      const int nedges = (blk[4] - '0');

      out.push_back(prong);
      out.push_back(nprongs);
      out.push_back(label);
      out.push_back(edge);
      out.push_back(nedges);
    }
  return out;
}


int main()
{
  using traintracks::traintrack;
  using ttauto::ttfoldgraph;

  // Issue #12 reproducer: two codings observed as topologically identical
  // in automaton exploration/plotting, but currently not deduplicated.
  const traintrack::intVec code29 = parse_compact_coding(
    "11111 11113 11123 13111 23111 11111 33111 11111 "
    "11133 13111 23111 11111 33111 11111");

  const traintrack::intVec code71 = parse_compact_coding(
    "11111 11133 11113 13111 23111 11111 33111 11111 "
    "11123 13111 23111 11111 33111 11111");

  const traintrack tt29(code29);
  const traintrack tt71(code71);

  // Guard that we are still reproducing the reported discrepancy now.
  assert(tt29.coding() != tt71.coding());
  assert(!(tt29 == tt71));

  // This is intentionally a reproducer test at this stage.  When the
  // canonicalization fix lands, this test should be updated to assert
  // equality and coding uniqueness instead.

  // Build the reported automaton and pin the currently duplicated
  // representatives by 0-based vertex index (issue lists them as 29 and 71
  // in 1-based indexing).
  traintrack seed(6,3);
  ttfoldgraph<traintrack> ttg(seed);
  assert(ttg.vertices() == 90);

  const int v28 = 28;
  const int v70 = 70;
  assert(v28 < ttg.vertices());
  assert(v70 < ttg.vertices());

  [[maybe_unused]] const traintrack& g28 = ttg.traintrack(v28);
  [[maybe_unused]] const traintrack& g70 = ttg.traintrack(v70);

  // Current behavior: both vertices are present separately and keep distinct
  // codings despite representing the same topology.
  assert(g28.coding() == code29);
  assert(g70.coding() == code71);
  assert(g28.coding() != g70.coding());
  assert(!(g28 == g70));

  return 0;
}
