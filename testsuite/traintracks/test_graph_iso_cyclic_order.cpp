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

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include "traintracks/traintrack.hpp"

#define REQUIRE(cond) \
  do { \
    if (!(cond)) { \
      std::cerr << "Requirement failed: " << #cond << " at " \
                << __FILE__ << ":" << __LINE__ << "\n"; \
      return 1; \
    } \
  } while (0)

int main()
{
  using traintracks::traintrack;

  auto parse_compact_coding = [](const char* text)
  {
    traintrack::intVec out;
    std::istringstream in(text);
    std::string blk;
    while (in >> blk)
      {
        if (blk.size() != 5)
          {
            std::cerr << "Bad coding block length: '" << blk << "'\n";
            std::exit(2);
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
  };

  // Three variants with the same underlying branch-incidence multigraph but
  // different local prong presentations in compact coding form.
  const traintrack tA(parse_compact_coding(
    "11111 11113 11123 13111 23111 11111 33111 11111 "
    "11133 13111 23111 11111 33111 11111"));

  const traintrack tB(parse_compact_coding(
    "11111 11133 11113 13111 23111 11111 33111 11111 "
    "11123 13111 23111 11111 33111 11111"));

  // A control track with different topology.
  const traintrack tC(5,3);

  REQUIRE(tA.coding() != tB.coding());
  REQUIRE(tA == tB);

  // Trivial self-equality and a basic negative control.
  REQUIRE(tA == tA);
  REQUIRE(!(tA == tC));
  REQUIRE(!(tB == tC));

  return 0;
}
