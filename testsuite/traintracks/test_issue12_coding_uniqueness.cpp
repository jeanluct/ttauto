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

#include <iostream>
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
      if (blk.size() != 5)
        {
          std::cerr << "Bad coding block length: '" << blk << "'\n";
          std::exit(2);
        }
      for (std::string::size_type i = 0; i < blk.size(); ++i)
        {
          if (!(blk[i] >= '1' && blk[i] <= '9'))
            {
              std::cerr << "Bad coding digit in block: '" << blk << "'\n";
              std::exit(2);
            }
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

static std::string to_compact_coding(const traintracks::traintrack::intVec& code)
{
  std::ostringstream out;
  for (std::size_t i = 0; i + 4 < code.size(); i += 5)
    {
      const int prong = code[i+0] + 1;
      const int nprongs = code[i+1];
      const int label = code[i+2] + 1;
      const int edge = code[i+3] + 1;
      const int nedges = code[i+4];
      out << prong << nprongs << label << edge << nedges;
      if (i + 5 < code.size()) out << " ";
    }
  return out.str();
}

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
  using ttauto::ttfoldgraph;

  const bool debug_issue12 = false;

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

  if (debug_issue12)
    {
      std::cout << "tt29==tt71: " << (tt29 == tt71) << "\n";
    }

  if (debug_issue12)
    {
      std::cout << "tt29 old coding:      " << to_compact_coding(tt29.coding()) << "\n";
      std::cout << "tt29 canonical coding:" << to_compact_coding(tt29.canonical_coding()) << "\n";
      std::cout << "tt71 old coding:      " << to_compact_coding(tt71.coding()) << "\n";
      std::cout << "tt71 canonical coding:" << to_compact_coding(tt71.canonical_coding()) << "\n";
    }

  // Orientation-preserving graph-isotopy equality should identify these.
  REQUIRE(tt29.coding() != tt71.coding());
  REQUIRE(tt29.canonical_coding() != tt71.canonical_coding());
  REQUIRE(tt29 == tt71);

  // Build the reported automaton and verify reflection pairing remains a
  // graph-level notion based on coding reversal.
  traintrack seed(6,3);
  ttfoldgraph<traintrack> ttg(seed);
  if (debug_issue12)
    {
      std::cout << "ttg vertices: " << ttg.vertices() << "\n";
    }
  REQUIRE(ttg.vertices() < 90);

  if (debug_issue12)
    {
      int ncanon_diff = 0;
      for (int i = 0; i < ttg.vertices(); ++i)
        {
          const traintrack& ti = ttg.traintrack(i);
          if (ti.coding() != ti.canonical_coding()) ++ncanon_diff;
        }
      std::cout << "vertices with coding!=canonical: " << ncanon_diff << "\n";
    }

  int npairs = 0;
  for (int i = 0; i < ttg.vertices(); ++i)
    {
      const traintrack& ti = ttg.traintrack(i);
      const int j = ttg.symmetric_double(i);
      if (j == -1) continue;

      REQUIRE(j >= 0 && j < ttg.vertices());
      REQUIRE(ttg.symmetric_double(j) == i);

      if (j <= i) continue;
      ++npairs;

      const traintrack& tj = ttg.traintrack(j);
      REQUIRE(ti.coding() == tj.coding(-1));
    }
  REQUIRE(npairs > 0);

  return 0;
}
