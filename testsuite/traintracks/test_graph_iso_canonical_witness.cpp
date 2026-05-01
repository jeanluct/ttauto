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
#include "traintracks/graph_iso.hpp"
#include "traintracks/traintrack.hpp"

#define REQUIRE(cond) \
  do { \
    if (!(cond)) { \
      std::cerr << "Requirement failed: " << #cond << " at " \
                << __FILE__ << ":" << __LINE__ << "\n"; \
      return 1; \
    } \
  } while (0)

static traintracks::traintrack::intVec parse_compact_coding(const char* text)
{
  traintracks::traintrack::intVec out;
  std::istringstream in(text);
  std::string blk;
  while (in >> blk)
    {
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
  using traintracks::graph_iso::canonical_witness;
  using traintracks::graph_iso::canonical_witness_oriented;

  const traintrack::intVec code29 = parse_compact_coding(
    "11111 11113 11123 13111 23111 11111 33111 11111 "
    "11133 13111 23111 11111 33111 11111");
  const traintrack::intVec code71 = parse_compact_coding(
    "11111 11133 11113 13111 23111 11111 33111 11111 "
    "11123 13111 23111 11111 33111 11111");

  const traintrack t29(code29);
  const traintrack t71(code71);

  // Structural isotopy remains the primary invariant regardless of
  // operator== mode.
  REQUIRE(traintracks::graph_iso::is_isotopic_oriented(t29,t71));

#if !defined(TRAINTRACKS_USE_GRAPH_ISO_EQUALITY) || TRAINTRACKS_USE_GRAPH_ISO_EQUALITY
  const bool expect_operator_eq = true;
#else
  const bool expect_operator_eq = false;
#endif
  REQUIRE((t29 == t71) == expect_operator_eq);

  // Canonical witness should collapse equivalent representatives to the same
  // multigon ordering/ranking and prong shifts. This is the key deterministic
  // ingredient used by canonical fold-index APIs.
  const canonical_witness w29 = canonical_witness_oriented(t29);
  const canonical_witness w71 = canonical_witness_oriented(t71);
  REQUIRE(w29.valid);
  REQUIRE(w71.valid);
  REQUIRE((int)w29.multigon_order.size() == t29.multigons());
  REQUIRE((int)w71.multigon_order.size() == t71.multigons());
  REQUIRE(w29.multigon_order == w71.multigon_order);
  REQUIRE(w29.multigon_rank == w71.multigon_rank);
  REQUIRE(w29.prong_shift == w71.prong_shift);

  // Canonical fold/cusp APIs should agree across structurally equal tracks,
  // independent of the equality mode selected for operator==.
  REQUIRE(t29.foldings() == t71.foldings());
  for (int f = 0; f < t29.foldings(); ++f)
    {
      REQUIRE(t29.fold_infinitesimal_index_canonical(f)
              == t71.fold_infinitesimal_index_canonical(f));
      REQUIRE(t29.fold_infinitesimal_generator_canonical(f,100)
              == t71.fold_infinitesimal_generator_canonical(f,100));
    }

  return 0;
}
