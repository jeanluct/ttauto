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
#include <vector>
#include <jlt/freeauto.hpp>
#include <jlt/mathmatrix.hpp>
#include "traintracks/graph_iso.hpp"
#include "traintracks/map.hpp"
#include "traintracks/traintrack.hpp"

namespace traintracks { class multigon; }

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
}

static int native_fold_index_from_canonical(const traintracks::traintrack& tt,
                                            const int fcanon)
{
  traintracks::multigon* mmc = 0;
  int pc = -1, ec = -1;
  tt.fold_cusp_location_canonical(fcanon,mmc,pc,ec);
  if (mmc == 0 || pc < 0 || ec < 0)
    {
      std::cerr << "Bad canonical cusp location.\n";
      std::exit(2);
    }

  auto locate_multigon_index = [&](const traintracks::multigon* needle)
  {
    for (int m = 0; m < tt.multigons(); ++m)
      {
        if (needle == &tt.Multigon(m)) return m;
      }
    return -1;
  };

  const int target_m = locate_multigon_index(mmc);
  if (target_m < 0)
    {
      std::cerr << "Canonical cusp multigon is not in track.\n";
      std::exit(2);
    }
  const int want_dir = (fcanon % 2);

  for (int f = want_dir; f < tt.foldings(); f += 2)
    {
      traintracks::multigon* mm = 0;
      int p = -1, e = -1;
      tt.fold_cusp_location(f,mm,p,e);
      if (mm == 0) continue;
      if (locate_multigon_index(mm) == target_m && p == pc && e == ec)
        return f;
    }

  std::cerr << "Could not map canonical fold index to native fold index.\n";
  std::exit(2);
}

// Build a relabeling-invariant fingerprint for permutation/permutation+1
// transition matrices under independent row/column relabelings.
//
// Performance rationale:
// - the previous implementation canonicalized dense matrices by exhaustive
//   bi-permutation search (factorial cost); deterministic but slow (~30s);
// - fold matrices in this code path are already decoded as mathmatrix_permplus1,
//   so we can fingerprint directly from sparse structure in O(1).
//
// Invariance rationale:
// - under independent row/column relabelings, all permutation matrices of a
//   fixed size are equivalent;
// - likewise, all proper permutation+1 matrices of a fixed size are equivalent;
// - therefore `(dim, is_perm)` is the complete class key in this quotient.
static std::vector<int> biperm_class_key(const traintracks::mathmatrix_permplus1& PM)
{
  std::vector<int> key;
  key.push_back(PM.dim());
  if (PM.is_perm())
    {
      key.push_back(0);
      return key;
    }

  // Non-permutation folds must expose one valid +1 location.
  const int p1r = PM.plus1_row();
  const int p1c = PM.plus1_col();
  if (p1r < 0 || p1c < 0)
    {
      std::cerr << "Missing +1 location in non-permutation fold matrix.\n";
      std::exit(2);
    }

  key.push_back(1);
  return key;
}

int main()
{
  using traintracks::traintrack;
  using traintracks::mathmatrix_permplus1;
  using traintracks::fold_transition_matrix;
  using traintracks::fold_traintrack_map;
  using traintracks::transition_matrix_from_map;

  // Same pair as issue #12: they differ in full cusp-slot-aware structure but
  // are useful as a deterministic stress case for canonical fold APIs.
  const traintrack::intVec code29 = parse_compact_coding(
    "11111 11113 11123 13111 23111 11111 33111 11111 "
    "11133 13111 23111 11111 33111 11111");
  const traintrack::intVec code71 = parse_compact_coding(
    "11111 11133 11113 13111 23111 11111 33111 11111 "
    "11123 13111 23111 11111 33111 11111");

  const traintrack t29(code29);
  const traintrack t71(code71);

  REQUIRE(!traintracks::graph_iso::is_isotopic_oriented(t29,t71));
  REQUIRE(t29.foldings() == t71.foldings());
  REQUIRE(t29.edges() == t71.edges());
  REQUIRE(t29.total_prongs() == t71.total_prongs());

  const int nmain = t29.edges();

  // For each canonical fold index f:
  // 1) canonical infinitesimal generator is valid on each representative,
  // 2) train-track map projected to main edges reproduces fold transition
  //    matrix for each representative.
  for (int f = 0; f < t29.foldings(); ++f)
    {
      const int g29 = t29.fold_infinitesimal_generator_canonical(f,nmain);
      const int g71 = t71.fold_infinitesimal_generator_canonical(f,nmain);
      REQUIRE(g29 >= nmain + 1);
      REQUIRE(g71 >= nmain + 1);

      const int f29 = native_fold_index_from_canonical(t29,f);
      const int f71 = native_fold_index_from_canonical(t71,f);

      traintracks::multigon* mm29 = 0;
      traintracks::multigon* mm71 = 0;
      int p29 = -1, e29 = -1;
      int p71 = -1, e71 = -1;
      t29.fold_cusp_location_canonical(f,mm29,p29,e29);
      t71.fold_cusp_location_canonical(f,mm71,p71,e71);
      REQUIRE(mm29 != 0);
      REQUIRE(mm71 != 0);
      REQUIRE(p29 >= 0);
      REQUIRE(p71 >= 0);
      REQUIRE(e29 >= 0);
      REQUIRE(e71 >= 0);

      const mathmatrix_permplus1 TM29 = fold_transition_matrix(t29,f29);
      const mathmatrix_permplus1 TM71 = fold_transition_matrix(t71,f71);

      const jlt::freeauto<int> AM29 = fold_traintrack_map(t29,f29);
      const jlt::freeauto<int> AM71 = fold_traintrack_map(t71,f71);
      const jlt::mathmatrix<int> TM29fromAM = transition_matrix_from_map(t29,AM29);
      const jlt::mathmatrix<int> TM71fromAM = transition_matrix_from_map(t71,AM71);
      REQUIRE(TM29fromAM == TM29.full());
      REQUIRE(TM71fromAM == TM71.full());

      REQUIRE(biperm_class_key(TM29)
              == biperm_class_key(traintracks::mathmatrix_permplus1(TM29fromAM)));
      REQUIRE(biperm_class_key(TM71)
              == biperm_class_key(traintracks::mathmatrix_permplus1(TM71fromAM)));
    }

  return 0;
}
