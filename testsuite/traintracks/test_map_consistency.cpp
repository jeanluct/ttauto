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

// One-fold train-track maps: agreement with the transition matrix, the
// side letter lives on the target multigon, and hand-checked words for the
// n=3 example (see devel/iss002/issue2_gates.tex, Section 4).

#include <algorithm>
#include <cstdlib>
#include <jlt/freeauto.hpp>
#include <jlt/stlio.hpp>
#include "check.hpp"
#include "oracles.hpp"
#include "traintracks/build.hpp"
#include "traintracks/fold_map.hpp"
#include "traintracks/map.hpp"
#include "traintracks/traintrack.hpp"
#include "ttauto/folding_path.hpp"
#include "ttauto/ttfoldgraph.hpp"


// Fold index (in traintrack::fold numbering) realising branch b at a vertex.
static int fold_index_of_branch(const traintracks::traintrack& tt, const int b)
{
  const jlt::mathmatrix<int> id = jlt::identity_matrix<int>(tt.edges());
  int nb = 0;
  for (int f = 0; f < tt.foldings(); ++f)
    {
      traintracks::traintrack t2(tt);
      jlt::freeauto<int> AM = t2.fold_traintrack_map(f);
      if (traintracks::transition_matrix_from_map(tt,AM) == id) continue;
      if (nb == b) return f;
      ++nb;
    }
  return -1;
}

static void check_vertex_fold_consistency(const ttauto::ttfoldgraph<traintracks::traintrack>& ttg,
                                          const int vertex)
{
  using traintracks::fold_map_data;
  using traintracks::mathmatrix_permplus1;
  using traintracks::traintrack;
  using traintracks::transition_matrix_from_map;
  using traintracks::ttnumbering;

  const traintrack& ttv = ttg.traintrack(vertex);
  const int nmain = ttv.edges();

  for (int b = 0; b < ttg.foldings(vertex); ++b)
    {
      const jlt::freeauto<int> AMstep = ttg.traintrack_map(vertex,b);
      const mathmatrix_permplus1 PM = ttg.transition_matrix(vertex,b);

      // Main-edge action agrees with the transition matrix.
      CHECK(transition_matrix_from_map(ttv,AMstep) == PM.full());
      // Branches are real folds, never permutations.
      CHECK(!PM.is_perm());

      // Redo the fold on a copy to get the fold record.
      const int f = fold_index_of_branch(ttv,b);
      CHECK(f >= 0);
      traintrack t2(ttv);
      fold_map_data fm;
      CHECK(t2.fold_with_map(f,fm));
      CHECK(t2 == ttg.traintrack(ttg.target_vertex(vertex,b)));
      const ttnumbering& A = fm.after;

      // The stored map is the fold record's map.
      const jlt::freeauto<int> AMfm = fm.to_freeauto();
      for (int g = 1; g <= nmain + fm.nsides; ++g)
        CHECK(AMstep[g] == AMfm[g]);

      // The record's matrix equals the unit-weight oracle (n folds), the
      // stored matrix, and the map's abelianisation.
      const jlt::mathmatrix<int> oracle = oracles::unit_weight_transition_matrix(ttv,f);
      CHECK(fm.transition_matrix_dense() == oracle);
      CHECK(fm.transition_matrix().full() == oracle);
      CHECK(PM.full() == oracle);
      CHECK(traintracks::fold_transition_matrix(ttv,f).full() == oracle);

      // Exactly one side letter appears in the main-edge images, in the
      // image of the moved edge, and it is a side of the target multigon:
      // it joins the two target prongs.
      int nsideletters = 0;
      for (int g = 1; g <= nmain; ++g)
        for (auto x : AMstep.get_action(g))
          if (A.is_side(x)) { ++nsideletters; CHECK(g == fm.moved + 1); CHECK(x == fm.side); }
      CHECK(nsideletters == 1);
      const int qa = A.tail_of(fm.side), qb = A.head_of(fm.side);
      CHECK(A.prong[qa].multigon == A.prong[fm.target_from].multigon);
      CHECK((qa == fm.target_from && qb == fm.target_to) ||
            (qa == fm.target_to && qb == fm.target_from));
      // The moved edge's image is a continuous path.
      CHECK(fm.moved_word.size() == 3);
      CHECK(A.head_of(fm.moved_word[0]) == A.tail_of(fm.moved_word[1]));
      CHECK(A.head_of(fm.moved_word[1]) == A.tail_of(fm.moved_word[2]));
      // Which copy comes first matches the record.
      CHECK((std::abs(fm.moved_word[0]) == std::abs(fm.edge_image[fm.moved])) == fm.moved_first);
    }
}


int main()
{
  using jlt::freeauto;
  using jlt::freeword;
  using traintracks::build_traintrack_list;
  using traintracks::transition_matrix_from_map;
  using traintracks::traintrack;
  using ttauto::folding_path;
  using ttauto::ttfoldgraph;

  typedef ttfoldgraph<traintrack> ttgraph;
  typedef jlt::vector<traintrack> ttVec;

  {
    // n=3: three punctured monogons, one of them (prong 1) carrying the two
    // edges 1: 0->1 and 2: 1->2 with a cusp between them.  The automaton has
    // a single vertex with two branches.
    const int n = 3;
    const int trk = 0;
    ttVec ttv = build_traintrack_list(n);
    ttgraph ttg(ttv[trk]);
    CHECK(ttg.vertices() == 1);
    CHECK(ttg.foldings(0) == 2);

    for (int f = 0; f < ttv[trk].foldings(); ++f)
      {
        traintracks::multigon* mmc = 0;
        int pc = -1, ec = -1;
        ttv[trk].fold_cusp_location(f,mmc,pc,ec);
        CHECK(mmc != 0 && pc >= 0 && ec >= 0);
      }

    check_vertex_fold_consistency(ttg,0);

    // Hand-checked words.  Branch 0 folds edge 1 onto edge 2: the end of
    // edge 1 slides along edge 2 to the monogon at prong 2 and once around
    // its puncture, which after renumbering is prong 1 (side letter 4).
    // Old edge 1 (0->1, head at the cusp) becomes new edge 1 (0->1), then
    // the loop at prong 1 traversed backwards, then new edge 2 (1->2).
    // Old edge 2 is reversed by the renumbering.
    const freeauto<int> AMf0 = ttg.traintrack_map(0,0);
    CHECK((AMf0[1] == freeword<int>({1,-4,2})));
    CHECK((AMf0[2] == freeword<int>({-2})));
    CHECK((AMf0[3] == freeword<int>({3})));
    CHECK((AMf0[4] == freeword<int>({5})));
    CHECK((AMf0[5] == freeword<int>({4})));

    // Branch 1 folds edge 2 onto edge 1, symmetrically.
    const freeauto<int> AMf1 = ttg.traintrack_map(0,1);
    CHECK((AMf1[1] == freeword<int>({-1})));
    CHECK((AMf1[2] == freeword<int>({1,-4,2})));
    CHECK((AMf1[3] == freeword<int>({4})));
    CHECK((AMf1[4] == freeword<int>({3})));
    CHECK((AMf1[5] == freeword<int>({5})));

    // Composed map for the path [1,0]: a continuous path, and its
    // abelianisation is the path's transition matrix.
    folding_path<traintrack> p(ttg,0);
    p.push_back(1);
    p.push_back(0);
    const jlt::mathmatrix<int> TMp = p.transition_matrix();
    const freeauto<int> AMp = p.traintrack_map();
    CHECK((AMp[1] == freeword<int>({-2,4,-1})));
    CHECK((AMp[2] == freeword<int>({1,-4,2,-5,-2})));
    CHECK(TMp == transition_matrix_from_map(ttv[trk],AMp));
  }

  {
    // Secondary deterministic scenario on a larger automaton.
    const int n = 4;
    const int trk = 1;
    ttVec ttv = build_traintrack_list(n);
    ttgraph ttg(ttv[trk]);

    for (int v = 0; v < ttg.vertices(); ++v)
      check_vertex_fold_consistency(ttg,v);

    folding_path<traintrack> p2(ttg,0);
    p2.push_back(0);
    p2.push_back(1);
    p2.push_back(0);
    CHECK(p2.transition_matrix() == transition_matrix_from_map(ttv[trk],p2.traintrack_map()));
  }

  return 0;
}
