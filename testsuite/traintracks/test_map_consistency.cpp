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

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <jlt/freeauto.hpp>
#include <jlt/stlio.hpp>
#include "traintracks/build.hpp"
#include "traintracks/map.hpp"
#include "traintracks/map_labels.hpp"
#include "traintracks/traintrack.hpp"
#include "ttauto/folding_path.hpp"
#include "ttauto/ttfoldgraph.hpp"


static void check_vertex_fold_consistency(const ttauto::ttfoldgraph<traintracks::traintrack>& ttg,
                                          const int vertex)
{
  using traintracks::mathmatrix_permplus1;
  using traintracks::transition_matrix_from_map;
  using traintracks::ttmap_labeler;

  // For each outgoing fold at this vertex, ensure the fold-level
  // train-track map and transition matrix agree on main-edge action.
  for (int f = 0; f < ttg.foldings(vertex); ++f)
    {
      const traintracks::traintrack& ttvf = ttg.traintrack(vertex);
      const ttmap_labeler labels(ttvf.edges(),ttvf.total_prongs());
      const jlt::freeauto<int> AMstep = ttg.traintrack_map(vertex,f);
      const mathmatrix_permplus1 PM = ttg.transition_matrix(vertex,f);

      assert(transition_matrix_from_map(ttvf,AMstep) == PM.full());

      // In non-permutation folds, exactly one positively oriented selected
      // infinitesimal generator should be injected into main-edge images.
      const int infg = std::abs(ttvf.fold_infinitesimal_generator(f,labels.nmain));
      int ninf = 0;
      int nneg = 0;
      int npos = 0;

      for (int g = 1; g <= labels.nmain; ++g)
        {
          for (auto img : AMstep.get_action(g))
            {
              if (img == -infg) { ++ninf; ++nneg; }
              else if (img == infg) { ++ninf; ++npos; }
            }
        }

      if (PM.is_perm())
        {
          assert(ninf == 0);
        }
      else
        {
          assert(ninf == 1);
          assert(npos == 1);
          assert(nneg == 0);
        }
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
    // Small hand-checked scenario used in issue #3 development.
    const int n = 3;
    const int trk = 0;
    ttVec ttv = build_traintrack_list(n);
    ttgraph ttg(ttv[trk]);

    for (int f = 0; f < ttv[trk].foldings(); ++f)
      {
        // Cusp-to-infinitesimal mapping must resolve to in-range indices.
        [[maybe_unused]] const int infix = ttv[trk].fold_infinitesimal_index(f);
        assert(infix >= 0);
        assert(infix < ttv[trk].total_prongs());

        traintracks::multigon* mmc = 0;
        int pc = -1, ec = -1;
        ttv[trk].fold_cusp_location(f,mmc,pc,ec);
        assert(mmc != 0);
        assert(pc >= 0);
        assert(ec >= 0);
      }

    check_vertex_fold_consistency(ttg,0);

    // Lock a known one-step map from the hand-worked example.
    const freeauto<int> AMf1 = ttg.traintrack_map(0,1);
    assert((AMf1[1] == freeword<int>({-1})));
    assert((AMf1[2] == freeword<int>({1,5,2})));

    const int v_after_f1 = ttg.target_vertex(0,1);
    const freeauto<int> AMf0_after_f1 = ttg.traintrack_map(v_after_f1,0);
    assert((AMf0_after_f1[1] == freeword<int>({1,5,2})));
    assert((AMf0_after_f1[2] == freeword<int>({-2})));

    // Lock composed map and matrix agreement for path [1,0].
    folding_path<traintrack> p(ttg,0);
    p.push_back(1);
    p.push_back(0);

    const jlt::mathmatrix<int> TMp = p.transition_matrix();
    const freeauto<int> AMp = p.traintrack_map();
    const jlt::mathmatrix<int> TMfromAMp = transition_matrix_from_map(ttv[trk],AMp);

    assert((AMp[1] == freeword<int>({-2,-5,-1})));
    assert((AMp[2] == freeword<int>({1,5,2,5,-2})));
    assert(TMp == TMfromAMp);
  }

  {
    // Secondary deterministic scenario on a larger automaton.
    const int n = 4;
    const int trk = 1;
    ttVec ttv = build_traintrack_list(n);
    ttgraph ttg(ttv[trk]);

    const int max_vertices_to_check = std::min(3,ttg.vertices());
    for (int v = 0; v < max_vertices_to_check; ++v)
      {
        check_vertex_fold_consistency(ttg,v);
      }

    folding_path<traintrack> p2(ttg,0);
    p2.push_back(0);
    p2.push_back(1);
    p2.push_back(0);

    const jlt::mathmatrix<int> TMp2 = p2.transition_matrix();
    const freeauto<int> AMp2 = p2.traintrack_map();
    const jlt::mathmatrix<int> TMfromAMp2 = transition_matrix_from_map(ttv[trk],AMp2);

    assert(TMp2 == TMfromAMp2);
  }

  return 0;
}
