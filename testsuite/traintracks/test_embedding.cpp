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

// The boundary walk of a train track (traintracks/embedding.hpp): the
// counting identities that say the walk really is the boundary region, and
// the puncture order of the two tracks drawn by hand in doc/ttauto.tex.

#include <iostream>
#include <set>
#include <vector>
#include "check.hpp"
#include "traintracks/build.hpp"
#include "traintracks/embedding.hpp"
#include "traintracks/traintrack.hpp"
#include "ttauto/ttfoldgraph.hpp"

using traintracks::traintrack;
using traintracks::tt_embedding;
using traintracks::ttnumbering;
using traintracks::outer_embedding;

static void check_walk(const traintrack& tt, const char* what)
{
  const ttnumbering N = tt.numbering();
  const tt_embedding emb = outer_embedding(N);

  // Every main edge is a bridge, since the multigons and main edges form a
  // tree, so the boundary region runs along both its sides; every side has
  // a complementary region on one side and the boundary region on the
  // other, so the walk meets it once.
  CHECK_MSG((int)emb.walk.size() == 2*tt.edges() + tt.total_prongs(),what);

  // Each main letter twice, each positive side letter once, no negative
  // side letter: those bound the polygons and the punctures.
  std::vector<int> mains(tt.edges(),0), sides(tt.total_prongs(),0);
  for (int i = 0; i < (int)emb.walk.size(); ++i)
    {
      const int d = emb.walk[i];
      if (N.is_main(d)) ++mains[N.edge_of(d)];
      else { CHECK_MSG(d > 0,what); ++sides[N.side_of(d)]; }
    }
  for (int e = 0; e < tt.edges(); ++e) CHECK_MSG(mains[e] == 2,what);
  for (int q = 0; q < tt.total_prongs(); ++q) CHECK_MSG(sides[q] == 1,what);

  // The corners between two main letters are the exterior cusps.
  CHECK_MSG((int)emb.cusp_after.size() == tt.cusps(),what);

  // Every puncture once, and position_of is its inverse.
  CHECK_MSG(emb.npunctures() == tt.punctures(),what);
  std::set<int> seen(emb.puncture_order.begin(),emb.puncture_order.end());
  CHECK_MSG((int)seen.size() == tt.punctures(),what);
  for (int i = 0; i < emb.npunctures(); ++i)
    CHECK_MSG(emb.position_of[emb.puncture_order[i]] == i+1,what);

  // Cutting the walk anywhere gives the same cyclic order of punctures.
  for (int i = 1; i < emb.npunctures(); ++i)
    {
      const tt_embedding e2
        = outer_embedding(N,N.side_letter(emb.puncture_order[i]));
      for (int j = 0; j < emb.npunctures(); ++j)
        CHECK_MSG(e2.puncture_order[j]
                  == emb.puncture_order[(i+j) % emb.npunctures()],what);
    }
}

int main()
{
  long nvertices = 0;

  // Every vertex of every automaton up to six punctures.
  for (int n = 3; n <= 6; ++n)
    {
      jlt::vector<traintrack> ttv = traintracks::build_traintrack_list(n);
      for (int s = 0; s < (int)ttv.size(); ++s)
        {
          ttauto::ttfoldgraph<traintrack> ttg(ttv[s]);
          for (int v = 0; v < ttg.vertices(); ++v)
            {
              check_walk(ttg.traintrack(v),"automaton vertex");
              ++nvertices;
            }
        }
    }

  // The two tracks of doc/ttauto.tex Fig. 2, drawn by hand with the
  // punctures in place.  Reading the pictures anticlockwise, which is
  // left to right along the real axis, the four punctures of the stratum
  // 1.1.1.1.(2) come in the order prong 0, 2, 3, 1, and those of
  // 1.1.1.1.3(1) in the order prong 0, 5, 6, 1.  This is what fixes the
  // sense of the walk, and so the handedness of every braid read off a
  // folding path.
  {
    jlt::vector<traintrack> ttv = traintracks::build_traintrack_list(4);
    const int want1[4] = {0,2,3,1};
    const int want2[4] = {0,5,6,1};
    const tt_embedding e1 = outer_embedding(ttv[0].numbering());
    const tt_embedding e2 = outer_embedding(ttv[1].numbering());
    for (int i = 0; i < 4; ++i)
      {
        CHECK_MSG(e1.puncture_order[i] == want1[i],"1.1.1.1.(2) by hand");
        CHECK_MSG(e2.puncture_order[i] == want2[i],"1.1.1.1.3(1) by hand");
      }
  }

  std::cout << "\ntest_embedding: boundary walk checked at " << nvertices
            << " automaton vertices\n";
  return 0;
}
