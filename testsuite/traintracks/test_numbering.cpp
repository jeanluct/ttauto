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

// Canonical prong/edge numbering (ttnumbering): structural invariants and
// agreement with the edge order used by weights().

#include <iostream>
#include <set>
#include <vector>
#include "check.hpp"
#include "traintracks/build.hpp"
#include "traintracks/traintrack.hpp"

using traintracks::traintrack;
using traintracks::ttnumbering;

static void check_numbering(const traintrack& tt, const char* what)
{
  const ttnumbering N = tt.numbering();

  CHECK_MSG(N.nprongs() == tt.total_prongs(), what);
  CHECK_MSG(N.nedges() == tt.edges(), what);
  CHECK_MSG(N.start_monogon == 0, what);
  CHECK_MSG(N.prong[0].multigon == 0 && N.prong[0].prong == 0, what);

  // prong_number is a bijection onto 0..nprongs-1 consistent with prong[].
  std::set<int> seen;
  for (int m = 0; m < tt.multigons(); ++m)
    {
      CHECK_MSG((int)N.prong_number[m].size() == tt.Multigon(m).prongs(), what);
      for (int p = 0; p < tt.Multigon(m).prongs(); ++p)
        {
          const int q = N.prong_number[m][p];
          CHECK_MSG(q >= 0 && q < N.nprongs(), what);
          CHECK_MSG(seen.insert(q).second, what);
          CHECK_MSG(N.prong[q].multigon == m && N.prong[q].prong == p, what);
          CHECK_MSG(N.prong[q].nprongs == tt.Multigon(m).prongs(), what);
          CHECK_MSG(N.prong[q].punctured == tt.Multigon(m).punctured(), what);
        }
    }

  // Sides stay inside their multigon; prongs of one multigon are numbered
  // consecutively in cycle order; monogon sides are loops.
  for (int q = 0; q < N.nprongs(); ++q)
    {
      const int q2 = N.side_to(q);
      CHECK_MSG(N.prong[q2].multigon == N.prong[q].multigon, what);
      if (N.prong[q].nprongs == 1) CHECK_MSG(q2 == q, what);
      else CHECK_MSG(q2 != q, what);
    }

  // Edges: valid distinct endpoints; letter helpers consistent.
  for (int e = 0; e < N.nedges(); ++e)
    {
      const int t = N.edge_tail[e], h = N.edge_head[e];
      CHECK_MSG(t >= 0 && t < N.nprongs() && h >= 0 && h < N.nprongs(), what);
      CHECK_MSG(N.prong[t].multigon != N.prong[h].multigon, what);
      CHECK_MSG(N.is_main(N.main_letter(e)) && !N.is_side(N.main_letter(e)), what);
      CHECK_MSG(N.tail_of(N.main_letter(e)) == t && N.head_of(N.main_letter(e)) == h, what);
      CHECK_MSG(N.tail_of(-N.main_letter(e)) == h && N.head_of(-N.main_letter(e)) == t, what);
      CHECK_MSG(N.edge_of(N.main_letter(e)) == e, what);
      // The edge object attached at the tail prong slot list contains it.
      const int m = N.prong[t].multigon, p = N.prong[t].prong;
      bool found = false;
      for (int s = 0; s < tt.Multigon(m).edges(p); ++s)
        if (tt.Multigon(m).Edge(p,s).get() == N.edge_ptr[e]) found = true;
      CHECK_MSG(found, what);
    }
  for (int q = 0; q < N.nprongs(); ++q)
    {
      CHECK_MSG(N.is_side(N.side_letter(q)) && !N.is_main(N.side_letter(q)), what);
      CHECK_MSG(N.side_of(N.side_letter(q)) == q, what);
      CHECK_MSG(N.tail_of(N.side_letter(q)) == q, what);
      CHECK_MSG(N.head_of(N.side_letter(q)) == N.side_to(q), what);
    }

  // Edge order agrees with weights(): set weight e+1 on the e-th weight
  // slot and read it back through the edge objects.
  {
    traintrack ttw(tt);
    const ttnumbering Nw = ttw.numbering();
    traintrack::dblVec wv(Nw.nedges());
    for (int e = 0; e < Nw.nedges(); ++e) wv[e] = e + 1;
    ttw.weights(wv.begin());
    for (int e = 0; e < Nw.nedges(); ++e)
      CHECK_MSG(Nw.edge_ptr[e]->weight() == e + 1, what);
  }

  // A copy has the same numbering, and so does a re-normalised copy.
  traintrack tt2(tt);
  CHECK_MSG(tt2.numbering() == N, what);
  tt2.normalise();
  CHECK_MSG(tt2.numbering() == N, what);
}

int main()
{
  int count = 0;
  for (int n = 3; n <= 5; ++n)
    {
      jlt::vector<traintrack> ttv = traintracks::build_traintrack_list(n);
      for (std::size_t trk = 0; trk < ttv.size(); ++trk)
        {
          check_numbering(ttv[trk],"initial track");
          // Also every neighbour reached by one fold.
          for (int f = 0; f < ttv[trk].foldings(); ++f)
            {
              traintrack t2(ttv[trk]);
              if (t2.fold(f)) { check_numbering(t2,"folded track"); ++count; }
            }
        }
    }
  CHECK(count > 0);
  std::cout << "test_numbering: OK (" << count << " folded tracks checked)\n";
  return 0;
}
