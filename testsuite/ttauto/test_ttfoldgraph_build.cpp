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

// Automaton construction (issue #14): it must not recurse once per
// vertex, the codings that index the vertices must be unique, and the
// graphs must be exactly the ones recorded here.

#include <cstddef>
#include <iostream>
#include <set>
#include <jlt/vector.hpp>
#include "traintracks/build.hpp"
#include "traintracks/traintrack.hpp"
#include "ttauto/ttfoldgraph.hpp"
#include "check.hpp"

#ifdef TTAUTO_HAVE_PTHREAD
#include <pthread.h>
#endif

using traintracks::traintrack;
typedef ttauto::ttfoldgraph<traintrack> ttgraph;

namespace {

// Vertex and branch counts of every automaton for n = 3..6.  These pin
// the construction order as well as the graph: the search, the gate
// test and the papers all quote vertex numbers, so a change here is a
// change to published material, not a test to be updated lightly.
struct expected { int n; int trk; int vertices; int branches; };

const expected graphs[] = {
  {3,0,1,2},
  {4,0,4,14}, {4,1,3,6},
  {5,0,11,50}, {5,1,21,76}, {5,2,3,6}, {5,3,9,18},
  {6,0,49,282}, {6,1,138,700}, {6,2,27,100}, {6,3,3,6},
  {6,4,110,412}, {6,5,21,42}, {6,6,28,56},
};

#ifdef TTAUTO_HAVE_PTHREAD
// The largest seven-puncture stratum, built on a thread of small stack.
const int small_stack_trk = 5;
const int small_stack_vertices = 1012;
int built_vertices = -1;

void* build_on_small_stack(void*)
{
  jlt::vector<traintrack> ttv = traintracks::build_traintrack_list(7);
  ttgraph ttg(ttv[small_stack_trk]);
  built_vertices = ttg.vertices();
  return 0;
}
#endif

} // namespace


int main()
{
#ifdef TTAUTO_HAVE_PTHREAD
  {
    // Construction must not use stack proportional to the vertex count.
    // Before the worklist replaced the recursion, this stratum needed
    // more than 384 KB; the loop builds it in under 32 KB.  Run it on a
    // 256 KB thread, which discriminates between the two with room to
    // spare on either side.
    const std::size_t stacksize = 256*1024;
    pthread_attr_t attr;
    CHECK(pthread_attr_init(&attr) == 0);
    CHECK(pthread_attr_setstacksize(&attr,stacksize) == 0);
    pthread_t th;
    CHECK(pthread_create(&th,&attr,build_on_small_stack,0) == 0);
    CHECK(pthread_join(th,0) == 0);
    pthread_attr_destroy(&attr);
    CHECK_MSG(built_vertices == small_stack_vertices,
              "n=7 stratum " << small_stack_trk+1 << " built "
              << built_vertices << " vertices");
    std::cout << "built n=7 stratum " << small_stack_trk+1 << " ("
              << built_vertices << " vertices) on a "
              << stacksize/1024 << " KB stack\n";
  }
#else
  std::cout << "pthread unavailable: skipping the small-stack build\n";
#endif

  int ngraphs = 0, nvertices = 0;
  for (const expected& e : graphs)
    {
      jlt::vector<traintrack> ttv = traintracks::build_traintrack_list(e.n);
      CHECK((int)ttv.size() > e.trk);
      ttgraph ttg(ttv[e.trk]);

      CHECK_MSG(ttg.vertices() == e.vertices,
                "n=" << e.n << " stratum " << e.trk+1 << ": "
                << ttg.vertices() << " vertices, expected " << e.vertices);

      int branches = 0;
      std::set<traintrack::intVec> codings;
      for (int v = 0; v < ttg.vertices(); ++v)
        {
          branches += ttg.foldings(v);
          codings.insert(ttg.traintrack(v).coding());
          for (int b = 0; b < ttg.foldings(v); ++b)
            {
              const int t = ttg.target_vertex(v,b);
              CHECK_MSG(t >= 0 && t < ttg.vertices(),
                        "n=" << e.n << " stratum " << e.trk+1
                        << ": branch " << b << " of vertex " << v
                        << " targets " << t);
            }
        }
      CHECK_MSG(branches == e.branches,
                "n=" << e.n << " stratum " << e.trk+1 << ": "
                << branches << " branches, expected " << e.branches);

      // Vertices are indexed by their coding while the graph is built,
      // so two vertices sharing a coding would be silently merged.  A
      // failure here would also mean issue #12 had returned.
      CHECK_MSG((int)codings.size() == ttg.vertices(),
                "n=" << e.n << " stratum " << e.trk+1 << ": only "
                << codings.size() << " distinct codings for "
                << ttg.vertices() << " vertices");

      ++ngraphs;
      nvertices += ttg.vertices();
    }
  CHECK(ngraphs > 0);

  std::cout << "test_ttfoldgraph_build: OK (" << ngraphs << " automata, "
            << nvertices << " vertices, codings distinct)\n";
  return 0;
}
