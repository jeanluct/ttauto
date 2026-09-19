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

// Exploratory companion of testsuite/traintracks/test_map_consistency.cpp:
// prints the canonical numbering, every one-step fold record and the
// composed map of a short path, and asserts the same consistency checks.

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <jlt/freeauto.hpp>
#include "traintracks/build.hpp"
#include "traintracks/fold_map.hpp"
#include "traintracks/map.hpp"
#include "traintracks/traintrack.hpp"
#include "ttauto/folding_path.hpp"
#include "ttauto/ttfoldgraph.hpp"

using jlt::freeauto;
using traintracks::fold_map_data;
using traintracks::traintrack;
using traintracks::ttnumbering;
using ttauto::folding_path;
using ttauto::ttfoldgraph;

static int fold_index_of_branch(const traintrack& tt, const int b)
{
  const jlt::mathmatrix<int> id = jlt::identity_matrix<int>(tt.edges());
  int nb = 0;
  for (int f = 0; f < tt.foldings(); ++f)
    {
      traintrack t2(tt);
      freeauto<int> AM = t2.fold_traintrack_map(f);
      if (traintracks::transition_matrix_from_map(tt,AM) == id) continue;
      if (nb == b) return f;
      ++nb;
    }
  return -1;
}

static void explore(const int n, const int trk, const int maxvertices)
{
  std::cout << "\n==== n=" << n << ", track " << trk << "\n";
  jlt::vector<traintrack> ttv = traintracks::build_traintrack_list(n);
  ttfoldgraph<traintrack> ttg(ttv[trk]);
  std::cout << "Folding automaton has " << ttg.vertices() << " vertices\n";

  for (int v = 0; v < std::min(maxvertices,ttg.vertices()); ++v)
    {
      const traintrack& tt = ttg.traintrack(v);
      std::cout << "\n-- vertex " << v << "\n";
      tt.print(std::cout);
      tt.numbering().print(std::cout);
      for (int b = 0; b < ttg.foldings(v); ++b)
        {
          const int f = fold_index_of_branch(tt,b);
          assert(f >= 0);
          traintrack t2(tt);
          fold_map_data fm;
          if (!t2.fold_with_map(f,fm)) { std::cerr << "fold failed\n"; std::exit(1); }
          std::cout << "branch " << b << " (fold index " << f << ") -> vertex "
                    << ttg.target_vertex(v,b) << "\n";
          fm.print(std::cout);
          const freeauto<int> AM = ttg.traintrack_map(v,b);
          std::cout << AM;
          // Consistency: stored map equals the record's map; matrix agrees.
          const freeauto<int> AMfm = fm.to_freeauto();
          for (int g = 1; g <= fm.nmain + fm.nsides; ++g) assert(AM[g] == AMfm[g]);
          assert(traintracks::transition_matrix_from_map(tt,AM)
                 == ttg.transition_matrix(v,b).full());
          // Continuity of every image in the target's numbering.
          const ttnumbering Nt = ttg.traintrack(ttg.target_vertex(v,b)).numbering();
          for (int g = 1; g <= fm.nmain; ++g)
            {
              const jlt::freeword<int> w = AM.get_action(g);
              int prev = 0;
              for (auto x : w)
                {
                  if (prev != 0) assert(Nt.head_of(prev) == Nt.tail_of(x));
                  prev = x;
                }
            }
        }
    }

  // A short path from vertex 0.
  folding_path<traintrack> p(ttg,0);
  for (int i = 0; i < 3; ++i) p.push_back(0);
  std::cout << "\nPath " << p << " composed map:\n" << p.traintrack_map();
  assert(p.transition_matrix()
         == traintracks::transition_matrix_from_map(ttv[trk],p.traintrack_map()));
}

int main()
{
  explore(3,0,10);
  explore(4,1,3);
  std::cout << "\ntest_ttmap: OK\n";
  return 0;
}
