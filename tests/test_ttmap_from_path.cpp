// <LICENSE
//   ttauto: a C++ library for building train track automata
//
//   https://github.com/jeanluct/ttauto
//
//   Copyright (C) 2010-2014  Jean-Luc Thiffeault   <jeanluc@math.wisc.edu>
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

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>

#include "folding_path.hpp"
#include "traintrack_build.hpp"
#include "traintrack.hpp"
#include "traintrack_map.hpp"
#include "ttfoldgraph.hpp"


template<class TrTr>
void check_path_map_matrix(const traintracks::ttfoldgraph<TrTr>& ttg,
                           const TrTr& tt,
                           const int initial_vertex,
                           const std::initializer_list<int>& folds,
                           const char* label)
{
  std::cout << "\n[ttmap-check] " << label << "\n";

  // Build an explicit folding path in the automaton.
  traintracks::folding_path<TrTr> p(ttg,initial_vertex);
  for (int f : folds) p.push_back(f);

  std::cout << "  folds:";
  for (int f : folds) std::cout << " " << f;
  std::cout << "\n";

  // Matrix tracked directly along the path.
  jlt::mathmatrix<int> TMpath = p.transition_matrix();

  // Train-track map (automorphism on main+infinitesimal generators) derived
  // by composing one-step fold maps along the same path.
  traintracks::freeauto<int> AMpath = p.traintrack_map();
  std::cout << "  train-track map:\n";
  std::cout << AMpath;

  // Project automorphism back to the main-edge transition matrix and require
  // exact agreement.
  jlt::mathmatrix<int> TMfromAM = traintracks::transition_matrix_from_map(tt,AMpath);

  std::cout << "  transition matrix (path):\n";
  TMpath.printMatrixForm(std::cout);
  std::cout << "  transition matrix (from map):\n";
  TMfromAM.printMatrixForm(std::cout);

  assert(TMpath == TMfromAM);
  std::cout << "  -> OK\n";
}


template<class TrTr>
void stress_random_paths(const traintracks::ttfoldgraph<TrTr>& ttg,
                         const TrTr& tt,
                         const int initial_vertex,
                         const int npaths,
                         const int max_len,
                         std::mt19937& rng,
                         const std::string& label)
{
  std::uniform_int_distribution<int> len_dist(1,max_len);

  for (int k = 0; k < npaths; ++k)
    {
      traintracks::folding_path<TrTr> p(ttg,initial_vertex);
      int v = initial_vertex;
      const int L = len_dist(rng);

      for (int step = 0; step < L; ++step)
        {
          const int nf = ttg.foldings(v);
          if (nf <= 0) break;
          std::uniform_int_distribution<int> fold_dist(0,nf-1);
          const int f = fold_dist(rng);
          p.push_back(f);
          v = ttg.target_vertex(v,f);
        }

      jlt::mathmatrix<int> TMpath = p.transition_matrix();
      traintracks::freeauto<int> AMpath = p.traintrack_map();
      jlt::mathmatrix<int> TMfromAM = traintracks::transition_matrix_from_map(tt,AMpath);

      assert(TMpath == TMfromAM);
    }

  std::cout << "[ttmap-random] " << label << " -> OK"
            << " (paths=" << npaths << ", max_len=" << max_len << ")\n";
}


int main(int argc, char** argv)
{
  using traintracks::traintrack;

  typedef traintracks::ttfoldgraph<traintrack> ttgraph;
  typedef jlt::vector<traintrack> ttVec;

  {
    // Small hand-checked scenario used in issue #3 work.
    const int n = 3;
    ttVec ttv = traintracks::build_traintrack_list(n);
    const int trk = 0;
    ttgraph ttg(ttv[trk]);

    // Extended from the fixed two-step example in TTMAP_EXAMPLE.md.
    check_path_map_matrix(ttg,ttv[trk],0,{1,0,1,0},
                          "n=3, trk=0, path [1,0,1,0]");
  }

  {
    // Secondary scenario on a larger automaton.
    const int n = 4;
    ttVec ttv = traintracks::build_traintrack_list(n);
    const int trk = 1;
    ttgraph ttg(ttv[trk]);

    // Extended nontrivial path used for consistency stress.
    check_path_map_matrix(ttg,ttv[trk],0,{0,1,0,1,0},
                          "n=4, trk=1, path [0,1,0,1,0]");
  }

  {
    // Randomized stress over all tracks for n=3..7 (path choice random only).
    // Optional argv[1] sets the RNG seed. Default keeps runs reproducible.
    unsigned seed = 123456u;
    if (argc > 1)
      {
        seed = static_cast<unsigned>(std::strtoul(argv[1],0,10));
      }
    std::mt19937 rng(seed);

    std::cout << "\n";

    for (int n = 3; n <= 7; ++n)
      {
        ttVec ttv = traintracks::build_traintrack_list(n);
        for (int trk = 0; trk < (int)ttv.size(); ++trk)
          {
            ttgraph ttg(ttv[trk]);
            std::string label = "n=" + std::to_string(n) +
                                ", trk=" + std::to_string(trk) +
                                " random paths";
            stress_random_paths(ttg,ttv[trk],0,20,10,rng,label);
          }
      }
  }

  std::cout << "\ntest_ttmap_from_path: ok\n";
}
