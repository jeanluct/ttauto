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
#include <iostream>

#include "folding_path.hpp"
#include "traintrack.hpp"
#include "traintracks_util.hpp"
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

  // Matrix tracked directly along the path (library convention uses transpose
  // for reporting/printing path transition matrices in this context).
  jlt::mathmatrix<int> TMpath = p.transition_matrix().transpose();

  // Train-track map (automorphism on main+infinitesimal generators) derived
  // by composing one-step fold maps along the same path.
  traintracks::free_auto<int> AMpath = p.traintrack_map();
  std::cout << "  train-track map:\n";
  std::cout << AMpath;

  // Project automorphism back to the main-edge transition matrix in the same
  // transposed convention and require exact agreement.
  jlt::mathmatrix<int> TMfromAM =
    traintracks::transition_matrix_from_map_transposed(tt,AMpath);

  std::cout << "  transition matrix (path):\n";
  TMpath.printMatrixForm(std::cout);
  std::cout << "  transition matrix (from map):\n";
  TMfromAM.printMatrixForm(std::cout);

  assert(TMpath == TMfromAM);
  std::cout << "  -> OK\n";
}


int main()
{
  using traintracks::traintrack;

  typedef traintracks::ttfoldgraph<traintrack> ttgraph;
  typedef jlt::vector<traintrack> ttVec;

  {
    // Small hand-checked scenario used in issue #3 work.
    const int n = 3;
    ttVec ttv = traintracks::ttbuild_list(n);
    const int trk = 0;
    ttgraph ttg(ttv[trk]);

    // Extended from the fixed two-step example in TTMAP_EXAMPLE.md.
    check_path_map_matrix(ttg,ttv[trk],0,{1,0,1,0},
                          "n=3, trk=0, path [1,0,1,0]");
  }

  {
    // Secondary scenario on a larger automaton.
    const int n = 4;
    ttVec ttv = traintracks::ttbuild_list(n);
    const int trk = 1;
    ttgraph ttg(ttv[trk]);

    // Extended nontrivial path used for consistency stress.
    check_path_map_matrix(ttg,ttv[trk],0,{0,1,0,1,0},
                          "n=4, trk=1, path [0,1,0,1,0]");
  }

  std::cout << "test_ttmap_from_path: ok\n";
}
