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

#include <cassert>
#include <iostream>
#include <list>

#include <jlt/freeauto.hpp>
#include <jlt/mathmatrix.hpp>
#include <jlt/vector.hpp>

#include "traintracks/build.hpp"
#include "traintracks/map.hpp"
#include "traintracks/traintrack.hpp"
#include "ttauto/folding_path.hpp"
#include "ttauto/ttfoldgraph.hpp"


int main()
{
  using traintracks::traintrack;
  using traintracks::build_traintrack_list;
  using ttauto::ttfoldgraph;
  using ttauto::folding_path;

  typedef ttfoldgraph<traintrack> ttgraph;

  const int n = 6;
  jlt::vector<traintrack> ttv = build_traintrack_list(n);

  // Issue #2 repro data from devel/iss002/n=6_5_bad_tt_data.m
  // Vertex cycle listed there is 1-based.
  const int bad_cycle_1based[] = {29, 46, 43, 71, 88, 85, 29};
  const int cycle_len = sizeof(bad_cycle_1based)/sizeof(bad_cycle_1based[0]);
  // Hardwired branch sequence for the bad example (0-based branch ids at
  // each step in the cycle above).
  const int bad_branches[] = {1, 0, 3, 2, 1, 2};

  // Hardwire the known bad graph: n=6, traintrack index 4, first pruned
  // subgraph (90 vertices) from devel/iss002/n=6_5_bad_tt_data.m.
  const int trk = 4;
  const int sgidx = 0;
  ttgraph full(ttv[trk]);
  std::list<ttgraph> sgs = ttauto::subgraphs(full);
  ttauto::prune_multihumps(sgs);
  assert((int)sgs.size() > sgidx);
  auto it = sgs.begin();
  std::advance(it,sgidx);
  const ttgraph& ttg = *it;

  std::cout << "Issue #2 bad-path reproducer\n";
  std::cout << "n=" << n << ", track index=" << trk
            << ", subgraph index=" << sgidx
            << ", vertices=" << ttg.vertices() << "\n";
  assert(ttg.vertices() == 90);

  int bad_cycle_0based[cycle_len];
  for (int i = 0; i < cycle_len; ++i)
    {
      bad_cycle_0based[i] = bad_cycle_1based[i]-1;
      assert(bad_cycle_0based[i] >= 0);
      assert(bad_cycle_0based[i] < ttg.vertices());
    }

  const int vstart = bad_cycle_0based[0];

  std::cout << "vertex cycle (1-based): ";
  for (int i = 0; i < cycle_len; ++i)
    {
      std::cout << bad_cycle_1based[i];
      if (i+1 < cycle_len) std::cout << " -> ";
    }
  std::cout << "\n";

  std::cout << "branch sequence (0-based): ";
  for (int i = 0; i < cycle_len-1; ++i)
    {
      const int v = bad_cycle_0based[i];
      const int vnext = bad_cycle_0based[i+1];

      const int b = bad_branches[i];
      assert(b >= 0 && b < ttg.foldings(v));
      assert(ttg.target_vertex(v,b) == vnext);
      if (i > 0) std::cout << ", ";
      std::cout << b;

      if (i == 0)
        {
          // Non-essential sanity print to make it obvious this is the
          // hardwired known-bad realization.
          std::cout << " [hardwired]";
        }
    }
  std::cout << "\n";

  folding_path<traintrack> p(ttg,vstart);
  for (int i = 0; i < cycle_len-1; ++i)
    {
      p.push_back(bad_branches[i]);
    }

  jlt::freeauto<int> AM = p.traintrack_map();
  jlt::mathmatrix<int> TM = p.transition_matrix();
  jlt::mathmatrix<int> TMfromAM =
    traintracks::transition_matrix_from_map(ttg.traintrack(vstart),AM);
  assert(TM == TMfromAM);

  std::cout << "\nComposed train-track map:\n";
  std::cout << AM;
  std::cout << "\nComposed transition matrix (path):\n";
  TM.printMatrixForm(std::cout);

  const bool primitive = TM.is_primitive();
  std::cout << "\nPrimitive transition matrix: "
            << (primitive ? "yes" : "no") << "\n";
  assert(primitive);

  std::cout << "\nTransition matrix from composed map:\n";
  TMfromAM.printMatrixForm(std::cout);
  std::cout << "\nOK: matrix(path) == matrix(from map).\n";

  return 0;
}
