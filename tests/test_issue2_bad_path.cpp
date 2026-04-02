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
#include <memory>
#include <sstream>
#include <vector>

#include <jlt/freeauto.hpp>
#include <jlt/mathmatrix.hpp>

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

  std::unique_ptr<ttgraph> ttg_ptr;
  int trk = -1;
  int sgidx = -1;

  for (int i = 0; i < (int)ttv.size(); ++i)
    {
      ttgraph full(ttv[i]);
      std::list<ttgraph> sgs = ttauto::subgraphs(full);
      ttauto::prune_multihumps(sgs);

      int j = 0;
      for (auto it = sgs.begin(); it != sgs.end(); ++it, ++j)
        {
          if (it->vertices() != 90) continue;

          bool cycle_ok = true;
          for (int k = 0; k < cycle_len; ++k)
            {
              int v = bad_cycle_1based[k]-1;
              if (v < 0 || v >= it->vertices())
                {
                  cycle_ok = false;
                  break;
                }
            }

          for (int k = 0; cycle_ok && k < cycle_len-1; ++k)
            {
              int v = bad_cycle_1based[k]-1;
              int vnext = bad_cycle_1based[k+1]-1;
              bool has_branch = false;
              for (int b = 0; b < it->foldings(v); ++b)
                {
                  if (it->target_vertex(v,b) == vnext)
                    {
                      has_branch = true;
                      break;
                    }
                }
              if (!has_branch) cycle_ok = false;
            }

          if (cycle_ok)
            {
              trk = i;
              sgidx = j;
              ttg_ptr.reset(new ttgraph(*it));
              break;
            }
        }

      if (ttg_ptr) break;
    }

  if (!ttg_ptr)
    {
      std::cerr << "Could not locate expected 90-vertex bad-stratum graph/cycle.\n";
      return 1;
    }

  const ttgraph& ttg = *ttg_ptr;

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

  std::vector<std::vector<int> > branch_options;
  int total_paths = 1;

  std::cout << "branch options per step (0-based):\n";
  for (int i = 0; i < cycle_len-1; ++i)
    {
      const int v = bad_cycle_0based[i];
      const int vnext = bad_cycle_0based[i+1];

      std::vector<int> options;
      for (int b = 0; b < ttg.foldings(v); ++b)
        {
          if (ttg.target_vertex(v,b) == vnext)
            {
              options.push_back(b);
            }
        }

      assert(!options.empty());
      branch_options.push_back(options);
      total_paths *= options.size();

      std::cout << "  step " << i << " (" << bad_cycle_1based[i]
                << " -> " << bad_cycle_1based[i+1] << "): ";
      for (int j = 0; j < (int)options.size(); ++j)
        {
          std::cout << options[j];
          if (j+1 < (int)options.size()) std::cout << ", ";
        }
      std::cout << "\n";
    }

  std::cout << "candidate branch-realizations: " << total_paths << "\n";

  bool all_tm_equal = true;
  bool all_map_strings_equal = true;
  bool first = true;
  jlt::mathmatrix<int> TM0;
  std::string AM0;

  std::vector<int> idx(branch_options.size(),0);
  for (int case_id = 1; ; ++case_id)
    {
      folding_path<traintrack> pcase(ttg,vstart);
      std::vector<int> branches;
      branches.reserve(branch_options.size());
      for (int i = 0; i < (int)branch_options.size(); ++i)
        {
          const int b = branch_options[i][idx[i]];
          branches.push_back(b);
          pcase.push_back(b);
        }

      jlt::freeauto<int> AM = pcase.traintrack_map();
      jlt::mathmatrix<int> TM = pcase.transition_matrix();
      jlt::mathmatrix<int> TMfromAM =
        traintracks::transition_matrix_from_map(ttg.traintrack(vstart),AM);
      assert(TM == TMfromAM);

      std::ostringstream oss;
      oss << AM;
      std::string AMstr = oss.str();

      if (first)
        {
          TM0 = TM;
          AM0 = AMstr;
          first = false;
        }
      else
        {
          if (!(TM == TM0)) all_tm_equal = false;
          if (AMstr != AM0) all_map_strings_equal = false;
        }

      std::cout << "\nCase " << case_id << " branch sequence: ";
      for (int i = 0; i < (int)branches.size(); ++i)
        {
          std::cout << branches[i];
          if (i+1 < (int)branches.size()) std::cout << ", ";
        }

      std::cout << "\nComposed train-track map:\n";
      std::cout << AM;
      std::cout << "\nComposed transition matrix (path):\n";
      TM.printMatrixForm(std::cout);

      int k = (int)idx.size()-1;
      while (k >= 0)
        {
          ++idx[k];
          if (idx[k] < (int)branch_options[k].size()) break;
          idx[k] = 0;
          --k;
        }
      if (k < 0) break;
    }

  std::cout << "\nAll candidate paths agree on transition matrix: "
            << (all_tm_equal ? "yes" : "no") << "\n";
  std::cout << "All candidate paths agree on full train-track map printout: "
            << (all_map_strings_equal ? "yes" : "no") << "\n";

  return 0;
}
