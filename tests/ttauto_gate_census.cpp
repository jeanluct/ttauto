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

// Gate-test census (issue #2): for every stratum with n punctures, search
// the main subautomaton as examples/ttauto_scan_strata.sh does (same
// per-stratum path lengths) with the Bestvina-Handel gate test on, and
// report every class rejected by the gate test alongside the accepted
// classes.  Also counts, on the n=5 first stratum, the closed paths of
// length <= 5 whose matrix is primitive but whose gates are disconnected.
//
// Usage: ttauto_gate_census [nmin [nmax]]   (defaults 3 6)

#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include "traintracks/build.hpp"
#include "traintracks/gates.hpp"
#include "traintracks/traintrack.hpp"
#include "ttauto/folding_path.hpp"
#include "ttauto/pAclass.hpp"
#include "ttauto/ttauto.hpp"
#include "ttauto/ttfoldgraph.hpp"

using traintracks::traintrack;
typedef ttauto::ttfoldgraph<traintrack> ttgraph;
typedef ttauto::ttauto<traintrack> ttsearch;

// Path lengths of examples/ttauto_scan_strata.sh (key "n:stratum", 1-based).
static int max_len_for(const int n, const int stratum)
{
  static const std::map<std::string,int> table = {
    {"6:3",5}, {"6:5",8}, {"6:7",8},
    {"7:6",5}, {"7:7",7}, {"7:11",8}, {"7:12",8},
  };
  std::ostringstream key; key << n << ":" << stratum;
  auto it = table.find(key.str());
  return (it == table.end() ? 4 : it->second);
}

static std::string polystr(const jlt::polynomial<int>& p)
{
  std::ostringstream os; os << p; return os.str();
}

int main(int argc, char** argv)
{
  const int nmin = (argc > 1 ? std::atoi(argv[1]) : 3);
  const int nmax = (argc > 2 ? std::atoi(argv[2]) : 6);

  std::cout.precision(6);
  std::cout.setf(std::ios::showpoint);

  int total_rejected_classes = 0, total_rejected_paths = 0, total_partial = 0;
  std::vector<std::string> summary;

  for (int n = nmin; n <= nmax; ++n)
    {
      jlt::vector<traintrack> ttv = traintracks::build_traintrack_list(n);
      for (std::size_t trk = 0; trk < ttv.size(); ++trk)
        {
          const int stratum = trk + 1;
          const int len = max_len_for(n,stratum);
          ttgraph full(ttv[trk]);
          std::list<ttgraph> sgs = ttauto::subgraphs(full);
          const ttgraph& ttg = *sgs.begin();

          // Silence the search's own chatter.
          std::ostringstream sink;
          std::streambuf* old = std::cout.rdbuf(sink.rdbuf());
          ttsearch tta(ttg);
          tta.max_pathlength(len);
          tta.search();
          std::cout.rdbuf(old);

          const ttsearch::pAlist& acc = tta.pA_list();
          const ttsearch::pAlist& rej = tta.rejected_pA_list();

          std::cout << "n=" << n << " stratum " << stratum << " (";
          ttv[trk].print_singularity_data(std::cout);
          std::cout << "), main subautomaton " << ttg.vertices()
                    << " vertices, path length <= " << len << ": "
                    << acc.size() << " classes accepted, "
                    << rej.size() << " classes rejected ("
                    << tta.gate_rejected() << " paths)\n";

          for (auto it = rej.begin(); it != rej.end(); ++it)
            {
              const bool partial = (acc.find(it->first) != acc.end());
              std::cout << "   REJECTED  lambda = " << it->second.dilatation()
                        << "  " << polystr(it->first)
                        << "  paths kept " << it->second.number_of_paths()
                        << (partial ? "  [class also has accepted paths]" : "")
                        << "\n";
              it->second.print_paths(std::cout,3);
              std::cout << "\n";
              ++total_rejected_classes;
              if (partial) ++total_partial;
              std::ostringstream os;
              os << "n=" << n << " stratum " << stratum << " lambda=" << it->second.dilatation()
                 << (partial ? " (partial)" : "");
              summary.push_back(os.str());
            }
          total_rejected_paths += tta.gate_rejected();
        }
    }

  std::cout << "\n==== Census summary (n=" << nmin << ".." << nmax << ")\n"
            << "rejected classes: " << total_rejected_classes
            << ", of which with some accepted representative: " << total_partial
            << "; rejected paths: " << total_rejected_paths << "\n";
  for (const std::string& s : summary) std::cout << "  " << s << "\n";

  // n=5 first stratum: primitive but disconnected closed paths of length <= 5.
  {
    jlt::vector<traintrack> ttv5 = traintracks::build_traintrack_list(5);
    ttgraph g5(ttv5[0]);
    int nclosed = 0, nprim = 0, nprim_disc = 0;
    std::vector<int> br;
    std::function<void(int,int,int)> rec = [&](int v0, int v, int depth) {
      if (depth > 0 && v == v0)
        {
          ++nclosed;
          ttauto::folding_path<traintrack> p(g5,v0);
          for (int b : br) p.push_back(b);
          if (p.transition_matrix().is_primitive())
            {
              ++nprim;
              if (!p.gates().connected) ++nprim_disc;
            }
        }
      if (depth == 5) return;
      for (int b = 0; b < g5.foldings(v); ++b)
        {
          br.push_back(b);
          rec(v0,g5.target_vertex(v,b),depth+1);
          br.pop_back();
        }
    };
    for (int v0 = 0; v0 < g5.vertices(); ++v0) rec(v0,v0,0);
    std::cout << "\nn=5 first stratum, all closed paths of length <= 5 from every vertex: "
              << nclosed << " closed, " << nprim << " with primitive matrix, "
              << nprim_disc << " primitive but gates disconnected\n";
  }
  return 0;
}
