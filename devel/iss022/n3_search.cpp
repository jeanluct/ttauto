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

// Certified complete search of the n=3 automaton to trace T, that is to
// dilatation (T + sqrt(T^2-4))/2.  Prints one line per stored accepted
// path:  dilatation | (unused) | path | braid word | verified.
// The summary (classes by characteristic polynomial, and "cut by cap",
// which must be 0 for the list to be complete) goes to stderr.  See
// devel/iss022/n3_benchmark.md.
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include "traintracks/braid.hpp"
#include "traintracks/build.hpp"
#include "traintracks/traintrack.hpp"
#include "ttauto/folding_path.hpp"
#include "ttauto/path_braid.hpp"
#include "ttauto/ttauto.hpp"
#include "ttauto/ttfoldgraph.hpp"
using traintracks::traintrack;
int main(int argc, char** argv)
{
  const int T = std::atoi(argv[1]);
  const double lam = 0.5*(T + std::sqrt((double)T*T - 4.0));
  jlt::vector<traintrack> ttv = traintracks::build_traintrack_list(3);
  const ttauto::ttfoldgraph<traintrack> g(ttv[0]);
  ttauto::ttauto<traintrack> tta(g);
  tta.max_dilatation(lam*(1+1e-9)).check_norms();
  tta.max_paths_to_save(10000000);
  std::ostringstream sink;
  std::streambuf* sv = std::cout.rdbuf(sink.rdbuf());
  tta.search();
  std::cout.rdbuf(sv);
  std::cerr << "T=" << T << " lambda<=" << lam << " classes(by charpoly) " << tta.pA_list().size()
            << " rejected " << tta.rejected_pA_list().size()
            << " cut by cap " << tta.path_length_exceeded() << "\n";
  for (auto it = tta.pA_list().begin(); it != tta.pA_list().end(); ++it)
    for (auto pt = it->second.paths().begin(); pt != it->second.paths().end(); ++pt)
      {
        const ttauto::folding_path<traintrack>& p = pt->first;
        bool verified = false;
        const traintracks::braidword b = ttauto::folding_path_braid(p,&verified);
        std::cout.precision(15); std::cout << it->second.dilatation() << "|";
        std::cout << "|" << p << "|";
        for (int x : b.word()) std::cout << x << " ";
        std::cout << "|" << verified << "\n";
      }
}
