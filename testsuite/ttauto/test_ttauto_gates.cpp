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

// The gate test inside the search (issue #2): on the main subautomaton of
// the n=6 stratum 3,3(2), the class with dilatation 2.01536 (the known
// spurious pseudo-Anosov) is reported with the gate test off, and with the
// gate test on it moves to the rejected list while every other class is
// unchanged.

#include <cmath>
#include <iostream>
#include <list>
#include <set>
#include <sstream>
#include <vector>
#include "check.hpp"
#include "traintracks/build.hpp"
#include "traintracks/traintrack.hpp"
#include "ttauto/pAclass.hpp"
#include "ttauto/ttauto.hpp"
#include "ttauto/ttfoldgraph.hpp"

using traintracks::traintrack;
typedef ttauto::ttfoldgraph<traintrack> ttgraph;
typedef ttauto::ttauto<traintrack> ttsearch;

// Polynomials of a class list as strings, for set comparison.
static std::set<std::string> polys(const ttsearch::pAlist& l)
{
  std::set<std::string> out;
  for (auto it = l.begin(); it != l.end(); ++it)
    {
      std::ostringstream os; os << it->first; out.insert(os.str());
    }
  return out;
}

static bool has_dilatation(const ttsearch::pAlist& l, const double lam)
{
  for (auto it = l.begin(); it != l.end(); ++it)
    if (std::fabs(it->second.dilatation() - lam) < 1e-6) return true;
  return false;
}

int main()
{
  jlt::vector<traintrack> ttv = traintracks::build_traintrack_list(6);
  ttgraph full(ttv[4]);
  std::list<ttgraph> sgs = ttauto::subgraphs(full);
  ttauto::prune_multihumps(sgs);
  const ttgraph& ttg = *sgs.begin();
  CHECK(ttg.vertices() == 90);

  const double lambda_bad = 2.01535718128;   // from pA_n=6_5_1_inv.m

  // Gate test off: the old behaviour.
  ttsearch off(ttg);
  off.check_gates(false).max_pathlength(7);
  off.search();
  CHECK(has_dilatation(off.pA_list(),lambda_bad));
  CHECK(off.rejected_pA_list().empty());
  CHECK(off.gate_rejected() == 0);

  // Gate test on (the default).
  ttsearch on(ttg);
  on.max_pathlength(7);
  on.search();
  CHECK(!has_dilatation(on.pA_list(),lambda_bad));
  CHECK(has_dilatation(on.rejected_pA_list(),lambda_bad));
  CHECK(on.gate_rejected() > 0);
  // The same candidates reach the gate test whether or not it is applied.
  CHECK(on.gate_candidates() == off.gate_candidates());
  CHECK(on.gate_candidates() > on.gate_rejected());
  CHECK(on.gate_rejection_rate() > 0.0 && on.gate_rejection_rate() < 1.0);

  // Every other class is unchanged, and nothing else was rejected.
  std::set<std::string> poff = polys(off.pA_list()), pon = polys(on.pA_list());
  std::set<std::string> prej = polys(on.rejected_pA_list());
  CHECK(prej.size() == 1);
  CHECK(poff.size() == pon.size() + 1);
  for (const std::string& q : pon) CHECK(poff.count(q) == 1);
  for (const std::string& q : prej) CHECK(poff.count(q) == 1 && pon.count(q) == 0);

  std::cout << "\ntest_ttauto_gates: classes with gates off " << poff.size()
            << ", with gates on " << pon.size() << ", rejected " << prej.size()
            << " (" << on.gate_rejected() << " of " << on.gate_candidates()
            << " candidate paths, " << 100*on.gate_rejection_rate() << "%)\n";
  std::cout << "Rejected classes:\n";
  on.print_rejected_pA_list(std::cout);
  std::cout << "test_ttauto_gates: OK\n";
  return 0;
}
