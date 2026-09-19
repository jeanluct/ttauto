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

// Length-bounded search on the n=5 first stratum: exact class list at
// path length 8, gate counters, and independent oracles for the two jlt
// routines the acceptance depends on (characteristic polynomial from
// principal minors; dilatation by power iteration on the stored matrix).

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>
#include "traintracks/build.hpp"
#include "traintracks/traintrack.hpp"
#include "ttauto/pAclass.hpp"
#include "ttauto/ttauto.hpp"
#include "ttauto/ttfoldgraph.hpp"
#include "check.hpp"
#include "oracles.hpp"

using traintracks::traintrack;
typedef jlt::polynomial<int> Poly;

// jlt's charpoly() is det(M - xI); the oracle computes det(xI - M).
static bool same_charpoly(const Poly& jltp, const Poly& oracle, const int n)
{
  const int sign = (n % 2 == 0 ? 1 : -1);
  for (int k = 0; k <= n; ++k)
    if (jltp[k] != sign*oracle[k]) return false;
  return true;
}

int main()
{
  const int n = 5;
  const int trk = 0;
  auto ttv = traintracks::build_traintrack_list(n);
  CHECK((int)ttv.size() > trk);

  ttauto::ttfoldgraph<traintrack> ttg(ttv[trk]);
  CHECK(ttg.vertices() == 11);
  CHECK(ttg.edges() == 4);

  ttauto::ttauto<traintrack> tta(ttg);
  tta.max_pathlength(8);
  tta.max_dilatation(3.0);
  tta.badword_length(0);
  tta.search();

  // Exact expectations (verified against examples/ttauto_scan_strata.md:
  // stratum 1 of n=5 has minimum 1.72208 at length 2).
  const auto& pal = tta.pA_list();
  CHECK(pal.size() == 8);
  CHECK(tta.gate_candidates() == 3555);
  CHECK(tta.gate_rejected() == 0);
  CHECK(tta.rejected_pA_list().empty());
  CHECK(tta.gate_rejection_rate() == 0.0);

  struct expected { double lambda; std::vector<int> coeffs; int shortest; int longest; };
  // coeffs are c_0..c_4 of jlt's det(M - xI).
  const std::vector<expected> exp = {
    {1.72208, { 1,-1,-1,-1, 1}, 2, 3},
    {2.08102, { 1,-1,-2,-1, 1}, 3, 5},
    {2.36921, { 1,-1,-3,-1, 1}, 4, 7},
    {2.61803, { 1,-1,-4,-1, 1}, 5, 8},
    {2.84054, { 1,-1,-5,-1, 1}, 6, 8},
    {2.29663, { 1,-2, 0,-2, 1}, 3, 5},
    {2.89005, { 1,-2,-2,-2, 1}, 4, 7},
    {2.96557, { 1,-3, 1,-3, 1}, 4, 6},
  };
  std::vector<bool> seen(exp.size(),false);
  double minlambda = 1e10;
  int npaths = 0;
  for (auto it = pal.begin(); it != pal.end(); ++it)
    {
      const auto& cls = it->second;
      CHECK(cls.number_of_paths() > 0);
      CHECK(cls.shortest() > 0);
      CHECK(cls.longest() >= cls.shortest());
      CHECK(cls.dilatation() > 1.0);
      CHECK(it->first == cls.polynomial());
      minlambda = std::min(minlambda,cls.dilatation());

      bool matched = false;
      for (size_t k = 0; k < exp.size(); ++k)
        {
          if (std::abs(cls.dilatation() - exp[k].lambda) > 1e-4) continue;
          bool same = true;
          for (int c = 0; c <= n-1; ++c) if (cls.polynomial()[c] != exp[k].coeffs[c]) same = false;
          if (!same) continue;
          CHECK(!seen[k]);
          seen[k] = true;
          matched = true;
          CHECK(cls.shortest() == exp[k].shortest);
          CHECK(cls.longest() == exp[k].longest);
        }
      CHECK_MSG(matched, "unexpected class with dilatation " << cls.dilatation());

      // Oracles on every stored representative: the class polynomial is
      // the characteristic polynomial of the stored matrix, the matrix is
      // primitive, and its spectral radius is the class dilatation.
      for (const auto& pm : cls.paths())
        {
          ++npaths;
          const auto& p = pm.first;
          const auto& TM = pm.second;
          CHECK(p.closed());
          CHECK(p.transition_matrix() == TM);
          CHECK(TM.is_primitive());
          CHECK(same_charpoly(cls.polynomial(),oracles::charpoly_by_minors(TM),ttg.edges()));
          CHECK(TM.charpoly() == cls.polynomial());
          const double rho = oracles::power_iteration_radius(TM);
          CHECK_MSG(std::abs(rho - cls.dilatation()) < 1e-8,
                    "power iteration " << rho << " vs findroot " << cls.dilatation());
          CHECK(std::abs(cls.polynomial()(cls.dilatation())) < 1e-6);
          CHECK(p.gates().connected);
        }
    }
  for (bool s : seen) CHECK(s);
  CHECK(std::abs(minlambda - 1.72208) < 1e-4);
  CHECK(npaths >= (int)pal.size());

  // Reporting helpers stay callable and non-empty.
  std::ostringstream out;
  tta.print_pA_list(out);
  CHECK(!out.str().empty());
  out.str("");
  tta.print_pA_list_MathematicaForm(out);
  CHECK(out.str().find("polynomial ->") != std::string::npos);
  out.str("");
  tta.print_rejected_pA_list(out);

  std::cout << "test_ttauto_search: OK (" << pal.size() << " classes, "
            << npaths << " representatives checked against oracles)\n";
  return 0;
}
