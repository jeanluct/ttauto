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

// Norm-bounded search mode (check_norms + max_dilatation), the mode
// examples/ttauto_min_example uses to reproduce the literature minima for
// n = 3, 4, 5 (Song-Ko-Los 2002; Ham-Song 2007).  This is the one search
// mode the length-bounded tests do not exercise.  With the
// Bestvina-Handel gate test on, every minimum is reproduced; on the n=4
// first stratum the gate test rejects a reducible candidate with the
// same dilatation as the true minimum.

#include <cmath>
#include <iostream>
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

struct stratum_case
{
  int n, trk;
  double dilmax;          // as in examples/ttauto_min_example
  double lambda;          // literature minimum
  std::vector<int> coeffs;// c_0..c_edges of jlt's det(M - xI)
  int shortest;
  int rejected_classes;   // classes moved to rejected_pA_list()
};

static void run_case(const stratum_case& c)
{
  auto ttv = traintracks::build_traintrack_list(c.n);
  CHECK((int)ttv.size() > c.trk);
  ttauto::ttfoldgraph<traintrack> ttg(ttv[c.trk]);
  ttauto::ttauto<traintrack> tta(ttg);
  tta.max_dilatation(c.dilmax).badword_length(0).check_norms();
  tta.search();

  const auto& pal = tta.pA_list();
  CHECK_MSG(pal.size() == 1, "n=" << c.n << " trk=" << c.trk << ": " << pal.size() << " classes");
  const auto& cls = pal.begin()->second;
  CHECK_MSG(std::abs(cls.dilatation() - c.lambda) < 1e-4,
            "n=" << c.n << " trk=" << c.trk << ": " << cls.dilatation());
  CHECK((int)c.coeffs.size() == ttg.edges()+1);
  for (int k = 0; k <= ttg.edges(); ++k)
    CHECK_MSG(cls.polynomial()[k] == c.coeffs[k],
              "n=" << c.n << " trk=" << c.trk << ": " << cls.polynomial());
  CHECK(cls.shortest() == c.shortest);
  CHECK(tta.gate_candidates() > 0);
  CHECK_MSG((int)tta.rejected_pA_list().size() == c.rejected_classes,
            "n=" << c.n << " trk=" << c.trk << ": " << tta.rejected_pA_list().size());
  CHECK((tta.gate_rejected() > 0) == (c.rejected_classes > 0));

  for (const auto& pm : cls.paths())
    {
      CHECK(pm.second.is_primitive());
      CHECK(std::abs(oracles::power_iteration_radius(pm.second) - c.lambda) < 1e-4);
      CHECK(pm.first.gates().connected);
    }
  for (auto it = tta.rejected_pA_list().begin(); it != tta.rejected_pA_list().end(); ++it)
    for (const auto& pm : it->second.paths())
      CHECK(!pm.first.gates().connected);
}

int main()
{
  const std::vector<stratum_case> cases = {
    // n=3: the golden-mean-squared braid sigma_1 sigma_2^-1.
    {3, 0, 2.62, 2.61803, { 1,-3, 1}, 2, 0},
    // n=4, stratum 1 (Song-Ko-Los): also a reducible candidate with the
    // same dilatation, characteristic polynomial (x-1)(x^2-3x+1), which
    // only the gate test removes.
    {4, 0, 2.62, 2.61803, {-1, 2, 2,-1}, 3, 1},
    {4, 1, 2.30, 2.29663, { 1,-2, 0,-2, 1}, 3, 0},
    // n=5 (Ham-Song).
    {5, 0, 1.73, 1.72208, { 1,-1,-1,-1, 1}, 2, 0},
    {5, 1, 1.73, 1.72208, {-1, 0, 2, 2, 0,-1}, 3, 0},
    {5, 2, 2.16, 2.15372, {-1, 2, 0, 0, 2,-1}, 3, 0},
    {5, 3, 2.02, 2.01536, { 1,-1, 0,-4, 0,-1, 1}, 4, 0},
  };
  for (const stratum_case& c : cases) run_case(c);

  // The rejected n=4 class has the reducible polynomial (x-1)(x^2-3x+1).
  {
    auto ttv = traintracks::build_traintrack_list(4);
    ttauto::ttfoldgraph<traintrack> ttg(ttv[0]);
    ttauto::ttauto<traintrack> tta(ttg);
    tta.max_dilatation(2.62).badword_length(0).check_norms();
    tta.search();
    CHECK(tta.rejected_pA_list().size() == 1);
    const auto& rej = tta.rejected_pA_list().begin()->second;
    const std::vector<int> want = {-1, 4,-4, 1};   // x^3 - 4x^2 + 4x - 1
    for (int k = 0; k <= 3; ++k) CHECK(rej.polynomial()[k] == want[k]);
    CHECK(std::abs(rej.dilatation() - 2.61803) < 1e-4);
    CHECK(std::abs(rej.polynomial()(1.0)) == 0);   // root at 1: reducible

    // With the gate test off the reducible class is reported as a pA.
    ttauto::ttauto<traintrack> ttb(ttg);
    ttb.check_gates(false).max_dilatation(2.62).badword_length(0).check_norms();
    ttb.search();
    CHECK(ttb.pA_list().size() == 2);
    CHECK(ttb.rejected_pA_list().empty());
    CHECK(ttb.gate_rejected() == 0);
  }

  std::cout << "test_min_dilatations: OK (" << cases.size() << " strata)\n";
  return 0;
}
