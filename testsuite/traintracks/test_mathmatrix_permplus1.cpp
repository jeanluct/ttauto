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

// mathmatrix_permplus1 against dense matrices: every one-fold transition
// matrix of the n=4 and n=5 (first stratum) automata round-trips through
// the sparse form, agrees with the unit-weight oracle, has a +1 entry
// exactly when the fold is legal, has the permutation order it claims, and
// multiplies dense matrices correctly on both sides.

#include <iostream>
#include <random>
#include <sstream>
#include <vector>
#include <jlt/mathmatrix.hpp>
#include <jlt/vector.hpp>
#include "traintracks/build.hpp"
#include "traintracks/map.hpp"
#include "traintracks/mathmatrix_permplus1.hpp"
#include "traintracks/traintrack.hpp"
#include "ttauto/ttfoldgraph.hpp"
#include "check.hpp"
#include "oracles.hpp"

using traintracks::mathmatrix_permplus1;
using traintracks::traintrack;
typedef jlt::mathmatrix<int> Mat;

static Mat dense_power(const Mat& P, const int k)
{
  Mat R(jlt::identity_matrix<int>(P.rows()));
  for (int i = 0; i < k; ++i) R = R*P;
  return R;
}

static void check_case(const mathmatrix_permplus1& pm, std::mt19937& rng)
{
  const Mat M = pm.full();
  const int n = M.rows();
  CHECK((int)pm.dim() == n && (int)pm.rows() == n && (int)pm.columns() == n);

  // Round trip through the dense constructor.
  mathmatrix_permplus1 rt(M);
  CHECK(rt.full() == M);
  CHECK(rt.row_perm() == pm.row_perm());
  CHECK(rt.column_perm() == pm.column_perm());
  CHECK(rt.plus1_row() == pm.plus1_row());
  CHECK(rt.plus1_col() == pm.plus1_col());
  CHECK(rt.is_perm() == pm.is_perm());

  // Row and column permutations are mutually inverse.
  for (int i = 0; i < n; ++i)
    {
      CHECK(pm.column_perm()[pm.row_perm()[i]] == i);
      CHECK(M(i,pm.row_perm()[i]) == 1);
    }

  // Column sums: all 1 for a permutation, one column summing to 2 otherwise.
  int ncol2 = 0;
  for (int j = 0; j < n; ++j)
    {
      int cs = 0;
      for (int i = 0; i < n; ++i) cs += M(i,j);
      if (cs == 2) ++ncol2; else CHECK(cs == 1);
    }
  CHECK(ncol2 == (pm.is_perm() ? 0 : 1));
  if (!pm.is_perm())
    {
      CHECK(M(pm.plus1_row(),pm.plus1_col()) == 1);
      CHECK(pm.plus1_row() != pm.row_perm()[pm.plus1_col()] ||
            pm.row_perm()[pm.plus1_row()] != pm.plus1_col());
    }
  CHECK(pm.is_identity() == (M == jlt::identity_matrix<int>(n)));

  // order(): the permutation part P satisfies P^order == I and no smaller
  // positive power does.
  Mat P(n,n);
  for (int i = 0; i < n; ++i) P(i,pm.row_perm()[i]) = 1;
  const int ord = pm.order();
  CHECK(ord >= 1 && ord <= n);
  CHECK(dense_power(P,ord) == jlt::identity_matrix<int>(n));
  for (int k = 1; k < ord; ++k)
    CHECK(dense_power(P,k) != jlt::identity_matrix<int>(n));

  // Products agree with dense multiplication on both sides.
  Mat A(n,n);
  std::uniform_int_distribution<int> pick(-5,5);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) A(i,j) = pick(rng);
  CHECK(pm * A == M * A);
  CHECK(A * pm == A * M);
  CHECK(pm * M == M * M);

  // Printers produce output.
  std::ostringstream os, om;
  os << pm;
  traintracks::printMathematicaForm(om,pm);
  CHECK(!os.str().empty() && !om.str().empty());
}

static int check_graph(const int n, const int trk, std::mt19937& rng)
{
  jlt::vector<traintrack> ttv = traintracks::build_traintrack_list(n);
  ttauto::ttfoldgraph<traintrack> ttg(ttv[trk]);
  int ncases = 0;
  for (int v = 0; v < ttg.vertices(); ++v)
    {
      const traintrack& tt = ttg.traintrack(v);
      // Every fold index, through the free function and the member.
      std::vector<Mat> fold_matrices;
      for (int f = 0; f < tt.foldings(); ++f)
        {
          const mathmatrix_permplus1 pm = traintracks::fold_transition_matrix(tt,f);
          CHECK(pm.full() == oracles::unit_weight_transition_matrix(tt,f));
          // (The member folds the track in place, so use a fresh copy.)
          traintrack ttm(tt);
          CHECK(ttm.fold_transition_matrix(f).full() == pm.full());
          traintrack ttf(tt);
          const bool legal = ttf.fold(f);
          CHECK(ttm == ttf);
          CHECK(pm.is_perm() == !legal);
          CHECK(pm.is_identity() == !legal);
          if (legal) fold_matrices.push_back(pm.full());
          check_case(pm,rng);
          ++ncases;
        }
      // Every stored branch matrix is one of the legal fold matrices and
      // abelianises the stored map.
      for (int b = 0; b < ttg.foldings(v); ++b)
        {
          const mathmatrix_permplus1& pm = ttg.transition_matrix(v,b);
          CHECK(!pm.is_perm());
          bool found = false;
          for (const Mat& F : fold_matrices) if (F == pm.full()) found = true;
          CHECK(found);
          CHECK(pm.full() == traintracks::transition_matrix_from_map(tt,ttg.traintrack_map(v,b)));
        }
    }
  return ncases;
}

int main()
{
  std::mt19937 rng(20260919);

  // Hand-built cases: a 6-cycle-free permutation and permutation+1.
  {
    const int n = 6;
    Mat Mperm(n,n);
    Mperm(0,4) = 1; Mperm(1,0) = 1; Mperm(2,1) = 1;
    Mperm(3,2) = 1; Mperm(4,3) = 1; Mperm(5,5) = 1;
    mathmatrix_permplus1 pperm(Mperm);
    CHECK(pperm.is_perm() && !pperm.is_identity());
    CHECK(pperm.order() == 5);
    check_case(pperm,rng);

    Mat Mplus1(Mperm);
    Mplus1(5,4) = 1;
    mathmatrix_permplus1 pplus1(Mplus1);
    CHECK(!pplus1.is_perm());
    CHECK(pplus1.plus1_row() == 5 && pplus1.plus1_col() == 4);
    check_case(pplus1,rng);

    mathmatrix_permplus1 pid(jlt::identity_matrix<int>(n));
    CHECK(pid.is_identity() && pid.order() == 1);
    check_case(pid,rng);

    // Construction from a row permutation and explicit +1 coordinates.
    jlt::vector<int> rp(n);
    for (int i = 0; i < n; ++i) rp[i] = (i+1) % n;
    mathmatrix_permplus1 pcyc(rp,2,0);
    CHECK(pcyc.order() == n);
    CHECK(pcyc.full()(2,0) == 1);
    check_case(pcyc,rng);
  }

  int ncases = 0;
  ncases += check_graph(4,0,rng);
  ncases += check_graph(4,1,rng);
  ncases += check_graph(5,0,rng);
  CHECK(ncases > 0);
  std::cout << "test_mathmatrix_permplus1: OK (" << ncases
            << " one-fold matrices checked)\n";
  return 0;
}
