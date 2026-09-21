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

// Braids read off closed folding paths (ttauto/path_braid.hpp).
//
// Every braid is checked against the path it came from by a route with
// nothing in common with the extraction: the growth of its action on
// Dynnikov coordinates must equal the Perron root of the path's transition
// matrix.  That equality only holds when the path really is pseudo-Anosov,
// so the sample is the closed paths with a primitive matrix and connected
// gates.
//
// folding_path_braid reports the outcome through its `verified` flag, so
// that a braid which does not check out is never mistaken for one that
// does.  Nothing fails in the range covered here.

#include <cmath>
#include <iostream>
#include <vector>
#include "check.hpp"
#include "traintracks/braid.hpp"
#include "traintracks/build.hpp"
#include "traintracks/traintrack.hpp"
#include "ttauto/folding_path.hpp"
#include "ttauto/path_braid.hpp"
#include "ttauto/ttfoldgraph.hpp"

using traintracks::braidword;
using traintracks::traintrack;
typedef ttauto::ttfoldgraph<traintrack> ttgraph;
typedef ttauto::folding_path<traintrack> path;

static bool primitive(const jlt::mathmatrix<int>& A)
{
  const int m = A.dim();
  jlt::mathmatrix<int> B(A);
  for (int k = 0; k < m*m+1; ++k)
    {
      bool allpos = true;
      for (int i = 0; i < m && allpos; ++i)
        for (int j = 0; j < m && allpos; ++j) if (B(i,j) <= 0) allpos = false;
      if (allpos) return true;
      B = B*A;
      for (int i = 0; i < m; ++i)
        for (int j = 0; j < m; ++j) if (B(i,j) > 1000000) B(i,j) = 1000000;
    }
  return false;
}

// Every closed path of length <= maxlen from every vertex.
static void walk(const ttgraph& ttg, const int v0, const int maxlen,
                 std::vector<int>& fp, long& nver, long& nunver)
{
  if (!fp.empty())
    {
      path p(ttg,v0);
      for (std::size_t i = 0; i < fp.size(); ++i) p.push_back(fp[i]);
      if (p.closed())
        {
          const jlt::mathmatrix<int> TM = p.transition_matrix();
          if (primitive(TM) && p.gates().connected
              && ttauto::detail::perron_root(TM) > 1.05)
            {
              bool verified = false;
              const braidword b = ttauto::folding_path_braid(p,&verified);
              CHECK(b.strings() == ttg.traintrack(v0).punctures());
              if (verified) ++nver; else ++nunver;
            }
        }
    }
  if ((int)fp.size() >= maxlen) return;

  int v = v0;
  for (std::size_t i = 0; i < fp.size(); ++i) v = ttg.target_vertex(v,fp[i]);
  for (int br = 0; br < ttg.foldings(v); ++br)
    {
      fp.push_back(br);
      walk(ttg,v0,maxlen,fp,nver,nunver);
      fp.pop_back();
    }
}

int main()
{
  // The braidword arithmetic itself.
  {
    const braidword d = braidword::delta(4,1);
    CHECK(d.length() == 3);
    CHECK(d.exponent_sum() == 3);
    // delta moves the puncture at position k to k-1, and 1 round to n.
    const std::vector<int> perm = d.permutation();
    CHECK(perm[0] == 4);
    for (int k = 2; k <= 4; ++k) CHECK(perm[k-1] == k-1);
    // beta_{i,j} = sigma_i ... sigma_{j-1} is the one-puncture case.
    const braidword b = braidword::block_swap(5,2,1,2,1);
    CHECK(b.length() == 2);
    CHECK(b.word()[0] == 2 && b.word()[1] == 3);
    CHECK(b.permutation()[1] == 4);
    // Growth, against values that do not depend on any of this.
    CHECK(std::fabs(braidword(3,{1,-2}).growth() - 2.6180339887) < 1e-8);
    CHECK(std::fabs(braidword(5,{1,2,3,4,1,2}).growth() - 1.7220838057) < 1e-8);
    // A braid that is not pseudo-Anosov grows polynomially, so the
    // iteration only approaches 1 from above.
    CHECK(braidword(4,{1,1,1}).growth() < 1.01);
    CHECK(std::fabs(braidword(3,{1,2}).growth() - 1) < 1e-8);
    const braidword w(5,{1,2,-3,4});
    CHECK(std::fabs((w*w.inverse()).reduce().length()) == 0);
    CHECK(w.exponent_sum() == 2);
  }

  // Three punctures: the one vertex, the two folds, and the classical
  // minimum.  The whole pipeline in its smallest instance.
  {
    jlt::vector<traintrack> ttv = traintracks::build_traintrack_list(3);
    ttgraph ttg(ttv[0]);
    CHECK(ttg.vertices() == 1);
    path p(ttg,0);
    p.push_back(0); p.push_back(1);
    bool verified = false;
    const braidword b = ttauto::folding_path_braid(p,&verified);
    CHECK(verified);
    CHECK(b.word() == std::vector<int>({-1,2}));
    CHECK(std::fabs(b.growth() - 2.6180339887) < 1e-8);
  }

  long nver = 0, nunver = 0;
  for (int n = 3; n <= 5; ++n)
    {
      jlt::vector<traintrack> ttv = traintracks::build_traintrack_list(n);
      for (int s = 0; s < (int)ttv.size(); ++s)
        {
          ttgraph ttg(ttv[s]);
          for (int v = 0; v < ttg.vertices(); ++v)
            {
              std::vector<int> fp;
              walk(ttg,v,4,fp,nver,nunver);
            }
        }
    }
  CHECK(nunver == 0);
  CHECK(nver > 300);

  // The canonical spurious pseudo-Anosov of issue #2: six punctures,
  // stratum 3 3 (2), the closed path 1,0,3,2,1,2 with dilatation 2.01536
  // that the gate test rejects.  CLAUDE.md gives its cycle as
  // {29,46,43,71,88,85} in the ninety-vertex subgraph, one-based; the same
  // path in the whole automaton starts at vertex 36, zero-based, and runs
  // through 57, 54, 86, 107, 104.  Taking it here avoids depending on how
  // subgraphs() and prune_multihumps order their output.  Its braid must have
  // that growth and fix exactly one puncture, since it is the
  // five-puncture minimum with an idle sixth.  The path runs through
  // cyclically symmetric vertices, so it also guards the rule that each
  // step continues from the automaton's copy of the track.
  {
    jlt::vector<traintrack> ttv = traintracks::build_traintrack_list(6);
    ttgraph ttg(ttv[4]);
    path p(ttg,36);
    const int branches[6] = {1,0,3,2,1,2};
    for (int i = 0; i < 6; ++i) p.push_back(branches[i]);
    CHECK(p.closed());
    CHECK(!p.gates().connected);
    CHECK(std::fabs(ttauto::detail::perron_root(p.transition_matrix())
                    - 2.01535718128) < 1e-8);
    const braidword b = ttauto::folding_path_braid(p);
    CHECK(b.strings() == 6);
    CHECK(b.word()
          == std::vector<int>({1,2,1,2,3,4,5,5,-4,-3,-5,-5,-4,-3,-2,-1}));
    CHECK(b.exponent_sum() == 0);
    CHECK(std::fabs(b.growth() - 2.01535718128) < 1e-8);
    // A five-puncture pseudo-Anosov with an idle sixth, so the permutation
    // is a five-cycle and one fixed point, as it is for the braid worked
    // out by hand in devel/iss002, which differs from this one by
    // conjugation.
    const std::vector<int> perm = b.permutation();
    int fixed = 0;
    for (int i = 0; i < 6; ++i) if (perm[i] == i+1) ++fixed;
    CHECK(fixed == 1);
  }

  // The gate-rejected class at four punctures: the closed path 1,1,2,2 from
  // vertex 0 of the stratum 1.1.1.1.(2), with characteristic polynomial
  // (x-1)(x^2-3x+1).  Its braid is the three-puncture minimum with the
  // fourth puncture left alone, so the permutation fixes one puncture and
  // the growth is phi^2 even though the class is not pseudo-Anosov.
  {
    jlt::vector<traintrack> ttv = traintracks::build_traintrack_list(4);
    ttgraph ttg(ttv[0]);
    path p(ttg,0);
    p.push_back(1); p.push_back(1); p.push_back(2); p.push_back(2);
    CHECK(p.closed());
    CHECK(!p.gates().connected);
    const braidword b = ttauto::folding_path_braid(p);
    CHECK(b.strings() == 4);
    CHECK(b.word() == std::vector<int>({1,1,2,3,3,-2,-3,-3,-2,-1}));
    CHECK(b.exponent_sum() == 0);
    CHECK(std::fabs(b.growth() - 2.6180339887) < 1e-8);
    const std::vector<int> perm = b.permutation();
    int fixed = 0;
    for (int i = 0; i < 4; ++i) if (perm[i] == i+1) ++fixed;
    CHECK(fixed == 1);
  }

  std::cout << "\ntest_braid_extraction: " << nver << " braids verified against"
            << " the Perron root of their path, " << nunver << " not\n";
  return 0;
}
