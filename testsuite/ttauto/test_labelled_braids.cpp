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

// The braid extraction, cross-checked against the labelled automaton.
//
// pure_braid() gives every puncture a distinct label, so a folding path
// closes up only if the homeomorphism it defines returns every puncture
// to itself, that is, only if its braid is pure.  The automaton decides
// that by comparing labelled codings, which is pure combinatorics; the
// braid's permutation comes from outer_embedding and fold_block_swap.
// The two share no code, so their agreement is worth asserting.
//
// Part B is the control.  Without it, "every labelled path gives a pure
// braid" and "no labelled vertex is cyclically symmetric" could both hold
// because the enumeration found nothing.
//
// Part C guards 653134f, the hardest bug of issue #4: at a cyclically
// symmetric track the automaton identifies the final track with the
// initial one only up to that symmetry, so folding_path_braids returns
// several candidates differing by a root of the full twist, and the right
// one has to be picked.  Note this needs an *unlabelled* automaton:
// labelling destroys the cyclic symmetry, so the labelled parts above
// cannot exercise that path at all.  test_braid_extraction.cpp guards it
// in aggregate, by asserting no pseudo-Anosov path goes unverified; this
// checks the structure directly.

#include <cmath>
#include <iostream>
#include <vector>
#include <jlt/mathmatrix.hpp>
#include <jlt/vector.hpp>
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

// is_cyclically_symmetric() is not const, and returns the order less one,
// so 0 means no symmetry.
static int cyclic_order_minus_one(const ttgraph& ttg, const int v)
{
  traintrack tt(ttg.traintrack(v));
  return tt.is_cyclically_symmetric();
}

static bool is_identity(const std::vector<int>& perm)
{
  for (int i = 0; i < (int)perm.size(); ++i) if (perm[i] != i+1) return false;
  return true;
}

// Every closed path of length <= maxlen from v0.  Counts the pure and
// the non-pure braids separately; in a labelled automaton the second
// count must stay zero.
static void walk(const ttgraph& ttg, const int v0, const int maxlen,
                 const bool labelled, std::vector<int>& fp,
                 long& npure, long& nmixed)
{
  if (!fp.empty())
    {
      path p(ttg,v0);
      for (std::size_t i = 0; i < fp.size(); ++i) p.push_back(fp[i]);
      if (p.closed())
        {
          if (labelled)
            {
              // No cyclic symmetry, so there is nothing to choose
              // between and the braid is unambiguous even when the path
              // is not pseudo-Anosov and so cannot be verified.
              CHECK(ttauto::folding_path_braids(p).size() == 1);
            }
          const braidword b = ttauto::folding_path_braid(p);
          CHECK(b.strings() == ttg.traintrack(v0).punctures());
          if (is_identity(b.permutation())) ++npure; else ++nmixed;
        }
    }
  if ((int)fp.size() >= maxlen) return;

  int v = v0;
  for (std::size_t i = 0; i < fp.size(); ++i) v = ttg.target_vertex(v,fp[i]);
  for (int br = 0; br < ttg.foldings(v); ++br)
    {
      fp.push_back(br);
      walk(ttg,v0,maxlen,labelled,fp,npure,nmixed);
      fp.pop_back();
    }
}

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

int main()
{
  using std::cout;
  using std::endl;

  //
  // Part A: in a labelled automaton every closed path gives a pure braid.
  //
  long npure = 0, nmixed = 0, nvertices = 0;
  for (int n = 3; n <= 4; ++n)
    {
      jlt::vector<traintrack> ttv = traintracks::build_traintrack_list(n);
      for (int s = 0; s < (int)ttv.size(); ++s)
        {
          traintrack tt(ttv[s]);
          tt.pure_braid();
          const ttgraph ttg(tt);

          // Distinct labels on every puncture leave no cyclic symmetry
          // for a rotation to preserve.
          for (int v = 0; v < ttg.vertices(); ++v)
            {
              CHECK_MSG(cyclic_order_minus_one(ttg,v) == 0, v);
              ++nvertices;
            }

          const int maxlen = (n == 3 ? 6 : 4);
          for (int v = 0; v < ttg.vertices(); ++v)
            {
              std::vector<int> fp;
              walk(ttg,v,maxlen,true,fp,npure,nmixed);
            }
        }
    }
  cout << "labelled: " << nvertices << " vertices, none cyclically symmetric"
       << endl;
  cout << "          " << npure << " closed paths, all giving pure braids"
       << endl;
  CHECK(nmixed == 0);
  CHECK(npure > 100);

  //
  // Part B: the control.  Unlabelled, the same strata do produce braids
  // that permute the punctures, and do have cyclically symmetric
  // vertices.  Both are what Part A depends on being absent.
  //
  long upure = 0, umixed = 0, usym = 0;
  for (int n = 3; n <= 4; ++n)
    {
      jlt::vector<traintrack> ttv = traintracks::build_traintrack_list(n);
      for (int s = 0; s < (int)ttv.size(); ++s)
        {
          const ttgraph ttg(ttv[s]);
          for (int v = 0; v < ttg.vertices(); ++v)
            if (cyclic_order_minus_one(ttg,v) > 0) ++usym;

          const int maxlen = (n == 3 ? 6 : 4);
          for (int v = 0; v < ttg.vertices(); ++v)
            {
              std::vector<int> fp;
              walk(ttg,v,maxlen,false,fp,upure,umixed);
            }
        }
    }
  cout << "unlabelled: " << upure << " pure and " << umixed
       << " non-pure closed paths, " << usym
       << " cyclically symmetric vertices" << endl;
  CHECK(umixed > 0);
  CHECK(usym > 0);

  //
  // Part C: at a cyclically symmetric vertex the braid really is
  // ambiguous, and the right candidate is the one that gets picked.
  //
  {
    jlt::vector<traintrack> ttv = traintracks::build_traintrack_list(4);
    const ttgraph ttg(ttv[0]);

    int vsym = -1, order = 0;
    for (int v = 0; v < ttg.vertices() && vsym < 0; ++v)
      if (cyclic_order_minus_one(ttg,v) > 0)
        { vsym = v; order = cyclic_order_minus_one(ttg,v) + 1; }
    CHECK(vsym >= 0);
    CHECK(order > 1);

    // The shortest pseudo-Anosov closed path at that vertex.
    bool found = false;
    for (int len = 1; len <= 6 && !found; ++len)
      {
        path p(ttg,vsym,len);
        do
          {
            if (!p.closed()) continue;
            const jlt::mathmatrix<int> TM = p.transition_matrix();
            if (!primitive(TM)) continue;
            if (!p.gates().connected) continue;
            const double lambda = ttauto::detail::perron_root(TM);
            if (lambda < 1.05) continue;

            // One candidate per power of the symmetry: the automaton
            // only matched the final track to the initial one up to it.
            const std::vector<braidword> all = ttauto::folding_path_braids(p);
            CHECK_MSG((int)all.size() == order, all.size());

            // They are genuinely different, so the choice matters.
            bool alldistinct = true;
            for (int i = 1; i < (int)all.size(); ++i)
              if (all[i].word() == all[0].word()) alldistinct = false;
            CHECK(alldistinct);

            // And the one chosen is the one whose growth matches the
            // path's own dilatation.  This is what 653134f fixed.
            bool verified = false;
            const braidword b = ttauto::folding_path_braid(p,&verified);
            CHECK(verified);
            CHECK(std::fabs(b.growth() - lambda) < 1e-6*lambda);
            CHECK(b.strings() == 4);

            cout << "symmetric vertex " << vsym << " of order " << order
                 << ": " << all.size() << " candidates, dilatation "
                 << lambda << endl;
            found = true;
            break;
          }
        while (++p);
      }
    CHECK(found);
  }

  cout << "\ntest_labelled_braids: all checks passed" << endl;
  return 0;
}
