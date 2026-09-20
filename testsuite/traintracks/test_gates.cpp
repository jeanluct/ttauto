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

// Bestvina-Handel gate test (issue #2): the known spurious pseudo-Anosov of
// the n=6, stratum 3,3(2) automaton is rejected at exactly the fixed
// 3-edge monogon, four genuine pseudo-Anosovs pass, the word-based and
// accumulator-based analyses agree, and the shape constraints of
// Bestvina-Handel Props. 3.3.3-3.3.4 hold on every connected primitive
// closed path of a small automaton.

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <vector>
#include <jlt/freeauto.hpp>
#include <jlt/mathmatrix.hpp>
#include <jlt/polynomial.hpp>
#include "check.hpp"
#include "traintracks/build.hpp"
#include "traintracks/gates.hpp"
#include "traintracks/map.hpp"
#include "traintracks/traintrack.hpp"
#include "ttauto/folding_path.hpp"
#include "ttauto/ttfoldgraph.hpp"

using jlt::freeauto;
using traintracks::fold_derivative;
using traintracks::gate_accumulator;
using traintracks::gate_analysis;
using traintracks::traintrack;
using traintracks::ttnumbering;
using ttauto::folding_path;
using ttauto::ttfoldgraph;

typedef ttfoldgraph<traintrack> ttgraph;

// Accumulator-based analysis of a closed path given by branch choices.
static gate_analysis analyse_branches(const ttgraph& ttg, const int v0,
                                      const std::vector<int>& br)
{
  gate_accumulator acc;
  int v = v0;
  for (int b : br)
    {
      const int vnext = ttg.target_vertex(v,b);
      acc.push_back(fold_derivative(ttg.traintrack_map(v,b),
                                    ttg.traintrack(vnext).numbering()));
      v = vnext;
    }
  CHECK(v == v0);
  return acc.analyse(ttg.traintrack(v0).numbering());
}

// Word-based analysis of the same path.
static gate_analysis analyse_words(const ttgraph& ttg, const int v0,
                                   const std::vector<int>& br)
{
  folding_path<traintrack> p(ttg,v0);
  for (int b : br) p.push_back(b);
  CHECK(p.closed());
  return traintracks::analyse_gates(ttg.traintrack(v0).numbering(),p.traintrack_map());
}

static void check_same(const gate_analysis& a, const gate_analysis& b)
{
  CHECK(a.connected == b.connected);
  CHECK(a.shapes_ok == b.shapes_ok);
  CHECK(a.vertices.size() == b.vertices.size());
  for (std::size_t i = 0; i < a.vertices.size(); ++i)
    {
      CHECK(a.vertices[i].name == b.vertices[i].name);
      CHECK(a.vertices[i].gates == b.vertices[i].gates);
      CHECK(a.vertices[i].joins == b.vertices[i].joins);
      CHECK(a.vertices[i].components == b.vertices[i].components);
    }
}

// Puncture corollary (note, Corollary 3.4): at a connected prong of a
// punctured multigon there is exactly one main gate, and the joins are
// exactly incoming-peripheral~main and main~outgoing-peripheral.
static void check_puncture_corollary(const gate_analysis& ga, const ttnumbering& N)
{
  for (const auto& v : ga.vertices)
    {
      if (v.prong < 0 || !v.connected) continue;
      CHECK(v.gates.size() == 3);
      // Distinct gate pairs joined: exactly two, each involving a
      // peripheral direction.
      std::map<int,int> gate_of;
      for (std::size_t g = 0; g < v.gates.size(); ++g)
        for (int d : v.gates[g]) gate_of[d] = g;
      std::set<std::pair<int,int> > pairs;
      for (const auto& t : v.joins)
        {
          const int a = gate_of[t.first], b = gate_of[t.second];
          pairs.insert(std::make_pair(std::min(a,b),std::max(a,b)));
          CHECK(N.is_side(t.first) || N.is_side(t.second));
        }
      CHECK(pairs.size() == 2);
      int nperiph = 0;
      for (const auto& g : v.gates)
        for (int d : g) if (N.is_side(d)) ++nperiph;
      CHECK(nperiph == 2);
    }
}

// Find every branch sequence realising a vertex cycle whose composed
// matrix matches M (row-major), as in pA_n=6_5_1_inv.m.
static std::vector<std::vector<int> > branches_for(const ttgraph& ttg,
                                                   const std::vector<int>& cyc,
                                                   const std::vector<int>& M)
{
  const int L = cyc.size() - 1;
  const int n = ttg.edges();
  std::vector<int> cur(L);
  std::vector<std::vector<int> > found;
  std::function<void(int)> rec = [&](int i) {
    if (i == L)
      {
        folding_path<traintrack> p(ttg,cyc[0]);
        for (int b : cur) p.push_back(b);
        const jlt::mathmatrix<int> A = p.transition_matrix();
        bool match = true;
        for (int r = 0; r < n; ++r)
          for (int c = 0; c < n; ++c)
            if (A(r,c) != M[r*n+c]) match = false;
        if (match) found.push_back(cur);
        return;
      }
    for (int b = 0; b < ttg.foldings(cyc[i]); ++b)
      if (ttg.target_vertex(cyc[i],b) == cyc[i+1]) { cur[i] = b; rec(i+1); }
  };
  rec(0);
  CHECK_MSG(!found.empty(), "no branch sequence realises the cycle with this matrix");
  return found;
}

// True if some unpunctured multigon has more gates than prongs, i.e. the
// derivative splits a prong into several gates.
static bool has_refined_unpunctured(const gate_analysis& ga, const ttnumbering& N)
{
  for (const auto& v : ga.vertices)
    if (v.prong < 0 && (int)v.gates.size() > N.prong[N.prong_number[v.multigon][0]].nprongs)
      return true;
  return false;
}

int main()
{
  // The 90-vertex main subautomaton of stratum 3,3(2) for n=6.
  jlt::vector<traintrack> ttv = traintracks::build_traintrack_list(6);
  ttgraph full(ttv[4]);
  std::list<ttgraph> sgs = ttauto::subgraphs(full);
  ttauto::prune_multihumps(sgs);
  const ttgraph& ttg = *sgs.begin();
  CHECK(ttg.vertices() == 90);
  const int n = ttg.edges();
  CHECK(n == 7);

  // ---- The bad path (0-based vertices, branches as in gates_check.cpp).
  {
    const std::vector<int> cyc = {28,45,42,70,87,84,28};
    const std::vector<int> M = {0,0,0,0,0,0,1, 1,1,1,0,0,0,0, 0,0,0,1,0,0,0,
                                1,2,0,0,0,0,0, 0,0,1,0,1,0,1, 0,0,1,0,2,0,0,
                                0,0,0,0,0,1,0};
    const std::vector<std::vector<int> > brs = branches_for(ttg,cyc,M);
    CHECK(brs.size() == 1);
    const std::vector<int>& br = brs[0];
    CHECK((br == std::vector<int>{1,0,3,2,1,2}));

    const gate_analysis ga = analyse_branches(ttg,cyc[0],br);
    const gate_analysis gw = analyse_words(ttg,cyc[0],br);
    check_same(ga,gw);
    std::cout << "Bad path:\n"; ga.print(std::cout);

    CHECK(!ga.connected);
    CHECK(!ga.shapes_ok);
    // Exactly one disconnected vertex: the prong of the punctured monogon
    // with three main edges, split into two main gates, two components.
    int nbad = 0;
    const ttnumbering N = ttg.traintrack(cyc[0]).numbering();
    for (const auto& v : ga.vertices)
      {
        if (v.connected) continue;
        ++nbad;
        CHECK(v.prong >= 0);
        CHECK(N.prong[v.prong].nprongs == 1 && N.prong[v.prong].punctured);
        CHECK(N.prong_letters[v.prong].size() == 3);
        CHECK(v.directions.size() == 5);
        CHECK(v.gates.size() == 4);
        CHECK(v.components == 2);
        // The map fixes this puncture: its loop is fixed by the derivative.
        CHECK(!v.shape_ok);
      }
    CHECK(nbad == 1);
  }

  // ---- Four control paths reported as pseudo-Anosov by the automaton.
  {
    struct control { std::vector<int> cyc; std::vector<int> M; };
    const std::vector<control> controls = {
      {{28,45,42,45,28},
       {2,2,1,0,0,0,0, 0,0,0,0,1,0,0, 0,0,0,0,0,1,0, 0,0,0,0,0,0,1,
        1,1,1,0,0,0,0, 0,0,0,1,0,0,0, 1,2,0,0,0,0,0}},
      {{2,37,39,19,22,20,12,2},
       {0,0,0,1,2,0,0, 1,1,0,0,1,0,0, 0,2,0,0,1,0,0, 0,0,0,0,0,1,0,
        0,0,0,0,0,0,1, 1,0,0,0,0,0,0, 0,0,1,0,2,0,0}},
      {{4,84,28,45,28,80,4},
       {0,0,0,0,1,0,0, 0,0,0,0,0,1,0, 0,0,1,1,0,0,1, 1,0,0,0,0,0,0,
        0,1,0,0,0,0,0, 0,0,4,1,0,0,2, 0,0,2,0,0,0,1}},
      {{21,42,45,28,45,42,21},
       {0,0,0,0,0,1,0, 0,0,4,2,0,0,1, 0,0,3,1,0,0,1, 1,0,0,0,0,0,0,
        0,1,0,0,0,0,0, 0,0,4,1,0,0,2, 0,0,0,0,1,0,0}},
    };
    bool refined_unpunctured = false;
    for (const control& c : controls)
      for (const std::vector<int>& br : branches_for(ttg,c.cyc,c.M))
        {
          const gate_analysis ga = analyse_branches(ttg,c.cyc[0],br);
          const gate_analysis gw = analyse_words(ttg,c.cyc[0],br);
          check_same(ga,gw);
          std::cout << "Control path starting at vertex " << c.cyc[0]+1 << ", branches";
          for (int b : br) std::cout << " " << b;
          std::cout << ":\n";
          ga.print(std::cout);
          CHECK(ga.connected);
          CHECK(ga.shapes_ok);
          const ttnumbering N = ttg.traintrack(c.cyc[0]).numbering();
          check_puncture_corollary(ga,N);
          if (has_refined_unpunctured(ga,N)) refined_unpunctured = true;
        }
    std::cout << (refined_unpunctured
                  ? "A control exhibits a refined prong at an unpunctured multigon.\n"
                  : "No control exhibits a refined prong at an unpunctured multigon.\n");
  }

  // ---- Broad sweep on a small automaton: every closed path of length <= 5
  // analyses without a fail-fast, and connected + primitive implies the
  // Bestvina-Handel shape pattern.
  {
    jlt::vector<traintrack> ttv5 = traintracks::build_traintrack_list(5);
    ttgraph g5(ttv5[0]);
    int nclosed = 0, nconnected = 0, naccepted = 0, nrefined = 0;
    std::vector<int> br;
    std::function<void(int,int,int)> rec = [&](int v0, int v, int depth) {
      if (depth > 0 && v == v0)
        {
          ++nclosed;
          const gate_analysis ga = analyse_branches(g5,v0,br);
          folding_path<traintrack> p(g5,v0);
          for (int b : br) p.push_back(b);
          const bool primitive = p.transition_matrix().is_primitive();
          if (ga.connected) ++nconnected;
          if (primitive && ga.connected) CHECK(ga.shapes_ok);
          if (primitive && ga.connected) check_puncture_corollary(ga,g5.traintrack(v0).numbering());
          if (primitive && ga.connected)
            {
              // Risk 1 of the plan: does the D-refinement of prongs at an
              // unpunctured multigon ever bite on an accepted path?
              ++naccepted;
              if (has_refined_unpunctured(ga,g5.traintrack(v0).numbering())) ++nrefined;
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
    CHECK(nclosed > 0);
    CHECK(naccepted > 0);
    std::cout << "n=5 sweep: " << nclosed << " closed paths, " << nconnected
              << " with connected gates, " << naccepted
              << " primitive and connected, of which " << nrefined
              << " have a refined prong at an unpunctured multigon\n";
  }

  // ---- The minimal example (n=4, stratum (2)): the 3-puncture golden-mean
  // pseudo-Anosov with an idle fourth puncture.  Nothing shorter fails for
  // n <= 4, and nothing fails on n=3 or on the n=4 stratum 3(1).
  {
    jlt::vector<traintrack> ttv4 = traintracks::build_traintrack_list(4);
    ttgraph g4(ttv4[0]);
    CHECK(g4.vertices() == 4);
    CHECK(g4.edges() == 3);

    const std::vector<int> br = {1,1,2,2};
    folding_path<traintrack> p(g4,0);
    for (int b : br) p.push_back(b);
    CHECK(p.closed());
    {
      const std::vector<int> vp = {0,2,0,3,0};
      CHECK(p.vertices().size() == vp.size());
      for (std::size_t i = 0; i < vp.size(); ++i) CHECK(p.vertices()[i] == vp[i]);
    }
    const jlt::mathmatrix<int> TM = p.transition_matrix();
    const std::vector<int> want = {0,0,1, 1,2,0, 0,1,2};
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) CHECK(TM(i,j) == want[3*i+j]);
    CHECK(TM.is_primitive());
    // Characteristic polynomial (x-1)(x^2-3x+1), spectral radius phi^2.
    const jlt::polynomial<int> cp = TM.charpoly();
    CHECK(cp(1) == 0);
    CHECK(std::abs(cp(2.6180339887)) < 1e-6);

    const gate_analysis ga = analyse_branches(g4,0,br);
    check_same(ga,analyse_words(g4,0,br));
    std::cout << "Minimal example (n=4 stratum (2), branches 1 1 2 2):\n";
    ga.print(std::cout);
    CHECK(!ga.connected);
    const ttnumbering N = g4.traintrack(0).numbering();
    int nbad = 0;
    for (const auto& v : ga.vertices)
      {
        if (v.connected) continue;
        ++nbad;
        CHECK(v.prong >= 0);
        CHECK(N.prong[v.prong].nprongs == 1 && N.prong[v.prong].punctured);
        CHECK(N.prong_letters[v.prong].size() == 3);
        CHECK(v.gates.size() == 4);
        CHECK(v.components == 2);
        // The idle puncture: its peripheral loop is fixed by the map.
        const int loop = N.side_letter(v.prong);
        const auto img = p.traintrack_map().get_action(loop);
        CHECK(img.size() == 1 && *img.begin() == loop);
      }
    CHECK(nbad == 1);

    // Exhaustive counts of primitive-but-disconnected closed paths by
    // length: none for n=3, none on n=4 stratum 3(1), and 0,0,0,4,10 on
    // n=4 stratum (2) for lengths 1..5.
    auto count_bad = [](const ttgraph& g, const int L) {
      int cnt = 0;
      std::vector<int> b;
      std::function<void(int,int,int)> rec = [&](int v0, int v, int d) {
        if (d == L)
          {
            if (v != v0) return;
            folding_path<traintrack> q(g,v0);
            for (int x : b) q.push_back(x);
            if (q.transition_matrix().is_primitive() && !q.gates().connected) ++cnt;
            return;
          }
        for (int x = 0; x < g.foldings(v); ++x)
          {
            b.push_back(x);
            rec(v0,g.target_vertex(v,x),d+1);
            b.pop_back();
          }
      };
      for (int v0 = 0; v0 < g.vertices(); ++v0) rec(v0,v0,0);
      return cnt;
    };
    jlt::vector<traintrack> ttv3 = traintracks::build_traintrack_list(3);
    ttgraph g3(ttv3[0]);
    ttgraph g4b(ttv4[1]);
    for (int L = 1; L <= 5; ++L)
      {
        CHECK_MSG(count_bad(g3,L) == 0, "n=3 length " << L);
        CHECK_MSG(count_bad(g4b,L) == 0, "n=4 stratum 3(1) length " << L);
      }
    const std::vector<int> expect4 = {0,0,0,4,10};
    for (int L = 1; L <= 5; ++L)
      CHECK_MSG(count_bad(g4,L) == expect4[L-1],
                "n=4 stratum (2) length " << L << ": " << count_bad(g4,L));
  }

  std::cout << "test_gates: OK\n";
  return 0;
}
