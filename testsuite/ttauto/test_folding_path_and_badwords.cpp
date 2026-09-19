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

// folding_path algebra (composition, subpaths, cyclic equality, hashing)
// and the content of the badwords table: every entry of badwords(ttg,L)
// at (v,l) is the square h*h of a closed path h of length l+1 from v whose
// transition matrix has the same zero pattern as its square; the length-1
// layer is complete; and the table is deterministic.

#include <iostream>
#include <list>
#include <unordered_map>
#include <jlt/mathmatrix.hpp>
#include "traintracks/traintrack.hpp"
#include "ttauto/badwords.hpp"
#include "ttauto/folding_path.hpp"
#include "ttauto/ttfoldgraph.hpp"
#include "check.hpp"

using traintracks::traintrack;
using ttauto::folding_path;
typedef folding_path<traintrack> fpath;
typedef jlt::mathmatrix<int> Mat;

static bool same_foldings(const fpath& a, const fpath& b)
{
  return a.initial_vertex() == b.initial_vertex() &&
    a.foldings() == b.foldings() && a.vertices() == b.vertices();
}

int main()
{
  const int n = 6;
  traintrack tt(n,3);
  ttauto::ttfoldgraph<traintrack> ttg(tt);
  CHECK(ttg.vertices() > 0);

  // ---- Path algebra on a fixed path.
  fpath p1(ttg,1);
  CHECK(p1.zerolength() && p1.length() == 0);
  CHECK(p1.initial_vertex() == 1 && p1.final_vertex() == 1);
  CHECK(p1.number_of_foldings() == ttg.foldings(1));
  p1.push_back(0);
  p1.push_back(0);
  p1.push_back(1);
  p1.push_back(1);
  CHECK(p1.length() == 4 && !p1.zerolength());
  CHECK((int)p1.vertices().size() == 5);
  for (int i = 0; i < 4; ++i)
    CHECK(p1.vertices()[i+1] == ttg.target_vertex(p1.vertices()[i],p1.foldings()[i]));
  CHECK(p1.closed());

  // Subpaths from either end recompose to the whole path.
  const fpath head = p1.subpath(2), tail = p1.subpath(-2);
  CHECK(head.length() == 2 && tail.length() == 2);
  CHECK(head.initial_vertex() == p1.initial_vertex());
  CHECK(tail.final_vertex() == p1.final_vertex());
  CHECK(head.final_vertex() == tail.initial_vertex());
  CHECK(same_foldings(head*tail,p1));
  fpath acc(head);
  acc *= tail;
  CHECK(same_foldings(acc,p1));
  CHECK(p1.ending_equals(tail));
  CHECK(!p1.ending_equals(head) || same_foldings(head,tail));
  CHECK(p1.subpath(0).zerolength());
  CHECK(same_foldings(p1.subpath(10),p1));

  // Matrix and map of a composite are the products of the pieces.
  CHECK(p1.transition_matrix() == tail.transition_matrix()*head.transition_matrix());

  // Cyclic equality and hashing of closed paths.
  fpath p2(p1);
  p2.cycle_path(-1);
  CHECK(p1 == p2);
  CHECK(!same_foldings(p1,p2) || p1.foldings()[0] == p1.foldings()[3]);
  CHECK(p2.closed());
  fpath p3(p1);
  p3.cycle_path(4);
  CHECK(same_foldings(p3,p1));

  typedef fpath::hash path_hash;
  std::unordered_map<fpath,int,path_hash> pl;
  pl[p1] = 1;
  pl[p2] = 2;
  CHECK(pl.size() == 1);
  CHECK(pl[p1] == 2);

  // An open path is compared directly, not cyclically.
  fpath q(head);
  fpath q2(head);
  CHECK(q == q2);
  q2.pop_back();
  CHECK(!(q == q2));
  CHECK(q2.length() == 1);

  // The initial vertex of an empty path can be changed; pushes then
  // follow the graph from there.
  fpath r(ttg,0);
  const int v0 = 2;
  r.initial_vertex(v0);
  CHECK(r.initial_vertex() == v0 && r.final_vertex() == v0);
  r.push_back(0);
  CHECK(r.vertices()[1] == ttg.target_vertex(v0,0));
  r.clear();
  CHECK(r.zerolength() && r.initial_vertex() == v0);

  // ---- badwords content.
  const int maxplen = 4;
  auto pbad = ttauto::badwords(ttg,maxplen);
  CHECK((int)pbad.rows() == ttg.vertices());
  CHECK((int)pbad.columns() == maxplen);

  int nbad = 0;
  for (int v = 0; v < ttg.vertices(); ++v)
    for (int l = 0; l < maxplen; ++l)
      for (const fpath& w : pbad(v,l))
        {
          ++nbad;
          CHECK((int)w.length() == 2*(l+1));
          CHECK(w.initial_vertex() == v);
          CHECK(w.closed());
          const fpath h = w.subpath(l+1);
          CHECK(h.closed());
          CHECK(same_foldings(h*h,w));
          const Mat TM = h.transition_matrix();
          CHECK(ttauto::pattern_equal(TM,TM*TM));
          CHECK(w.transition_matrix() == TM*TM);
        }
  CHECK(nbad > 0);

  // The length-1 layer is complete: every self-loop whose matrix pattern
  // is idempotent appears exactly once.
  for (int v = 0; v < ttg.vertices(); ++v)
    {
      int nloops = 0;
      for (int b = 0; b < ttg.foldings(v); ++b)
        {
          if (ttg.target_vertex(v,b) != v) continue;
          fpath h(ttg,v);
          h.push_back(b);
          const Mat TM = h.transition_matrix();
          if (!ttauto::pattern_equal(TM,TM*TM)) continue;
          ++nloops;
          int nfound = 0;
          for (const fpath& w : pbad(v,0))
            if (same_foldings(w.subpath(1),h)) ++nfound;
          CHECK(nfound == 1);
        }
      CHECK(nloops == (int)pbad(v,0).size());
    }

  // Deterministic: recomputing gives the same table.
  auto pbad2 = ttauto::badwords(ttg,maxplen);
  for (int v = 0; v < ttg.vertices(); ++v)
    for (int l = 0; l < maxplen; ++l)
      {
        CHECK(pbad(v,l).size() == pbad2(v,l).size());
        auto i2 = pbad2(v,l).begin();
        for (const fpath& w : pbad(v,l)) CHECK(same_foldings(w,*i2++));
      }

  std::cout << "test_folding_path_and_badwords: OK (" << nbad << " bad words)\n";
  return 0;
}
