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

// The automaton lifted to (vertex, permutation) states, which the
// completion bound of issue #23 needs (ttauto/completion_bound.hpp).
//
// Three checks.  Every fold matrix reassembles as P (I + e_{a'b}) from
// its decomposition.  The n=3 automaton, whose two folds carry no
// permutation, gives R = {identity} and two unit positions.  And on the
// n=4 automata, R(v) and U(v) agree exactly with a brute-force
// enumeration of the paths from v back to the start, made independently
// by multiplying the fold matrices themselves.

#include <iostream>
#include <set>
#include <vector>
#include <jlt/mathmatrix.hpp>
#include "check.hpp"
#include "traintracks/build.hpp"
#include "traintracks/traintrack.hpp"
#include "ttauto/completion_bound.hpp"
#include "ttauto/ttfoldgraph.hpp"

using traintracks::traintrack;
typedef ttauto::ttfoldgraph<traintrack> ttgraph;
typedef jlt::mathmatrix<int> Mat;
using ttauto::perm;
using ttauto::position;

static Mat perm_matrix(const perm& p)
{
  Mat P(p.size(),p.size());
  for (std::size_t c = 0; c < p.size(); ++c) P(p[c],c) = 1;
  return P;
}

// Read a permutation matrix back as a map, or exit if it is not one.
static perm as_perm(const Mat& P)
{
  const int n = P.dim();
  perm p(n,-1);
  for (int c = 0; c < n; ++c)
    for (int i = 0; i < n; ++i)
      if (P(i,c) == 1) { CHECK(p[c] == -1); p[c] = i; }
      else CHECK(P(i,c) == 0);
  for (int c = 0; c < n; ++c) CHECK(p[c] >= 0);
  return p;
}

// All paths from v back to v0 of length at most len, by depth-first
// search.  Along a path from v the frame permutation S (a matrix) is the
// product of the permutation parts so far, found by stripping the extra 1
// from each fold matrix; a fold F = P + e_{ab} contributes the unit
// S^{-1} P^{-1} e_{ab} S, located from the matrix product itself.
static void brute(const ttgraph& g, const int v0, const int w, const Mat& S,
                  const int len, std::vector<std::pair<Mat,std::set<position> > >& paths,
                  std::set<position> units)
{
  if (w == v0) paths.push_back(std::make_pair(S,units));
  if (len == 0) return;
  for (int f = 0; f < g.foldings(w); ++f)
    {
      const Mat F = g.transition_matrix(w,f).full();
      const int n = F.dim();
      // Strip the extra 1: it is where a row sum is 2 and a column sum is 2.
      Mat P(F);
      for (int a = 0; a < n; ++a)
        for (int b = 0; b < n; ++b)
          {
            int rs = 0, cs = 0;
            for (int k = 0; k < n; ++k) { rs += F(a,k); cs += F(k,b); }
            if (F(a,b) == 1 && rs == 2 && cs == 2)
              {
                Mat T(F); T(a,b) = 0;
                bool isperm = true;
                for (int i = 0; i < n && isperm; ++i)
                  {
                    int r = 0, c = 0;
                    for (int k = 0; k < n; ++k) { r += T(i,k); c += T(k,i); }
                    if (r != 1 || c != 1) isperm = false;
                  }
                if (isperm) P = T;
              }
          }
      const Mat E = F - P;                         // the extra 1, e_ab
      const Mat Pinv = perm_matrix(ttauto::inverse(as_perm(P)));
      const Mat Sinv = perm_matrix(ttauto::inverse(as_perm(S)));
      const Mat U = Sinv * Pinv * E * S;           // S^{-1} P^{-1} e_ab S
      std::set<position> units2(units);
      int found = 0;
      for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
          if (U(i,j) != 0) { CHECK(U(i,j) == 1); units2.insert(position(i,j)); ++found; }
      CHECK(found == 1);
      brute(g,v0,g.target_vertex(w,f),P*S,len-1,paths,units2);
    }
}

int main()
{
  using std::cout;
  using std::endl;

  // Every fold matrix is P (I + e_{a'b}).
  long long nfolds = 0;
  for (int n = 3; n <= 5; ++n)
    {
      jlt::vector<traintrack> ttv = traintracks::build_traintrack_list(n);
      for (int s = 0; s < (int)ttv.size(); ++s)
        {
          const ttgraph g(ttv[s]);
          for (int v = 0; v < g.vertices(); ++v)
            for (int f = 0; f < g.foldings(v); ++f)
              {
                const ttauto::fold_parts fp = ttauto::decompose_fold(g,v,f);
                const int d = g.edges();
                Mat IE(d,d);
                for (int i = 0; i < d; ++i) IE(i,i) = 1;
                IE(fp.unit.first,fp.unit.second) += 1;
                CHECK(perm_matrix(fp.p) * IE == g.transition_matrix(v,f).full());
                ++nfolds;
              }
        }
    }
  cout << nfolds << " folds for n = 3..5 reassemble as P (I + e)" << endl;

  // Lemma 5.5 along every path of length up to 6 from every vertex, n = 4
  // and 5: the frame S and the units from advance_frame() rebuild the
  // path's matrix, S (I + e_k) ... (I + e_1) = F_k ... F_1.
  long long npaths = 0;
  for (int n = 4; n <= 5; ++n)
    {
      jlt::vector<traintrack> ttv = traintracks::build_traintrack_list(n);
      for (int s = 0; s < (int)ttv.size(); ++s)
        {
          const ttgraph g(ttv[s]);
          const int d = g.edges();
          Mat I(d,d);
          for (int i = 0; i < d; ++i) I(i,i) = 1;
          for (int v = 0; v < g.vertices(); ++v)
            {
              // Depth-first over paths from v, carrying the true product
              // M, the frame s and N = the product of the (I + e) factors.
              struct node { int w; Mat M; perm s; Mat N; int len; };
              std::vector<node> stack(1,node{v,I,ttauto::identity_perm(d),I,0});
              while (!stack.empty())
                {
                  const node x = stack.back(); stack.pop_back();
                  CHECK(perm_matrix(x.s) * x.N == x.M);
                  ++npaths;
                  if (x.len == 6) continue;
                  for (int f = 0; f < g.foldings(x.w); ++f)
                    {
                      node y{g.target_vertex(x.w,f),
                             g.transition_matrix(x.w,f).full() * x.M,
                             x.s, Mat(), x.len+1};
                      const position u =
                        ttauto::advance_frame(ttauto::decompose_fold(g,x.w,f),y.s);
                      Mat IE(I);
                      IE(u.first,u.second) += 1;
                      y.N = IE * x.N;
                      stack.push_back(y);
                    }
                }
            }
        }
    }
  cout << npaths << " paths for n = 4, 5 rebuild as S (I + e) ... (I + e)" << endl;

  // n = 3: no permutations, two folds.
  {
    jlt::vector<traintrack> ttv = traintracks::build_traintrack_list(3);
    const ttgraph g(ttv[0]);
    CHECK(g.vertices() == 1);
    const ttauto::completion_sets<traintrack> cs(g,0);
    CHECK(cs.completion_perms(0).size() == 1);
    CHECK(*cs.completion_perms(0).begin() == ttauto::identity_perm(2));
    CHECK(cs.unit_positions(0).size() == 2);
    CHECK(ttauto::fold_group_order(g) == 1);
  }

  // n = 4: the sets against brute force, from every start vertex.  The
  // brute-force length grows by two until the sets stop changing, which
  // is when every (vertex, permutation) state has been reached.
  jlt::vector<traintrack> ttv = traintracks::build_traintrack_list(4);
  for (int s = 0; s < (int)ttv.size(); ++s)
    {
      const ttgraph g(ttv[s]);
      for (int v0 = 0; v0 < g.vertices(); ++v0)
        {
          const ttauto::completion_sets<traintrack> cs(g,v0);
          int len = 0;
          for (int v = 0; v < g.vertices(); ++v)
            {
              std::set<perm> R, Rprev;
              std::set<position> U, Uprev;
              for (len = 4; len <= 24; len += 2)
                {
                  Rprev = R; Uprev = U;
                  R.clear(); U.clear();
                  std::vector<std::pair<Mat,std::set<position> > > paths;
                  Mat I(g.edges(),g.edges());
                  for (int i = 0; i < g.edges(); ++i) I(i,i) = 1;
                  brute(g,v0,v,I,len,paths,std::set<position>());
                  for (const auto& pu : paths)
                    {
                      R.insert(as_perm(pu.first));
                      U.insert(pu.second.begin(),pu.second.end());
                    }
                  if (len > 4 && R == Rprev && U == Uprev && !R.empty()) break;
                  if (len > 8 && R.empty()) break;   // v cannot return to v0
                }
              CHECK(len <= 24);
              CHECK_MSG(R == cs.completion_perms(v),
                        "n=4 stratum " << s << " start " << v0 << " vertex " << v);
              CHECK_MSG(U == cs.unit_positions(v),
                        "n=4 stratum " << s << " start " << v0 << " vertex " << v);
            }
          cout << "n=4 stratum " << s << ", start " << v0 << ": "
               << cs.lifted_states() << " lifted states, sets match brute force"
               << endl;
        }
    }

  cout << "\ntest_completion_bound: all checks passed" << endl;
  return 0;
}
