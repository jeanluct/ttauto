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

#ifndef TTAUTO_COMPLETION_BOUND_HPP
#define TTAUTO_COMPLETION_BOUND_HPP

// What the completion bound of issue #23 needs to know about an automaton
// (devel/iss023/pruning_bounds.tex, Section 5.2).
//
// A fold matrix is F = P + e_ab, a permutation plus a single 1.  Written
// as F = P (I + e_{a'b}) with P(a') = a, the permutations of a path can be
// moved to the left (Lemma 5.5), so a path's matrix is Q N: Q the product
// of its folds' permutation parts, N a product of factors I + e, one per
// fold, their units written in the labels of the path's first vertex.
//
// For a search from a start vertex v0, and a prefix ending at v with
// matrix Q_A N_A, a completion (a path from v back to v0) contributes a
// permutation product r and a set of units u, both written relative to v.
// The closed path then has total permutation r Q_A, and the units become
// Q_A^{-1} u Q_A in the labels of v0.  This header computes, for every v,
//
//   R(v): the permutation products r of all paths from v back to v0, and
//   U(v): the positions of all the units u met on such paths,
//
// which do not depend on the prefix.  They are found on the automaton
// lifted to (vertex, permutation) states.
//
// Permutations are maps c -> p[c], with P e_c = e_{p[c]}.

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <set>
#include <utility>
#include <vector>
#include "ttauto/ttfoldgraph.hpp"

namespace ttauto {

typedef std::vector<int> perm;
typedef std::pair<int,int> position;	// (row, column) of a single 1

inline perm identity_perm(const int n)
{
  perm p(n);
  for (int c = 0; c < n; ++c) p[c] = c;
  return p;
}

// p o q: first q, then p.
inline perm compose(const perm& p, const perm& q)
{
  perm r(q.size());
  for (std::size_t c = 0; c < q.size(); ++c) r[c] = p[q[c]];
  return r;
}

inline perm inverse(const perm& p)
{
  perm r(p.size());
  for (std::size_t c = 0; c < p.size(); ++c) r[p[c]] = (int)c;
  return r;
}

// The fold F = P (I + e_{a'b}) as its permutation P and the position
// (a',b) of its unit.
struct fold_parts
{
  perm p;
  position unit;
};

template<class TrTr>
fold_parts decompose_fold(const ttfoldgraph<TrTr>& ttg, const int v,
			  const int f)
{
  const typename ttfoldgraph<TrTr>::Matpp1& F = ttg.transition_matrix(v,f);
  if (F.is_perm())
    {
      std::cerr << "A fold matrix without its extra 1 in ";
      std::cerr << "ttauto::decompose_fold\n";
      std::exit(1);
    }
  // full() has F(j, row_perm()[j]) = 1, so P(row_perm()[j]) = j, which is
  // column_perm(); and P^{-1}(a) = row_perm()[a].
  fold_parts fp;
  const int n = (int)F.dim();
  fp.p.resize(n);
  for (int c = 0; c < n; ++c) fp.p[c] = F.column_perm()[c];
  fp.unit = position(F.row_perm()[F.plus1_row()],F.plus1_col());
  return fp;
}

// One step of Lemma 5.5.  A path so far has matrix S N, S the product of
// its permutation parts (the map s).  Appending the fold F = P (I + e')
// gives F S N = (P S) (I + S^{-1} e' S) N: return the position of the new
// unit S^{-1} e' S, and replace s by P o s.
inline position advance_frame(const fold_parts& fp, perm& s)
{
  const perm si = inverse(s);
  const position u(si[fp.unit.first],si[fp.unit.second]);
  s = compose(fp.p,s);
  return u;
}

template<class TrTr>
class completion_sets
{
public:
  completion_sets(const ttfoldgraph<TrTr>& ttg_, const int v0_)
    : ttg(ttg_), v0(v0_), n(ttg_.edges()), R(ttg_.vertices()),
      U(ttg_.vertices()), nlifted(0)
  {
    const int V = ttg.vertices();
    folds.resize(V);
    for (int v = 0; v < V; ++v)
      for (int f = 0; f < ttg.foldings(v); ++f)
	folds[v].push_back(decompose_fold(ttg,v,f));

    find_R();
    find_U();
  }

  int start_vertex() const { return v0; }

  // Permutation products of the paths from v back to the start.
  const std::set<perm>& completion_perms(const int v) const { return R[v]; }

  // Positions of the units met on those paths, relative to v.
  const std::set<position>& unit_positions(const int v) const { return U[v]; }

  // Number of (vertex, permutation) states found, over all vertices: the
  // size of the automaton lifted to permutations, restricted to states
  // that can return to the start.
  long long lifted_states() const { return nlifted; }

  const fold_parts& fold(const int v, const int f) const { return folds[v][f]; }

private:
  const ttfoldgraph<TrTr>& ttg;
  const int v0;
  const int n;
  std::vector<std::vector<fold_parts> > folds;
  std::vector<std::set<perm> > R;
  std::vector<std::set<position> > U;
  long long nlifted;

  // Backwards from (v0, identity): if r is the product of a path from v
  // to v0, and the fold f at w leads to v with permutation P, then the
  // path w -> v -> ... -> v0 has product r o P, the fold coming first.
  void find_R()
  {
    const int V = ttg.vertices();
    std::vector<std::vector<std::pair<int,int> > > into(V);   // (w, f)
    for (int w = 0; w < V; ++w)
      for (int f = 0; f < ttg.foldings(w); ++f)
	into[ttg.target_vertex(w,f)].push_back(std::make_pair(w,f));

    std::vector<std::pair<int,perm> > queue;
    R[v0].insert(identity_perm(n));
    queue.push_back(std::make_pair(v0,identity_perm(n)));
    for (std::size_t k = 0; k < queue.size(); ++k)
      {
	const int v = queue[k].first;
	const perm r = queue[k].second;
	for (const std::pair<int,int>& wf : into[v])
	  {
	    const perm r2 = compose(r,folds[wf.first][wf.second].p);
	    if (R[wf.first].insert(r2).second)
	      queue.push_back(std::make_pair(wf.first,r2));
	  }
      }
    nlifted = (long long)queue.size();
  }

  // Forwards from (v, identity) for each v: at a state (w, s), with s the
  // permutation product so far, the fold f at w contributes the unit
  // s^{-1} e' s, at position (s^{-1}(a'), s^{-1}(b)), provided the fold
  // leads somewhere from which the start can be reached.
  void find_U()
  {
    const int V = ttg.vertices();
    std::vector<bool> back(V,false);
    for (int v = 0; v < V; ++v) back[v] = !R[v].empty();

    for (int v = 0; v < V; ++v)
      {
	if (!back[v]) continue;
	std::set<std::pair<int,perm> > seen;
	std::vector<std::pair<int,perm> > queue;
	queue.push_back(std::make_pair(v,identity_perm(n)));
	seen.insert(queue.back());
	// Once every off-diagonal position is in, U(v) cannot grow.
	const std::size_t all = (std::size_t)n*(n-1);
	for (std::size_t k = 0; k < queue.size() && U[v].size() < all; ++k)
	  {
	    const int w = queue[k].first;
	    for (int f = 0; f < ttg.foldings(w); ++f)
	      {
		const int w2 = ttg.target_vertex(w,f);
		if (!back[w2]) continue;
		perm s2 = queue[k].second;
		U[v].insert(advance_frame(folds[w][f],s2));
		std::pair<int,perm> next(w2,s2);
		if (seen.insert(next).second) queue.push_back(next);
	      }
	  }
      }
  }
};

// Order of the group generated by the permutation parts of all the folds
// of the automaton.
template<class TrTr>
long long fold_group_order(const ttfoldgraph<TrTr>& ttg)
{
  const int n = ttg.edges();
  std::vector<perm> gens;
  for (int v = 0; v < ttg.vertices(); ++v)
    for (int f = 0; f < ttg.foldings(v); ++f)
      gens.push_back(decompose_fold(ttg,v,f).p);
  std::set<perm> group;
  std::vector<perm> queue(1,identity_perm(n));
  group.insert(queue[0]);
  for (std::size_t k = 0; k < queue.size(); ++k)
    for (const perm& g : gens)
      {
	const perm h = compose(g,queue[k]);
	if (group.insert(h).second) queue.push_back(h);
      }
  return (long long)group.size();
}

} // namespace ttauto

#endif // TTAUTO_COMPLETION_BOUND_HPP
