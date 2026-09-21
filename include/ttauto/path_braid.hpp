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

#ifndef TTAUTO_PATH_BRAID_HPP
#define TTAUTO_PATH_BRAID_HPP

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include "traintracks/braid.hpp"
#include "traintracks/embedding.hpp"
#include "traintracks/fold_map.hpp"
#include "traintracks/map.hpp"
#include "ttauto/folding_path.hpp"
#include "ttauto/ttfoldgraph.hpp"

namespace ttauto {

// The braid of a closed folding path.
//
// A closed path defines a homeomorphism of the punctured disc: perform the
// folds, then identify the final track with the initial one by the
// canonical numbering (doc/ttauto.tex).  Reading that homeomorphism as a
// braid needs the punctures to have positions, which is what a proper
// embedding gives (traintracks/embedding.hpp).  Each fold either leaves the
// track properly embedded, when the edge it moves ends up alongside an
// unpunctured multigon, or wraps the moved edge once around a puncture, and
// then the punctures have to be moved to undo it.  The braid is the product
// of those motions, finished by the rotation that brings the initial
// track's root monogon back to position 1 for the final identification.
//
// The braid is only defined up to conjugacy, since the puncture put at
// position 1 is a choice; the choice here is the root monogon of the
// canonical coding.

namespace detail {

// Fold index of a branch.  ttfoldgraph does not record it: build_graph
// numbers branches by their rank among the folds with a non-identity
// transition matrix and throws the fold index away, and a subgraph
// renumbers them again when it drops the branches that leave its block.
// So identify the fold by what it does, its transition matrix and the
// coding of the track it produces.  Two branches of one vertex can agree
// on both -- they are still different folds, at different cusps -- so
// break the tie by rank, which works for a subgraph too since taking a
// subgraph keeps the order of the branches it keeps.
template<class TrTr>
int fold_index_of_branch(const ttfoldgraph<TrTr>& ttg, const int v,
                         const int branch)
{
  const TrTr& tt = ttg.traintrack(v);
  const jlt::mathmatrix<int> TM(ttg.transition_matrix(v,branch).full());
  const typename TrTr::intVec target
    = ttg.traintrack(ttg.target_vertex(v,branch)).coding();

  int rank = 0;
  for (int b = 0; b < branch; ++b)
    if (ttg.transition_matrix(v,b).full() == TM
        && ttg.traintrack(ttg.target_vertex(v,b)).coding() == target) ++rank;

  int seen = 0;
  for (int f = 0; f < tt.foldings(); ++f)
    {
      if (traintracks::fold_transition_matrix(tt,f).full() != TM) continue;
      TrTr folded(tt);
      if (!folded.fold(f)) continue;
      if (folded.coding() != target) continue;
      if (seen == rank) return f;
      ++seen;
    }

  std::cerr << "No fold for branch " << branch << " at vertex " << v
            << " in ttauto::detail::fold_index_of_branch.\n";
  std::exit(1);
}

} // namespace detail

// Every braid a closed path can stand for.  Normally one; at a
// cyclically symmetric initial vertex the automaton identifies the final
// track with the initial one only up to the track's own symmetry, which is
// a root of the full twist, so there is one candidate per symmetry.  They
// differ by a central-ish factor and all have the same permutation.
template<class TrTr>
std::vector<traintracks::braidword>
folding_path_braids(const folding_path<TrTr>& p)
{
  using traintracks::braidword;
  using traintracks::outer_embedding;
  using traintracks::tt_embedding;
  using traintracks::ttnumbering;

  if (!p.closed())
    {
      std::cerr << "Path is not closed in ttauto::folding_path_braids.\n";
      std::exit(1);
    }

  const ttfoldgraph<TrTr>& ttg = p.graph();
  int v = p.initial_vertex();
  TrTr tt(ttg.traintrack(v));

  ttnumbering num = tt.numbering();
  tt_embedding emb = outer_embedding(num);
  const int np = emb.npunctures();

  braidword b(np);

  for (int k = 0; k < (int)p.length(); ++k)
    {
      const int branch = p.foldings()[k];
      const int f = detail::fold_index_of_branch(ttg,v,branch);

      traintracks::fold_map_data fm;
      const ttnumbering before = tt.numbering();
      if (!tt.fold_with_map(f,fm))
        {
          std::cerr << "Fold " << f << " failed at vertex " << v
                    << " in ttauto::folding_path_braids.\n";
          std::exit(1);
        }

      const traintracks::fold_swap fs =
        traintracks::fold_block_swap(before,fm,emb);
      b *= traintracks::fold_braid(before,fm,emb);

      emb = outer_embedding(fm.after,fs.next_cut_dart);
      v = ttg.target_vertex(v,branch);

      // Continue from the automaton's own copy of the track, not from the
      // one just folded.  The two have the same coding, but when the track
      // has a symmetry they need not be the same physical track, and it is
      // the automaton's copy that the next branch index refers to.
      tt = ttg.traintrack(v);
    }

  // The automaton identifies the final track with the initial one by the
  // canonical numbering, which puts the root monogon, prong 0, at position
  // 1.  Rotate until it is.  Rotating all the way round is the full twist,
  // which is central and which the automaton cannot see anyway, since
  // nothing here records a framing at the boundary of the disc; so take
  // whichever way round is shorter.
  int a = emb.position_of[0] - 1;
  if (2*a > np) a -= np;
  b *= (a >= 0 ? traintracks::rotation_braid(np,a)
                : traintracks::rotation_braid(np,-a).inverse());
  b.reduce();

  // The initial track's cyclic symmetry, if any, is a root of the full
  // twist and leaves that identification ambiguous by its powers.
  TrTr tt0(ttg.traintrack(p.initial_vertex()));
  const int order = tt0.is_cyclically_symmetric() + 1;

  std::vector<braidword> all;
  for (int j = 0; j < order; ++j)
    {
      braidword c = b;
      c *= traintracks::rotation_braid(np,j*(np/order));
      all.push_back(c.reduce());
    }
  return all;
}

namespace detail {

// Perron root of the path's transition matrix, by power iteration.  The
// matrix of a pseudo-Anosov path is primitive, so this converges.
inline double perron_root(const jlt::mathmatrix<int>& A)
{
  const int m = A.dim();
  std::vector<double> x(m,1.0), y(m);
  double lam = 1;
  for (int it = 0; it < 4000; ++it)
    {
      for (int i = 0; i < m; ++i)
        { y[i] = 0; for (int j = 0; j < m; ++j) y[i] += A(i,j)*x[j]; }
      double s = 0;
      for (int i = 0; i < m; ++i) s += y[i]*y[i];
      s = std::sqrt(s);
      if (!(s > 0)) return 0;
      for (int i = 0; i < m; ++i) y[i] /= s;
      double t1 = 0, t2 = 0;
      for (int i = 0; i < m; ++i)
        { double t = 0; for (int j = 0; j < m; ++j) t += A(i,j)*y[j];
          t1 += t*y[i]; t2 += y[i]*y[i]; }
      lam = t1/t2;
      x = y;
    }
  return lam;
}

} // namespace detail

// The braid of a closed folding path.
//
// The answer is checked: the braid's growth, computed from its action on
// Dynnikov coordinates and so by a route with nothing in common with the
// rest of this, must equal the Perron root of the path's transition
// matrix.  That is what picks the right candidate at a cyclically
// symmetric vertex, and `verified` reports whether any candidate passed.
// A braid that did not verify is still returned, since it is right up to
// the ambiguity, but it should not be trusted.
//
// The check only bites for a path whose transition matrix is primitive and
// whose gates are connected, since only then is the growth of the track
// map the dilatation of the braid; otherwise verified is left false.
template<class TrTr>
traintracks::braidword folding_path_braid(const folding_path<TrTr>& p,
                                          bool* verified = 0)
{
  const std::vector<traintracks::braidword> all = folding_path_braids(p);
  if (verified) *verified = false;

  const double lambda = detail::perron_root(p.transition_matrix());
  if (lambda > 1 + 1e-8)
    {
      for (std::size_t j = 0; j < all.size(); ++j)
        if (std::fabs(all[j].growth() - lambda) < 1e-6*lambda)
          {
            if (verified) *verified = true;
            return all[j];
          }
    }
  return all.front();
}

} // namespace ttauto

#endif // TTAUTO_PATH_BRAID_HPP
