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

#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <unordered_set>

#include "traintracks/coding.hpp"
#include "traintracks/traintrack.hpp"
#include "traintracks/util.hpp"

namespace traintracks {

namespace {

struct coding_visit_key
{
  // Directed entry state into a multigon during coding traversal.
  const multigon* mm;
  int pin;
  int ein;
  int dir;

  bool operator==(const coding_visit_key& other) const
  {
    return (mm == other.mm && pin == other.pin &&
            ein == other.ein && dir == other.dir);
  }
};

struct coding_visit_key_hash
{
  std::size_t operator()(const coding_visit_key& key) const
  {
    std::size_t h = static_cast<std::size_t>(reinterpret_cast<std::uintptr_t>(key.mm));
    // Boost-style hash_combine mix: golden-ratio constant plus shifted
    // accumulator bits to reduce collisions when folding in each field.
    h ^= static_cast<std::size_t>(key.pin + 0x9e3779b9 + (h << 6) + (h >> 2));
    h ^= static_cast<std::size_t>(key.ein + 0x9e3779b9 + (h << 6) + (h >> 2));
    h ^= static_cast<std::size_t>(key.dir + 0x9e3779b9 + (h << 6) + (h >> 2));
    return h;
  }
};

void recursive_coding_guarded(
  const multigon& mm,
  int pin,
  int ein,
  detail::coding_engine::coding_vec& code,
  int dir,
  std::unordered_set<coding_visit_key,coding_visit_key_hash>& visited)
{
  // If this directed entry state was already expanded, stop this branch.
  const coding_visit_key key{&mm,pin,ein,dir};
  if (!visited.insert(key).second)
    {
      return;
    }

  int p = pin, e = ein;

  if (mm.edges() == 1)
    {
      // Uncusped monogon: emit its marker block and terminate this branch.
      detail::coding_block(0,1,mm.label(),0,1).append_to(code);
      return;
    }

  do
    {
      int prong = traintracks::mod(dir*(p-pin),mm.prongs());
      int nprongs = mm.prongs();
      int label = mm.label();
      int nedges = mm.edges(p);
      int edge = (dir == 1 ? e : nedges-1-e);

      detail::coding_block(prong,nprongs,label,edge,nedges).append_to(code);

      int pout, eout;
      multigon* ed = mm.Edge(p,e)->target_multigon(&mm,pout,eout);

      // Skip recursing back through the incoming edge to this multigon.
      if (!(p == pin && e == ein))
        {
          recursive_coding_guarded(*ed,pout,eout,code,dir,visited);
        }
      mm.cycle_edges(p,e,dir);
    }
  while (!(p == pin && e == ein));
}

} // namespace

// Canonicalize multigon order and coding root for invariant comparisons.
void traintrack::normalise()
{
  // Canonicalization has two phases:
  // 1) normalize each multigon's local representation;
  // 2) order multigons globally by multigon strict-order (sort), then rotate
  //    the initial uncusped monogon to the lexicographically minimising
  //    coding start (minimise_coding).
  for (int m = 0; m < multigons(); ++m)
    {
      Multigon(m).normalise();
    }
  sort();
  detail::coding_engine::minimise_coding(*this);
  isnormalised = true;
}

namespace detail {

// Emit coding rooted at uncusped monogon mono in orientation dir.
coding_engine::coding_vec
coding_engine::coding_from_monogon(const traintrack& tt,
				   int mono,
				   int dir)
{
  // Build a canonical coding sequence beginning at uncusped monogon mono,
  // then traverse adjacent multigons in orientation dir.
  if (tt.Multigon(mono).edges() > 1)
    {
      std::cerr << "Not an uncusped monogon in traintrack::traintrack::coding_from_monogon.\n";
      std::exit(1);
    }

  coding_vec code;

  // Start the coding vector with the monogon.
  // Uncusped monogons are marked with a coding_block() sequence.
  coding_block(0,1,tt.Multigon(mono).label(),0,1).append_to(code);

  // Recurse down and compute coding.
  // Start by finding the multigon the edge is attached to, and which prong.
  int pmono, pemono;
  multigon* egmono =
    tt.Multigon(mono).Edge(0,0)->target_multigon(&tt.Multigon(mono),pmono,pemono);

  recursive_coding(*egmono,pmono,pemono,code,dir);
  return code;
}

// Normalize monogon root by lexicographically minimising coding.
coding_engine::coding_vec coding_engine::minimise_coding(traintrack& tt)
{
  // Among uncusped monogon starts, choose lexicographically minimal coding
  // and rotate track so the minimiser is in slot 0.
  int mono = 0, monomin = 0;
  coding_vec codemin = coding_from_monogon(tt,mono);

  // Loop over monogons with only one edge (uncusped).
  while (tt.Multigon(++mono).edges() <= 1)
    {
      coding_vec code = coding_from_monogon(tt,mono);
      if (code < codemin) { codemin = code; monomin = mono; }
    }
  // Move the minimising monogon to the first slot.
  tt.swap(0,monomin);

  return codemin;
}

// Read-only canonical coding over all uncusped monogon starts.
coding_engine::coding_vec coding_engine::coding(const traintrack& tt, int dir)
{
  // Read-only canonical coding: find the minimal coding over all uncusped
  // monogon starts, preserving track order.
  tt.require_normalised("traintrack::coding");

  int mono = 0;
  coding_vec codemin = coding_from_monogon(tt,mono,dir);

  // Loop over monogons with only one edge (uncusped).
  while (tt.Multigon(++mono).edges() <= 1)
    {
      coding_vec code = coding_from_monogon(tt,mono,dir);
      if (code < codemin) { codemin = code; }
    }

  return codemin;
}

// DFS over edge-neighbor graph emitting fixed-length coding blocks.
void coding_engine::recursive_coding(const multigon& mm,
				     int pin,
				     int ein,
				     coding_vec& code,
				     int dir)
{
  // Depth-first edge traversal that emits coding blocks in canonical order.
  // Guard revisits of the same directed entry state to avoid runaway recursion
  // on cyclic traversals.
  std::unordered_set<coding_visit_key,coding_visit_key_hash> visited;
  recursive_coding_guarded(mm,pin,ein,code,dir,visited);
}

// Detect cyclic coding matches and recover induced branch permutation.
mathmatrix_permplus1 coding_engine::cyclic_symmetry(traintrack& tt)
{
  // Detect cyclic symmetry by comparing codings from uncusped monogon starts,
  // then derive the induced branch permutation from transported weights.
  tt.require_normalised("traintrack::cyclic_symmetry");

  int mono = 0, nmatch = 0;
  coding_vec code = coding_from_monogon(tt,mono);
  mathmatrix_permplus1 perm(jlt::identity_matrix<int>(tt.edges()));

  // Loop over monogons with only one edge (uncusped) and save codings.
  while (tt.Multigon(++mono).edges() <= 1)
    {
      if (code == coding_from_monogon(tt,mono))
	{
	  ++nmatch;
	  // The first time we have a match, compute the permutation matrix.
	  if (nmatch == 1)
	    {
	      traintrack::dblVec w(tt.edges()), w2(tt.edges());
	      jlt::mathmatrix<int> M(tt.edges(),tt.edges());
	      for (int i = 0; i < tt.edges(); ++i)
		{
		  // Set initial weights.
		  w[i] = 1;
		  tt.weights(w.begin());
		  w[i] = 0;
		  // Find where the weights are in terms of the new labels.
		  w2 = tt.weights(mono);
		  int j = 0;
		  for (j = 0; j < tt.edges(); ++j) if (w2[j] != 0) break;
		  // Make permutation matrix.
		  M(j,i) = 1;
		}
	      perm = mathmatrix_permplus1(M);
	    }
	}
    }

  if (nmatch)
    {
      // Sanity check: perm should be such that perm^order=id, where
      // order = nmatch+1. But it should not be id for a smaller power.
      if (perm.order() != nmatch+1)
	{
	  std::cerr << "Bad permutation in traintrack::traintrack::cyclic_symmetry().\n";
	  std::exit(1);
	}
    }

  return perm;
}

// Stream canonical coding as contiguous printable blocks.
std::ostream& coding_engine::print_coding(const traintrack& tt,
					  std::ostream& strm,
					  int dir)
{
  // Pretty-printer for coding blocks in canonical orientation dir.
  int print_length = coding_block::length;
  // If we're not labeling multiprongs, don't print the label, which
  // means the coding blocks are shorter.
  if (!traintrack::label_multiprongs) --print_length;

  coding_vec code = coding(tt,dir);
  for (int i = 0; i < (int)code.size(); i += coding_block::length)
    {
      if (traintrack::label_multiprongs)
        strm << code[i]+1 << code[i+1] << code[i+2]+1 << code[i+3]+1 << code[i+4];
      else
        strm << code[i]+1 << code[i+1] << code[i+3]+1 << code[i+4];
      if (i != ((int)code.size() - print_length)) strm << " ";
    }
  return strm;
}

} // namespace detail

// Public wrapper that delegates coding generation to coding_engine.
traintrack::intVec traintrack::coding(const int dir) const
{
  return detail::coding_engine::coding(*this,dir);
}

// Public wrapper that delegates symmetry detection to coding_engine.
mathmatrix_permplus1 traintrack::cyclic_symmetry()
{
  return detail::coding_engine::cyclic_symmetry(*this);
}

// Public wrapper that delegates coding formatting to coding_engine.
std::ostream& traintrack::print_coding(std::ostream& strm,
				       const int dir) const
{
  return detail::coding_engine::print_coding(*this,strm,dir);
}

} // namespace traintracks
