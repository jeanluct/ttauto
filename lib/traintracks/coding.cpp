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
#include <iostream>
#include <map>
#include <vector>

#include "traintracks/coding.hpp"
#include "traintracks/traintrack.hpp"
#include "traintracks/util.hpp"

namespace traintracks {

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
  int p = pin, e = ein;

  if (mm.edges() == 1)
    {
      // This is an uncusped monogon so we won't recurse. Record it
      // here and continue. Uncusped monogons are marked with a
      // coding_block() sequence.
      coding_block(0,1,mm.label(),0,1).append_to(code);
      return;
    }

  do
    {
      // Prong number relative to entry prong into multigon.
      // The entry prong is labeled 0, and the other prongs
      // clockwise from 0 (anticlockwise if dir = -1).
      int prong = traintracks::mod(dir*(p-pin),mm.prongs());
      // Number of prongs in outgoing multigon.
      int nprongs = mm.prongs();
      // Label of the multigon.
      int label = mm.label();
      // Number of edges in the outgoing prong.
      int nedges = mm.edges(p);
      // The outgoing edge. Number anticlockwise if dir = -1.
      int edge = (dir == 1 ? e : nedges-1-e);

      // Record the block in the coding corresponding to this edge.
      coding_block(prong,nprongs,label,edge,nedges).append_to(code);

      // Find the next multigon down.
      int pout, eout;
      multigon* ed = mm.Edge(p,e)->target_multigon(&mm,pout,eout);

      // Don't recurse down the entry edge.
      if (!(p == pin && e == ein))
	{
	  recursive_coding(*ed,pout,eout,code,dir);
	}
      // Increment the edge and prong (decrement if dir = -1).
      mm.cycle_edges(p,e,dir);
    }
  while (!(p == pin && e == ein));
}

// Canonical prong/edge numbering rooted at uncusped monogon mono.
ttnumbering coding_engine::numbering(const traintrack& tt, int mono)
{
  if (mono < 0 || mono >= tt.multigons() || tt.Multigon(mono).edges() != 1)
    {
      std::cerr << "Not an uncusped monogon in traintracks::coding_engine::numbering.\n";
      std::exit(1);
    }

  ttnumbering num;
  num.start_monogon = mono;
  num.prong_number.resize(tt.multigons());
  for (int m = 0; m < tt.multigons(); ++m)
    num.prong_number[m].assign(tt.Multigon(m).prongs(),-1);

  // Monogon mono: prong number 0, edge number 0 oriented away from it.
  number_prongs(tt,mono,0,num);
  const multigon& mm0 = tt.Multigon(mono);
  num.edge_tail.push_back(num.prong_number[mono][0]);
  num.edge_head.push_back(-1);
  num.edge_ptr.push_back(mm0.Edge(0,0).get());

  int pmono, pemono;
  multigon* egmono = mm0.Edge(0,0)->target_multigon(&mm0,pmono,pemono);

  if (egmono->edges() > 1)
    {
      recursive_numbering(tt,*egmono,pmono,pemono,0,num);
    }
  else
    {
      // Degenerate two-monogon track.
      const int mj = tt.multigon_index(egmono);
      number_prongs(tt,mj,pmono,num);
      num.edge_head[0] = num.prong_number[mj][pmono];
    }

  // Signed main letters at each prong in slot order.
  {
    std::map<const edge*,int> number_of;
    for (int e = 0; e < num.nedges(); ++e) number_of[num.edge_ptr[e]] = e;
    num.prong_letters.assign(num.nprongs(),std::vector<int>());
    for (int q = 0; q < num.nprongs(); ++q)
      {
        const multigon& mm = tt.Multigon(num.prong[q].multigon);
        const int p = num.prong[q].prong;
        for (int e = 0; e < mm.edges(p); ++e)
          {
            auto it = number_of.find(mm.Edge(p,e).get());
            if (it == number_of.end())
              {
                std::cerr << "Edge not numbered in traintracks::coding_engine::numbering.\n";
                std::exit(1);
              }
            const int en = it->second;
            num.prong_letters[q].push_back(num.edge_tail[en] == q ? en+1 : -(en+1));
          }
      }
  }

  // Every prong must have been numbered exactly once, every edge must
  // have a head, and every cusp must have been met.
  if (num.nprongs() != tt.total_prongs() || num.nedges() != tt.edges() ||
      num.ncusps() != tt.cusps())
    {
      std::cerr << "Incomplete walk in traintracks::coding_engine::numbering.\n";
      std::exit(1);
    }
  for (int e = 0; e < num.nedges(); ++e)
    {
      if (num.edge_head[e] < 0)
        {
          std::cerr << "Edge without head in traintracks::coding_engine::numbering.\n";
          std::exit(1);
        }
    }

  return num;
}

void coding_engine::number_prongs(const traintrack& tt, int mi, int pin,
                                  ttnumbering& num)
{
  const multigon& mm = tt.Multigon(mi);
  const int k = mm.prongs();
  if (num.prong_number[mi][pin] >= 0)
    {
      std::cerr << "Multigon visited twice in traintracks::coding_engine::numbering.\n";
      std::exit(1);
    }
  for (int j = 0; j < k; ++j)
    {
      const int pp = traintracks::mod(pin+j,k);
      num.prong_number[mi][pp] = num.prong.size();
      ttnumbering::prong_info info;
      info.multigon = mi;
      info.prong = pp;
      info.nprongs = k;
      info.punctured = mm.punctured();
      num.prong.push_back(info);
    }
}

// The depth-first walk that defines the canonical order of edges, prongs
// and cusps (the one the coding also uses): the entry edge has already
// been numbered by the caller; number the other edges of mm in cycle_edges
// order, descending into every non-monogon target.
void coding_engine::recursive_numbering(const traintrack& tt,
                                        const multigon& mm,
                                        int pin,
                                        int ein,
                                        int entry_edge,
                                        ttnumbering& num)
{
  const int mi = tt.multigon_index(&mm);
  number_prongs(tt,mi,pin,num);
  num.edge_head[entry_edge] = num.prong_number[mi][pin];

  int p = pin, e = ein;
  // The cusp at the entry slot comes first in fold order.
  if (e < mm.edges(p)-1) num.cusp.push_back(std::make_pair(num.prong_number[mi][p],e));
  mm.cycle_edges(p,e);

  do
    {
      const edge* E = mm.Edge(p,e).get();
      const int en = num.nedges();
      num.edge_tail.push_back(num.prong_number[mi][p]);
      num.edge_head.push_back(-1);
      num.edge_ptr.push_back(E);
      // Cusp at this slot, before descending into the child.
      if (e < mm.edges(p)-1) num.cusp.push_back(std::make_pair(num.prong_number[mi][p],e));

      int pout, eout;
      multigon* ed = E->target_multigon(&mm,pout,eout);

      if (ed->edges() > 1)
        {
          recursive_numbering(tt,*ed,pout,eout,en,num);
        }
      else
        {
          const int mj = tt.multigon_index(ed);
          number_prongs(tt,mj,pout,num);
          num.edge_head[en] = num.prong_number[mj][pout];
        }
      mm.cycle_edges(p,e);
    }
  while (!(p == pin && e == ein));
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
  // If we're not labeling multigons, don't print the label, which
  // means the coding blocks are shorter.
  if (!traintrack::label_multigons) --print_length;

  coding_vec code = coding(tt,dir);
  for (int i = 0; i < (int)code.size(); i += coding_block::length)
    {
      if (traintrack::label_multigons)
        strm << code[i]+1 << code[i+1] << code[i+2]+1 << code[i+3]+1 << code[i+4];
      else
        strm << code[i]+1 << code[i+1] << code[i+3]+1 << code[i+4];
      if (i != ((int)code.size() - print_length)) strm << " ";
    }
  return strm;
}

} // namespace detail

// Public wrapper that delegates coding generation to coding_engine.
//
// ttnumbering methods
//

int ttnumbering::side_to(const int q) const
{
  const prong_info& info = prong[q];
  const int pp = traintracks::mod(info.prong+1,info.nprongs);
  return prong_number[info.multigon][pp];
}

bool ttnumbering::is_main(const int letter) const
{
  return (letter != 0 && std::abs(letter) <= nedges());
}

bool ttnumbering::is_side(const int letter) const
{
  const int a = std::abs(letter);
  return (a > nedges() && a <= nedges() + nprongs());
}

int ttnumbering::edge_of(const int letter) const
{
  if (!is_main(letter))
    {
      std::cerr << "Not a main letter in traintracks::ttnumbering::edge_of.\n";
      std::exit(1);
    }
  return std::abs(letter) - 1;
}

int ttnumbering::side_of(const int letter) const
{
  if (!is_side(letter))
    {
      std::cerr << "Not a side letter in traintracks::ttnumbering::side_of.\n";
      std::exit(1);
    }
  return std::abs(letter) - nedges() - 1;
}

int ttnumbering::tail_of(const int letter) const
{
  if (is_main(letter))
    {
      const int e = edge_of(letter);
      return (letter > 0 ? edge_tail[e] : edge_head[e]);
    }
  const int q = side_of(letter);
  return (letter > 0 ? q : side_to(q));
}

void ttnumbering::fold_cusp(const int f, int& m, int& p, int& slot) const
{
  if (f < 0 || f >= foldings())
    {
      std::cerr << "Illegal folding index in traintracks::ttnumbering::fold_cusp.\n";
      std::exit(1);
    }
  const std::pair<int,int>& c = cusp[f/2];
  m = prong[c.first].multigon;
  p = prong[c.first].prong;
  slot = c.second;
}

int ttnumbering::side_from_prev(const int q) const
{
  const prong_info& info = prong[q];
  const int pp = traintracks::mod(info.prong-1,info.nprongs);
  return prong_number[info.multigon][pp];
}

std::vector<int> ttnumbering::directions_at_prong(const int q) const
{
  std::vector<int> d;
  if (prong[q].punctured) d.push_back(-side_letter(side_from_prev(q)));
  d.insert(d.end(),prong_letters[q].begin(),prong_letters[q].end());
  if (prong[q].punctured) d.push_back(side_letter(q));
  return d;
}

bool ttnumbering::operator==(const ttnumbering& o) const
{
  // Edge identity (edge_ptr) is deliberately excluded: two copies of the
  // same track have the same numbering.
  if (start_monogon != o.start_monogon) return false;
  if (prong.size() != o.prong.size()) return false;
  for (std::size_t q = 0; q < prong.size(); ++q)
    {
      if (prong[q].multigon != o.prong[q].multigon ||
          prong[q].prong != o.prong[q].prong ||
          prong[q].nprongs != o.prong[q].nprongs ||
          prong[q].punctured != o.prong[q].punctured) return false;
    }
  return (prong_number == o.prong_number &&
          edge_tail == o.edge_tail && edge_head == o.edge_head &&
          prong_letters == o.prong_letters && cusp == o.cusp);
}

bool ttnumbering::same_labels(const ttnumbering& o) const
{
  if (prong.size() != o.prong.size()) return false;
  for (std::size_t q = 0; q < prong.size(); ++q)
    {
      if (prong[q].nprongs != o.prong[q].nprongs ||
          prong[q].punctured != o.prong[q].punctured) return false;
    }
  return (edge_tail == o.edge_tail && edge_head == o.edge_head &&
          prong_letters == o.prong_letters && cusp == o.cusp);
}

std::ostream& ttnumbering::print(std::ostream& strm) const
{
  strm << "prongs (number: multigon,prong; k; punctured):\n";
  for (int q = 0; q < nprongs(); ++q)
    {
      strm << "  " << q << ": (" << prong[q].multigon << "," << prong[q].prong
           << "); k=" << prong[q].nprongs
           << (prong[q].punctured ? "; punctured" : "")
           << "; side " << side_letter(q) << " -> prong " << side_to(q) << "\n";
    }
  strm << "cusps (fold index/2: prong number, slot):\n";
  for (int c = 0; c < ncusps(); ++c)
    strm << "  " << c << ": prong " << cusp[c].first << " slot " << cusp[c].second << "\n";
  strm << "edges (number: tail prong -> head prong):\n";
  for (int e = 0; e < nedges(); ++e)
    {
      strm << "  " << main_letter(e) << ": " << edge_tail[e] << " -> "
           << edge_head[e] << "\n";
    }
  return strm;
}

ttnumbering traintrack::numbering() const
{
  require_normalised("traintrack::numbering");
  return detail::coding_engine::numbering(*this,0);
}

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
