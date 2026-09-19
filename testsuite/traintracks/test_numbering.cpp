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

// Canonical prong/edge numbering (ttnumbering): structural invariants and
// agreement with the edge order used by weights().

#include <iostream>
#include <set>
#include <vector>
#include "check.hpp"
#include "traintracks/build.hpp"
#include "traintracks/coding.hpp"
#include "traintracks/traintrack.hpp"
#include "ttauto/ttfoldgraph.hpp"

using traintracks::traintrack;
using traintracks::ttnumbering;

static void check_numbering(const traintrack& tt, const char* what)
{
  const ttnumbering N = tt.numbering();

  CHECK_MSG(N.nprongs() == tt.total_prongs(), what);
  CHECK_MSG(N.nedges() == tt.edges(), what);
  CHECK_MSG(N.start_monogon == 0, what);
  CHECK_MSG(N.prong[0].multigon == 0 && N.prong[0].prong == 0, what);

  // prong_number is a bijection onto 0..nprongs-1 consistent with prong[].
  std::set<int> seen;
  for (int m = 0; m < tt.multigons(); ++m)
    {
      CHECK_MSG((int)N.prong_number[m].size() == tt.Multigon(m).prongs(), what);
      for (int p = 0; p < tt.Multigon(m).prongs(); ++p)
        {
          const int q = N.prong_number[m][p];
          CHECK_MSG(q >= 0 && q < N.nprongs(), what);
          CHECK_MSG(seen.insert(q).second, what);
          CHECK_MSG(N.prong[q].multigon == m && N.prong[q].prong == p, what);
          CHECK_MSG(N.prong[q].nprongs == tt.Multigon(m).prongs(), what);
          CHECK_MSG(N.prong[q].punctured == tt.Multigon(m).punctured(), what);
        }
    }

  // Sides stay inside their multigon; prongs of one multigon are numbered
  // consecutively in cycle order; monogon sides are loops.
  for (int q = 0; q < N.nprongs(); ++q)
    {
      const int q2 = N.side_to(q);
      CHECK_MSG(N.prong[q2].multigon == N.prong[q].multigon, what);
      if (N.prong[q].nprongs == 1) CHECK_MSG(q2 == q, what);
      else CHECK_MSG(q2 != q, what);
    }

  // Edges: valid distinct endpoints; letter helpers consistent.
  for (int e = 0; e < N.nedges(); ++e)
    {
      const int t = N.edge_tail[e], h = N.edge_head[e];
      CHECK_MSG(t >= 0 && t < N.nprongs() && h >= 0 && h < N.nprongs(), what);
      CHECK_MSG(N.prong[t].multigon != N.prong[h].multigon, what);
      CHECK_MSG(N.is_main(N.main_letter(e)) && !N.is_side(N.main_letter(e)), what);
      CHECK_MSG(N.tail_of(N.main_letter(e)) == t && N.head_of(N.main_letter(e)) == h, what);
      CHECK_MSG(N.tail_of(-N.main_letter(e)) == h && N.head_of(-N.main_letter(e)) == t, what);
      CHECK_MSG(N.edge_of(N.main_letter(e)) == e, what);
      // The edge object attached at the tail prong slot list contains it.
      const int m = N.prong[t].multigon, p = N.prong[t].prong;
      bool found = false;
      for (int s = 0; s < tt.Multigon(m).edges(p); ++s)
        if (tt.Multigon(m).Edge(p,s).get() == N.edge_ptr[e]) found = true;
      CHECK_MSG(found, what);
    }
  for (int q = 0; q < N.nprongs(); ++q)
    {
      CHECK_MSG(N.is_side(N.side_letter(q)) && !N.is_main(N.side_letter(q)), what);
      CHECK_MSG(N.side_of(N.side_letter(q)) == q, what);
      CHECK_MSG(N.tail_of(N.side_letter(q)) == q, what);
      CHECK_MSG(N.head_of(N.side_letter(q)) == N.side_to(q), what);
    }

  // Edge order agrees with weights(): set weight e+1 on the e-th weight
  // slot and read it back through the edge objects.
  {
    traintrack ttw(tt);
    const ttnumbering Nw = ttw.numbering();
    traintrack::dblVec wv(Nw.nedges());
    for (int e = 0; e < Nw.nedges(); ++e) wv[e] = e + 1;
    ttw.weights(wv.begin());
    for (int e = 0; e < Nw.nedges(); ++e)
      CHECK_MSG(Nw.edge_ptr[e]->weight() == e + 1, what);
  }

  // A copy has the same numbering, and so does a re-normalised copy.
  traintrack tt2(tt);
  CHECK_MSG(tt2.numbering() == N, what);
  tt2.normalise();
  CHECK_MSG(tt2.numbering() == N, what);
}

// ---- Oracles: the walks that weights() and fold() used before they were
// expressed through the numbering, written against the public API.

static void oracle_weight_walk(const traintrack& tt, const traintracks::multigon& mm,
                               int pin, int ein, std::vector<const traintracks::edge*>& order)
{
  int p = pin, e = ein;
  mm.cycle_edges(p,e);
  do
    {
      order.push_back(mm.Edge(p,e).get());
      int pout, eout;
      traintracks::multigon* ed = mm.Edge(p,e)->target_multigon(&mm,pout,eout);
      if (ed->edges() > 1) oracle_weight_walk(tt,*ed,pout,eout,order);
      mm.cycle_edges(p,e);
    }
  while (!(p == pin && e == ein));
}

static std::vector<const traintracks::edge*> oracle_edge_order(const traintrack& tt, int mono)
{
  std::vector<const traintracks::edge*> order;
  const traintracks::multigon& m0 = tt.Multigon(mono);
  order.push_back(m0.Edge(0,0).get());
  int pmono, pemono;
  traintracks::multigon* eg = m0.Edge(0,0)->target_multigon(&m0,pmono,pemono);
  if (eg->edges() > 1) oracle_weight_walk(tt,*eg,pmono,pemono,order);
  return order;
}

// Cusps as (multigon pointer, prong, slot) in the order of the old
// recursive_find_cusp: entry slot first, then slots in cycle order,
// descending before moving on.
typedef std::vector<std::pair<const traintracks::multigon*,std::pair<int,int> > > cusp_list;

static void oracle_cusp_walk(const traintracks::multigon& mm, int pin, int ein, cusp_list& out)
{
  int p = pin, e = ein;
  do
    {
      if (e < mm.edges(p)-1) out.push_back(std::make_pair(&mm,std::make_pair(p,e)));
      int pout, eout;
      traintracks::multigon* ed = mm.Edge(p,e)->target_multigon(&mm,pout,eout);
      if (!(p == pin && e == ein) && ed->edges() > 1) oracle_cusp_walk(*ed,pout,eout,out);
      mm.cycle_edges(p,e);
    }
  while (!(p == pin && e == ein));
}

static cusp_list oracle_cusps(const traintrack& tt)
{
  cusp_list out;
  const traintracks::multigon& m0 = tt.Multigon(0);
  int pmono, pemono;
  traintracks::multigon* eg = m0.Edge(0,0)->target_multigon(&m0,pmono,pemono);
  if (eg->edges() > 1) oracle_cusp_walk(*eg,pmono,pemono,out);
  return out;
}

static void check_against_oracles(const traintrack& tt, const char* what)
{
  // Edge order from every uncusped monogon, via numbering() and weights().
  for (int mono = 0; mono < tt.multigons(); ++mono)
    {
      if (tt.Multigon(mono).edges() != 1) continue;
      const std::vector<const traintracks::edge*> order = oracle_edge_order(tt,mono);
      const ttnumbering N = traintracks::detail::coding_engine::numbering(tt,mono);
      CHECK_MSG((int)order.size() == N.nedges(), what);
      for (int e = 0; e < N.nedges(); ++e) CHECK_MSG(N.edge_ptr[e] == order[e], what);
      traintrack ttw(tt);
      traintrack::dblVec wv(N.nedges());
      for (int e = 0; e < N.nedges(); ++e) wv[e] = 10 + e;
      ttw.weights(wv.begin());
      const traintrack::dblVec back = ttw.weights(mono);
      const std::vector<const traintracks::edge*> order0 = oracle_edge_order(ttw,0);
      const std::vector<const traintracks::edge*> orderm = oracle_edge_order(ttw,mono);
      for (int e = 0; e < N.nedges(); ++e)
        {
          // weights(mono)[e] is the weight of the e-th edge of the mono walk,
          // which was set to 10 + (its position in the monogon-0 walk).
          int pos0 = -1;
          for (int k = 0; k < (int)order0.size(); ++k) if (order0[k] == orderm[e]) pos0 = k;
          CHECK_MSG(pos0 >= 0, what);
          CHECK_MSG(back[e] == 10 + pos0, what);
        }
    }

  // Cusp order via numbering() and fold_cusp_location().
  const cusp_list cl = oracle_cusps(tt);
  const ttnumbering N = tt.numbering();
  CHECK_MSG((int)cl.size() == N.ncusps() && N.ncusps() == tt.cusps(), what);
  CHECK_MSG(N.foldings() == tt.foldings(), what);
  for (int f = 0; f < tt.foldings(); ++f)
    {
      int m, p, slot;
      N.fold_cusp(f,m,p,slot);
      CHECK_MSG(&tt.Multigon(m) == cl[f/2].first, what);
      CHECK_MSG(p == cl[f/2].second.first && slot == cl[f/2].second.second, what);
      traintracks::multigon* mmc = 0;
      int pc = -1, ec = -1;
      tt.fold_cusp_location(f,mmc,pc,ec);
      CHECK_MSG(mmc == cl[f/2].first && pc == p && ec == slot, what);
    }
}

// The labels are a function of the coding.  For every vertex of an
// automaton and every legal fold from it: the fold result carries the same
// labels as the track rebuilt from its coding and as the stored target
// vertex (this is what lets ttfoldgraph compose one-step maps and the gate
// accumulator across vertices without a transport step); and at a
// cyclically symmetric vertex, the numbering started from any of the
// minimising monogons has the same labels, exactly order() of them.
static void check_labels_from_coding(const int n, const int trk,
                                     int& nfolds, int& ncyclic)
{
  jlt::vector<traintrack> ttv = traintracks::build_traintrack_list(n);
  ttauto::ttfoldgraph<traintrack> ttg(ttv[trk]);
  for (int v = 0; v < ttg.vertices(); ++v)
    {
      const traintrack& tt = ttg.traintrack(v);
      const ttnumbering N0 = tt.numbering();
      CHECK(N0.same_labels(N0));
      CHECK(N0 == N0);

      // Alternative start monogons: the same labels iff the coding from
      // there is minimal too, i.e. exactly cyclic_symmetry().order() of
      // the uncusped monogons give the same labels.
      traintrack tsym(tt);
      const int order = tsym.cyclic_symmetry().order();
      int nsame = 0;
      for (int m = 0; m < tt.multigons(); ++m)
        {
          if (tt.Multigon(m).edges() != 1) continue;
          const ttnumbering Nm = traintracks::detail::coding_engine::numbering(tt,m);
          if (Nm.same_labels(N0)) ++nsame;
          if (m != 0 && Nm.same_labels(N0)) CHECK(Nm != N0);  // bookkeeping differs
        }
      CHECK_MSG(nsame == order, "n=" << n << " trk=" << trk << " v=" << v
                << ": " << nsame << " starts share the labels, order " << order);
      if (order > 1) ++ncyclic;

      for (int f = 0; f < tt.foldings(); ++f)
        {
          traintrack t0(tt);
          if (!t0.fold(f)) continue;
          ++nfolds;
          const ttnumbering Nf = t0.numbering();
          traintrack tcode(t0.coding());
          CHECK_MSG(tcode.numbering().same_labels(Nf),
                    "fold result vs rebuilt from coding, n=" << n << " v=" << v << " f=" << f);
          int target = -1;
          for (int w = 0; w < ttg.vertices() && target < 0; ++w)
            if (ttg.traintrack(w) == t0) target = w;
          CHECK(target >= 0);
          CHECK_MSG(ttg.traintrack(target).numbering().same_labels(Nf),
                    "fold result vs stored vertex, n=" << n << " v=" << v << " f=" << f);
        }
    }
}

int main()
{
  int count = 0;
  for (int n = 3; n <= 5; ++n)
    {
      jlt::vector<traintrack> ttv = traintracks::build_traintrack_list(n);
      for (std::size_t trk = 0; trk < ttv.size(); ++trk)
        {
          check_numbering(ttv[trk],"initial track");
          check_against_oracles(ttv[trk],"initial track vs oracles");
          // Also every neighbour reached by one fold.
          for (int f = 0; f < ttv[trk].foldings(); ++f)
            {
              traintrack t2(ttv[trk]);
              if (t2.fold(f))
                {
                  check_numbering(t2,"folded track");
                  check_against_oracles(t2,"folded track vs oracles");
                  ++count;
                }
            }
        }
    }
  CHECK(count > 0);

  int nfolds = 0, ncyclic = 0;
  check_labels_from_coding(4,0,nfolds,ncyclic);
  check_labels_from_coding(5,0,nfolds,ncyclic);
  check_labels_from_coding(6,2,nfolds,ncyclic);
  check_labels_from_coding(6,4,nfolds,ncyclic);
  CHECK(nfolds > 0);
  CHECK_MSG(ncyclic > 0, "no cyclically symmetric vertex visited");

  std::cout << "test_numbering: OK (" << count << " folded tracks checked; "
            << nfolds << " fold results carry the labels of their coding, "
            << ncyclic << " cyclically symmetric vertices)\n";
  return 0;
}
