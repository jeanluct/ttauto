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

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <vector>
#include "traintracks/graph_iso.hpp"
#include "traintracks/traintrack.hpp"

namespace traintracks {
namespace graph_iso {

// Isotopy approach used here:
// - Convert each train track to a prong-incidence multigraph where endpoints are
//   (multigon,prong) slots and branch multiplicity is preserved.
// - Solve orientation-preserving isomorphism by backtracking over multigon and
//   prong correspondences.
// - Preserve multigon labels and puncture flags in candidate matching.
//
// This is not a canonical-coding sequence algorithm; it is a direct structural
// matcher. It is in the same family as Ullmann/VF2-style recursive pruning,
// specialized to train-track prong structure and cyclic order constraints.

namespace {

struct endpoint
{
  int m;
  int p;

  bool operator==(const endpoint& o) const
  {
    return (m == o.m && p == o.p);
  }
};

using endpoint_pair = std::pair<endpoint,endpoint>;

// Order endpoints deterministically in maps.
inline bool endpoint_less(const endpoint& a, const endpoint& b)
{
  if (a.m != b.m) return a.m < b.m;
  return a.p < b.p;
}

// Order endpoint pairs deterministically in maps.
inline bool endpoint_pair_less(const endpoint_pair& a,
                               const endpoint_pair& b)
{
  if (!(a.first == b.first)) return endpoint_less(a.first,b.first);
  return endpoint_less(a.second,b.second);
}

// Normalize endpoint pair ordering.
void normalize_endpoint_pair(endpoint_pair& ep)
{
  if (endpoint_less(ep.second,ep.first)) std::swap(ep.first,ep.second);
}

using endpoint_count_map = std::map<endpoint_pair,int,
                                    bool(*)(const endpoint_pair&,const endpoint_pair&)>;

// Collect branch multiplicities between prong endpoints.
endpoint_count_map collect_branch_endpoints(const traintrack& tt)
{
  endpoint_count_map counts(endpoint_pair_less);

  std::map<const edge*,endpoint> first;

  for (int m = 0; m < tt.multigons(); ++m)
    {
      const multigon& mm = tt.Multigon(m);
      for (int p = 0; p < mm.prongs(); ++p)
        {
          for (int e = 0; e < mm.edges(p); ++e)
            {
#ifdef TRAINTRACKS_NO_SHARED_PTR
              const edge* eg = mm.Edge(p,e);
#else
              const edge* eg = mm.Edge(p,e).get();
#endif
              auto it = first.find(eg);
              if (it == first.end())
                {
                  first[eg] = endpoint{m,p};
                }
              else
                {
                  endpoint_pair ep(it->second,endpoint{m,p});
                  normalize_endpoint_pair(ep);
                  ++counts[ep];
                  first.erase(it);
                }
            }
        }
    }

  if (!first.empty())
    {
      std::cerr << "Unpaired branch endpoint in graph_iso::collect_branch_endpoints.\n";
      std::exit(1);
    }

  return counts;
}

// Generate orientation-preserving prong maps for a k-prong multigon.
//   For k>=3 this is exactly the k cyclic rotations.
std::vector<std::vector<int> > all_oriented_cycles(const int k)
{
  std::vector<std::vector<int> > out;
  if (k <= 0) return out;

  out.reserve(k);
  for (int shift = 0; shift < k; ++shift)
    {
      std::vector<int> rot(k);
      for (int i = 0; i < k; ++i) rot[i] = (i + shift) % k;
      out.push_back(rot);
    }

  return out;
}

// generate candidate prong maps for a pair of multigons.
std::vector<std::vector<int> > mapping_options_for_multigon(const multigon& a,
                                                             const multigon& b)
{
  const int k = a.prongs();
  if (k != b.prongs()) return {};

  if (k == 1)
    {
      return {std::vector<int>{0}};
    }
  if (k == 2)
    {
      return {std::vector<int>{0,1},std::vector<int>{1,0}};
    }
  return all_oriented_cycles(k);
}

struct iso_search
{
  const traintrack& a;
  const traintrack& b;
  std::vector<int> mb;
  std::vector<int> used_b;
  std::vector<std::vector<int> > prong_map;
  endpoint_count_map cnt_a;
  endpoint_count_map cnt_b;

  iso_search(const traintrack& aa, const traintrack& bb)
    : a(aa), b(bb), mb(aa.multigons(),-1), used_b(bb.multigons(),0),
      prong_map(aa.multigons()), cnt_a(collect_branch_endpoints(aa)),
      cnt_b(collect_branch_endpoints(bb)) {}

  // Prune search using already-fixed multigon/prong correspondences.
  bool check_partial_edges() const
  {
    endpoint_count_map mapped(endpoint_pair_less);

    for (const auto& kv : cnt_a)
      {
        const endpoint_pair& epa = kv.first;
        const int c = kv.second;

        const int m1 = mb[epa.first.m];
        const int m2 = mb[epa.second.m];
        if (m1 < 0 || m2 < 0) continue;

        const int p1 = prong_map[epa.first.m][epa.first.p];
        const int p2 = prong_map[epa.second.m][epa.second.p];

        endpoint_pair epb(endpoint{m1,p1},endpoint{m2,p2});
        normalize_endpoint_pair(epb);
        mapped[epb] += c;
      }

    for (const auto& kv : mapped)
      {
        auto it = cnt_b.find(kv.first);
        if (it == cnt_b.end()) return false;
        if (it->second < kv.second) return false;
      }

    return true;
  }

  // Check whether all multigons have been assigned.
  bool done() const
  {
    for (int m = 0; m < a.multigons(); ++m)
      {
        if (mb[m] < 0) return false;
      }
    return true;
  }

  // Verify full edge-multiplicity agreement for a complete mapping.
  bool check_full_edges() const
  {
    endpoint_count_map mapped(endpoint_pair_less);
    for (const auto& kv : cnt_a)
      {
        const endpoint_pair& epa = kv.first;
        const int c = kv.second;
        const int m1 = mb[epa.first.m];
        const int m2 = mb[epa.second.m];
        const int p1 = prong_map[epa.first.m][epa.first.p];
        const int p2 = prong_map[epa.second.m][epa.second.p];
        endpoint_pair epb(endpoint{m1,p1},endpoint{m2,p2});
        normalize_endpoint_pair(epb);
        mapped[epb] += c;
      }
    return (mapped == cnt_b);
  }

  // Choose next unassigned multigon (largest prong count first).
  int choose_next_m() const
  {
    int best = -1;
    for (int m = 0; m < a.multigons(); ++m)
      {
        if (mb[m] >= 0) continue;
        if (best < 0) { best = m; continue; }
        if (a.Multigon(m).prongs() > a.Multigon(best).prongs()) best = m;
      }
    return best;
  }

  // Recursively search for an orientation-preserving isomorphism.
  bool rec()
  {
    if (done()) return check_full_edges();

    const int m = choose_next_m();
    const multigon& ma = a.Multigon(m);

    for (int n = 0; n < b.multigons(); ++n)
      {
        if (used_b[n]) continue;
        const multigon& mbb = b.Multigon(n);
        if (ma.prongs() != mbb.prongs()) continue;
        if (ma.label() != mbb.label()) continue;
        if (ma.punctured() != mbb.punctured()) continue;

        std::vector<std::vector<int> > opts = mapping_options_for_multigon(ma,mbb);
        for (const auto& pmap : opts)
          {
            mb[m] = n;
            used_b[n] = 1;
            prong_map[m] = pmap;

            if (check_partial_edges() && rec()) return true;

            mb[m] = -1;
            used_b[n] = 0;
            prong_map[m].clear();
          }
      }

    return false;
  }
};

struct witness_search
{
  struct mg_sig
  {
    int prongs;
    int punct;
    int label;

    bool operator==(const mg_sig& o) const
    {
      return (prongs == o.prongs && punct == o.punct && label == o.label);
    }

    bool operator<(const mg_sig& o) const
    {
      if (prongs != o.prongs) return prongs < o.prongs;
      if (punct != o.punct) return punct < o.punct;
      return label < o.label;
    }
  };

  const traintrack& t;
  std::vector<mg_sig> sig_by_m;
  std::vector<mg_sig> sorted_sigs;
  std::vector<int> order;
  std::vector<int> used;
  std::vector<int> shift;

  std::vector<int> best_cert;
  std::vector<int> best_order;
  std::vector<int> best_shift;
  bool have_best;

  witness_search(const traintrack& tt)
    : t(tt), sig_by_m(tt.multigons()), sorted_sigs(tt.multigons()),
      order(tt.multigons(),-1), used(tt.multigons(),0),
      shift(tt.multigons(),0), have_best(false) {}

  void init_signatures()
  {
    for (int m = 0; m < t.multigons(); ++m)
      {
        sig_by_m[m] = mg_sig{t.Multigon(m).prongs(),
                             t.Multigon(m).punctured() ? 1 : 0,
                             t.Multigon(m).label()};
        sorted_sigs[m] = sig_by_m[m];
      }
    std::sort(sorted_sigs.begin(),sorted_sigs.end());
  }

  std::vector<int> build_certificate() const
  {
    std::vector<int> cert;
    cert.reserve(4*t.multigons() + 8*t.edges());

    std::vector<mg_sig> sig_at_slot(t.multigons());
    for (int m = 0; m < t.multigons(); ++m)
      {
        sig_at_slot[order[m]] = sig_by_m[m];
      }

    for (int s = 0; s < t.multigons(); ++s)
      {
        cert.push_back(sig_at_slot[s].prongs);
        cert.push_back(sig_at_slot[s].punct);
        cert.push_back(sig_at_slot[s].label);
      }

    endpoint_count_map cnt = collect_branch_endpoints(t);
    std::vector<std::vector<int> > edgesig;
    edgesig.reserve(cnt.size());
    for (const auto& kv : cnt)
      {
        int m1 = order[kv.first.first.m];
        int p1 = (kv.first.first.p + shift[kv.first.first.m]) %
                 t.Multigon(kv.first.first.m).prongs();
        int m2 = order[kv.first.second.m];
        int p2 = (kv.first.second.p + shift[kv.first.second.m]) %
                 t.Multigon(kv.first.second.m).prongs();
        if (m2 < m1 || (m2 == m1 && p2 < p1))
          {
            std::swap(m1,m2);
            std::swap(p1,p2);
          }
        edgesig.push_back({m1,p1,m2,p2,kv.second});
      }
    std::sort(edgesig.begin(),edgesig.end());
    for (const auto& e : edgesig)
      {
        cert.insert(cert.end(),e.begin(),e.end());
      }

    return cert;
  }

  void maybe_record_best()
  {
    std::vector<int> cert = build_certificate();
    std::vector<int> ord = order;

    if (!have_best || cert < best_cert || (cert == best_cert && ord < best_order))
      {
        best_cert = cert;
        best_order = ord;
        best_shift = shift;
        have_best = true;
      }
  }

  void rec(const int slot)
  {
    if (slot >= t.multigons())
      {
        maybe_record_best();
        return;
      }

    const mg_sig& req = sorted_sigs[slot];

    for (int m = 0; m < t.multigons(); ++m)
      {
        if (used[m]) continue;
        if (!(sig_by_m[m] == req)) continue;

        used[m] = 1;
        order[m] = slot;

        const multigon& mm = t.Multigon(m);
        const int k = mm.prongs();
        if (k <= 0)
          {
            used[m] = 0;
            order[m] = -1;
            continue;
          }

        for (int s = 0; s < k; ++s)
          {
            shift[m] = s;
            rec(slot+1);
          }

        used[m] = 0;
        order[m] = -1;
      }
  }
};

} // namespace

// Test orientation-preserving train-track isotopy.
bool is_isotopic_oriented(const traintrack& lhs, const traintrack& rhs)
{
  if (lhs.multigons() != rhs.multigons()) return false;
  if (lhs.edges() != rhs.edges()) return false;

  iso_search s(lhs,rhs);
  return s.rec();
}

graph_iso::canonical_witness canonical_witness_oriented(const traintrack& tt)
{
  witness_search ws(tt);
  ws.init_signatures();
  ws.rec(0);

  graph_iso::canonical_witness w;
  w.valid = ws.have_best;
  if (!w.valid) return w;

  w.multigon_rank = ws.best_order;
  w.prong_shift = ws.best_shift;
  w.multigon_order.assign(tt.multigons(),-1);
  for (int m = 0; m < tt.multigons(); ++m)
    {
      const int r = w.multigon_rank[m];
      if (r >= 0 && r < tt.multigons()) w.multigon_order[r] = m;
    }
  return w;
}

} // namespace graph_iso
} // namespace traintracks
