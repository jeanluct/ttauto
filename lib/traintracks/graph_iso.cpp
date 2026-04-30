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
#include <set>
#include <vector>
#include "traintracks/graph_iso.hpp"
#include "traintracks/traintrack.hpp"

namespace traintracks {
namespace graph_iso {

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

inline bool endpoint_less(const endpoint& a, const endpoint& b)
{
  if (a.m != b.m) return a.m < b.m;
  return a.p < b.p;
}

inline bool endpoint_equal(const endpoint& a, const endpoint& b)
{
  return (a.m == b.m && a.p == b.p);
}

inline bool endpoint_pair_less(const endpoint_pair& a,
                               const endpoint_pair& b)
{
  if (!endpoint_equal(a.first,b.first)) return endpoint_less(a.first,b.first);
  return endpoint_less(a.second,b.second);
}

void normalize_endpoint_pair(endpoint_pair& ep)
{
  if (endpoint_less(ep.second,ep.first)) std::swap(ep.first,ep.second);
}

using endpoint_count_map = std::map<endpoint_pair,int,
                                    bool(*)(const endpoint_pair&,const endpoint_pair&)>;

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

std::vector<std::vector<int> > all_oriented_cycles(const int k)
{
  std::vector<std::vector<int> > out;
  if (k <= 0) return out;

  std::vector<int> perm(k);
  for (int i = 0; i < k; ++i) perm[i] = i;

  do
    {
      bool ok = true;
      for (int i = 0; i < k; ++i)
        {
          int a = perm[i];
          int b = perm[(i+1)%k];
          if (((a+1)%k) != b) { ok = false; break; }
        }
      if (ok) out.push_back(perm);
    }
  while (std::next_permutation(perm.begin(),perm.end()));

  return out;
}

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

  bool done() const
  {
    for (int m = 0; m < a.multigons(); ++m)
      {
        if (mb[m] < 0) return false;
      }
    return true;
  }

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

} // namespace

bool is_isotopic_oriented(const traintrack& lhs, const traintrack& rhs)
{
  if (lhs.multigons() != rhs.multigons()) return false;
  if (lhs.edges() != rhs.edges()) return false;

  iso_search s(lhs,rhs);
  return s.rec();
}

} // namespace graph_iso
} // namespace traintracks
