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

#include <cassert>
#include <algorithm>
#include <iostream>
#include <list>
#include <map>
#include <vector>

#include <jlt/freeauto.hpp>
#include <jlt/mathmatrix.hpp>
#include <jlt/vector.hpp>

#include "traintracks/build.hpp"
#include "traintracks/map.hpp"
#include "traintracks/traintrack.hpp"
#include "ttauto/folding_path.hpp"
#include "ttauto/ttfoldgraph.hpp"


// Derivative map on oriented generators.
//
// For each signed generator a, D(a) is the first oriented generator appearing
// in the reduced word AM(a). This captures the one-step "turn" dynamics on
// directions, not the full edge-path image.
//
// Important: D is not a free-group automorphism in general. In particular,
// D(-a) need not equal -D(a), because D only keeps first-letter data.
void print_derivative_map(const jlt::freeauto<int>& AM,
                          const int ngen,
                          std::ostream& strm = std::cout)
{
  strm << "\nDerivative map D on signed main generators:\n";
  for (int a = 1; a <= ngen; ++a)
    {
      jlt::freeword<int> img_pos = AM.get_action(a);
      jlt::freeword<int> img_neg = AM.get_action(-a);
      assert(!img_pos.empty());
      assert(!img_neg.empty());

      strm << "  D(" << a << ") = " << *img_pos.begin();
      strm << ", D(" << -a << ") = " << *img_neg.begin() << "\n";
    }
}

int multigon_index_from_ptr(const traintracks::traintrack& tt,
                            const traintracks::multigon* mm)
{
  for (int m = 0; m < tt.multigons(); ++m)
    {
      if (&tt.Multigon(m) == mm) return m;
    }
  std::cerr << "Could not resolve multigon pointer to index.\n";
  std::exit(1);
}

struct endpoint
{
  int m;
  int p;
};

typedef std::pair<int,int> vertex_key;

std::vector<std::pair<endpoint,endpoint> >
main_edge_end_vertices(traintracks::traintrack tt)
{
  const int n = tt.edges();
  std::vector<std::pair<endpoint,endpoint> > ends(
      n,std::make_pair(endpoint{-1,-1},endpoint{-1,-1}));
  const traintracks::traintrack& ctt = tt;

  for (int ei = 0; ei < n; ++ei)
    {
      traintracks::traintrack::dblVec wv(n,0.0);
      wv[ei] = 1.0;
      tt.weights(wv.begin());

      bool found = false;
      for (int m = 0; m < ctt.multigons() && !found; ++m)
        {
          for (int p = 0; p < ctt.Multigon(m).prongs() && !found; ++p)
            {
              for (int e = 0; e < ctt.Multigon(m).edges(p) && !found; ++e)
                {
                  auto ep = ctt.Multigon(m).Edge(p,e);
                  if (ep->weight() != 1.0) continue;
                  const traintracks::multigon* mm0 = ep->ending_multigon(0);
                  const traintracks::multigon* mm1 = ep->ending_multigon(1);
                  ends[ei] = std::make_pair(
                    endpoint{multigon_index_from_ptr(tt,mm0),ep->ending_prong(0)},
                    endpoint{multigon_index_from_ptr(tt,mm1),ep->ending_prong(1)});
                  found = true;
                }
            }
        }
      assert(found);
    }

  return ends;
}

std::vector<std::pair<endpoint,endpoint> >
infinitesimal_end_vertices(const traintracks::traintrack& tt)
{
  const int ninf = tt.total_prongs();
  std::vector<std::pair<endpoint,endpoint> > ends(
      ninf,std::make_pair(endpoint{-1,-1},endpoint{-1,-1}));

  int ix = 0;
  for (int m = 0; m < tt.multigons(); ++m)
    {
      const int k = tt.Multigon(m).prongs();
      for (int p = 0; p < k; ++p)
        {
          // Canonical infinitesimal orientation: from prong p to next prong
          // (clockwise cyclic order on the multigon boundary).
          ends[ix] = std::make_pair(endpoint{m,p},endpoint{m,(p+1)%k});
          ++ix;
        }
    }
  assert(ix == ninf);

  return ends;
}

void print_vertex_directed_generators(const traintracks::traintrack& tt,
                                      std::ostream& strm = std::cout)
{
  const int nmain = tt.edges();
  const int ninf = tt.total_prongs();
  const traintracks::ttmap_labeler labels(nmain,ninf);

  const std::vector<std::pair<endpoint,endpoint> > medges =
    main_edge_end_vertices(tt);
  const std::vector<std::pair<endpoint,endpoint> > iedges =
    infinitesimal_end_vertices(tt);

  std::map<std::pair<int,int>,std::vector<int> > by_start;

  // Main directed generators.
  for (int e = 0; e < nmain; ++e)
    {
      const int g = labels.main_gen(e);
      by_start[std::make_pair(medges[e].first.m,medges[e].first.p)].push_back(g);
      by_start[std::make_pair(medges[e].second.m,medges[e].second.p)].push_back(-g);
    }

  // Infinitesimal directed generators.
  for (int i = 0; i < ninf; ++i)
    {
      const int g = labels.infinitesimal_gen(i);
      by_start[std::make_pair(iedges[i].first.m,iedges[i].first.p)].push_back(g);
      by_start[std::make_pair(iedges[i].second.m,iedges[i].second.p)].push_back(-g);
    }

  strm << "\nDirected generators by vertex (multigon,prong):\n";
  for (auto& kv : by_start)
    {
      std::vector<int>& dirs = kv.second;
      std::sort(dirs.begin(),dirs.end());
      strm << "  vertex (" << kv.first.first+1 << "," << kv.first.second+1
           << "): ";
      for (int i = 0; i < (int)dirs.size(); ++i)
        {
          strm << dirs[i];
          if (i+1 < (int)dirs.size()) strm << " ";
        }
      strm << "\n";
    }
}

vertex_key start_vertex_of_generator(
    const int g,
    const int nmain,
    const std::vector<std::pair<endpoint,endpoint> >& medges,
    const std::vector<std::pair<endpoint,endpoint> >& iedges)
{
  const int gi = std::abs(g);
  if (gi <= nmain)
    {
      const int ei = gi-1;
      if (g > 0)
        return std::make_pair(medges[ei].first.m,medges[ei].first.p);
      return std::make_pair(medges[ei].second.m,medges[ei].second.p);
    }

  const int ii = gi - nmain - 1;
  assert(ii >= 0 && ii < (int)iedges.size());
  if (g > 0)
    return std::make_pair(iedges[ii].first.m,iedges[ii].first.p);
  return std::make_pair(iedges[ii].second.m,iedges[ii].second.p);
}

struct gate_partition
{
  std::map<vertex_key,std::vector<std::vector<int> > > gates;
  std::map<vertex_key,std::map<int,int> > gate_of_dir;
};

gate_partition compute_gate_partition(const traintracks::traintrack& tt,
                                      const jlt::freeauto<int>& AM)
{
  gate_partition gp;

  const int nmain = tt.edges();
  const int ninf = tt.total_prongs();
  const traintracks::ttmap_labeler labels(nmain,ninf);
  const int ngen = labels.num_generators();

  const std::vector<std::pair<endpoint,endpoint> > medges =
    main_edge_end_vertices(tt);
  const std::vector<std::pair<endpoint,endpoint> > iedges =
    infinitesimal_end_vertices(tt);

  // Build derivative map D on all signed generators +/-1..+/-ngen.
  std::map<int,int> D;
  for (int a = 1; a <= ngen; ++a)
    {
      int dpos = *AM.get_action(a).begin();
      int dneg = *AM.get_action(-a).begin();
      assert(std::abs(dpos) <= ngen && dpos != 0);
      assert(std::abs(dneg) <= ngen && dneg != 0);
      D[a] = dpos;
      D[-a] = dneg;
    }

  auto eventually_meet = [&](const int a, const int b) -> bool {
    int x = a, y = b;
    for (int k = 1; k <= 2*ngen+1; ++k)
      {
        x = D[x];
        y = D[y];
        if (x == y) return true;
      }
    return false;
  };

  std::map<vertex_key,std::vector<int> > by_start;
  for (int g = 1; g <= ngen; ++g)
    {
      by_start[start_vertex_of_generator(g,nmain,medges,iedges)].push_back(g);
      by_start[start_vertex_of_generator(-g,nmain,medges,iedges)].push_back(-g);
    }

  for (auto& kv : by_start)
    {
      const vertex_key v = kv.first;
      std::vector<int>& dirs = kv.second;
      std::sort(dirs.begin(),dirs.end());

      const int nd = dirs.size();
      std::vector<int> parent(nd);
      for (int i = 0; i < nd; ++i) parent[i] = i;

      auto find = [&](int x) {
        while (parent[x] != x)
          {
            parent[x] = parent[parent[x]];
            x = parent[x];
          }
        return x;
      };

      auto unite = [&](int a, int b) {
        int ra = find(a), rb = find(b);
        if (ra != rb) parent[rb] = ra;
      };

      for (int i = 0; i < nd; ++i)
        {
          for (int j = i+1; j < nd; ++j)
            {
              if (eventually_meet(dirs[i],dirs[j])) unite(i,j);
            }
        }

      std::map<int,std::vector<int> > classes;
      for (int i = 0; i < nd; ++i)
        {
          classes[find(i)].push_back(dirs[i]);
        }

      int gid = 0;
      for (auto& gk : classes)
        {
          std::vector<int>& gate = gk.second;
          std::sort(gate.begin(),gate.end());
          gp.gates[v].push_back(gate);
          for (int x : gate) gp.gate_of_dir[v][x] = gid;
          ++gid;
        }
    }

  return gp;
}

void print_gate_candidates(const traintracks::traintrack& tt,
                           const jlt::freeauto<int>& AM,
                           std::ostream& strm = std::cout)
{
  gate_partition gp = compute_gate_partition(tt,AM);

  strm << "\nGates by vertex (multigon,prong) on all generators "
       << "(main+infinitesimal):\n";
  for (auto& kv : gp.gates)
    {
      const int m = kv.first.first;
      const int p = kv.first.second;
      std::vector<std::vector<int> >& gates = kv.second;

      strm << "  vertex (" << m+1 << "," << p+1 << "): ";
      bool first_gate = true;
      for (auto& gate : gates)
        {
          if (!first_gate) strm << ", ";
          strm << "{";
          for (int i = 0; i < (int)gate.size(); ++i)
            {
              strm << gate[i];
              if (i+1 < (int)gate.size()) strm << " ";
            }
          strm << "}";
          first_gate = false;
        }
      strm << "\n";
    }
}


int main()
{
  using traintracks::traintrack;
  using traintracks::build_traintrack_list;
  using ttauto::ttfoldgraph;
  using ttauto::folding_path;

  typedef ttfoldgraph<traintrack> ttgraph;

  const int n = 6;
  jlt::vector<traintrack> ttv = build_traintrack_list(n);

  // Issue #2 repro data from devel/iss002/n=6_5_bad_tt_data.m
  // Vertex cycle listed there is 1-based.
  const int bad_cycle_1based[] = {29, 46, 43, 71, 88, 85, 29};
  const int cycle_len = sizeof(bad_cycle_1based)/sizeof(bad_cycle_1based[0]);
  // Hardwired branch sequence for the bad example (0-based branch ids at
  // each step in the cycle above).
  const int bad_branches[] = {1, 0, 3, 2, 1, 2};

  // Hardwire the known bad graph: n=6, traintrack index 4, first pruned
  // subgraph (90 vertices) from devel/iss002/n=6_5_bad_tt_data.m.
  const int trk = 4;
  const int sgidx = 0;
  ttgraph full(ttv[trk]);
  std::list<ttgraph> sgs = ttauto::subgraphs(full);
  ttauto::prune_multihumps(sgs);
  assert((int)sgs.size() > sgidx);
  auto it = sgs.begin();
  std::advance(it,sgidx);
  const ttgraph& ttg = *it;

  std::cout << "Issue #2 bad-path reproducer\n";
  std::cout << "n=" << n << ", track index=" << trk
            << ", subgraph index=" << sgidx
            << ", vertices=" << ttg.vertices() << "\n";
  assert(ttg.vertices() == 90);

  int bad_cycle_0based[cycle_len];
  for (int i = 0; i < cycle_len; ++i)
    {
      bad_cycle_0based[i] = bad_cycle_1based[i]-1;
      assert(bad_cycle_0based[i] >= 0);
      assert(bad_cycle_0based[i] < ttg.vertices());
    }

  const int vstart = bad_cycle_0based[0];

  std::cout << "vertex cycle (1-based): ";
  for (int i = 0; i < cycle_len; ++i)
    {
      std::cout << bad_cycle_1based[i];
      if (i+1 < cycle_len) std::cout << " -> ";
    }
  std::cout << "\n";

  std::cout << "branch sequence (0-based): ";
  for (int i = 0; i < cycle_len-1; ++i)
    {
      const int v = bad_cycle_0based[i];
      const int vnext = bad_cycle_0based[i+1];

      const int b = bad_branches[i];
      assert(b >= 0 && b < ttg.foldings(v));
      assert(ttg.target_vertex(v,b) == vnext);
      if (i > 0) std::cout << ", ";
      std::cout << b;

      if (i == 0)
        {
          // Non-essential sanity print to make it obvious this is the
          // hardwired known-bad realization.
          std::cout << " [hardwired]";
        }
    }
  std::cout << "\n";

  folding_path<traintrack> p(ttg,vstart);
  for (int i = 0; i < cycle_len-1; ++i)
    {
      p.push_back(bad_branches[i]);
    }

  jlt::freeauto<int> AM = p.traintrack_map();
  jlt::mathmatrix<int> TM = p.transition_matrix();
  jlt::mathmatrix<int> TMfromAM =
    traintracks::transition_matrix_from_map(ttg.traintrack(vstart),AM);
  assert(TM == TMfromAM);

  std::cout << "\nComposed train-track map:\n";
  std::cout << AM;
  std::cout << "\nComposed transition matrix (path):\n";
  TM.printMatrixForm(std::cout);

  const bool primitive = TM.is_primitive();
  std::cout << "\nPrimitive transition matrix: "
            << (primitive ? "yes" : "no") << "\n";
  assert(primitive);

  std::cout << "\nTransition matrix from composed map:\n";
  TMfromAM.printMatrixForm(std::cout);

  std::cout << "\nOK: matrix(path) == matrix(from map).\n";

  // The transition matrix acts on main edges only, so for D we print only
  // signed main generators +/-1..+/-edges().
  print_derivative_map(AM,ttg.edges());
  print_vertex_directed_generators(ttg.traintrack(vstart));
  print_gate_candidates(ttg.traintrack(vstart),AM);

  return 0;
}
