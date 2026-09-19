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

// Bestvina-Handel gate test.  See include/traintracks/gates.hpp.

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <vector>

#include "traintracks/gates.hpp"

namespace traintracks {

namespace {

void fail(const char* msg)
{
  std::cerr << msg << " in traintracks::gates.\n";
  std::exit(1);
}

turn make_turn(int a, int b)
{
  if (a > b) std::swap(a,b);
  return turn(a,b);
}

// A letter is a direction if it is a main letter or a peripheral side.
bool is_direction(const ttnumbering& N, const int letter)
{
  if (N.is_main(letter)) return true;
  if (!N.is_side(letter)) return false;
  return N.prong[N.side_of(letter)].punctured;
}

// A letter is "real" (an edge of the Bestvina-Handel graph G) if it is a
// direction; sides of unpunctured multigons are infinitesimal.
inline bool is_real(const ttnumbering& N, const int letter)
{
  return is_direction(N,letter);
}

struct dsu
{
  std::vector<int> parent;
  explicit dsu(int n) : parent(n) { for (int i = 0; i < n; ++i) parent[i] = i; }
  int find(int x) { while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; }
  void unite(int a, int b) { a = find(a); b = find(b); if (a != b) parent[b] = a; }
};

} // anonymous namespace


int bh_vertex_of(const ttnumbering& N, const int letter)
{
  const int q = N.tail_of(letter);
  if (N.prong[q].punctured) return q;
  return -(N.prong[q].multigon + 1);
}


fold_derivative::fold_derivative(const jlt::freeauto<int>& AM,
                                 const ttnumbering& target)
  : ngen(AM.numgens()), D(2*AM.numgens()+1,0)
{
  if (ngen != target.nedges() + target.nprongs())
    fail("Map and numbering have different generator counts");

  for (int a = 1; a <= ngen; ++a)
    {
      const jlt::freeword<int> wp = AM.get_action(a);
      const jlt::freeword<int> wn = AM.get_action(-a);
      if (wp.empty() || wn.empty()) fail("Empty image word");
      D[a + ngen] = *wp.begin();
      D[-a + ngen] = *wn.begin();
    }

  // Turns taken by the images of main edges: consecutive real letters.
  for (int e = 0; e < target.nedges(); ++e)
    {
      const jlt::freeword<int> w = AM.get_action(target.main_letter(e));
      int prev = 0;        // previous letter of any kind
      int prev_real = 0;   // previous main or peripheral letter
      for (auto x : w)
        {
          if (prev != 0 && target.head_of(prev) != target.tail_of(x))
            fail("Image word is not a continuous path");
          prev = x;
          if (!is_real(target,x)) continue;
          // Arrive along prev_real (direction -prev_real at the vertex),
          // leave along x; infinitesimal sides in between are dropped.
          if (prev_real != 0) turns.insert(make_turn(-prev_real,x));
          prev_real = x;
        }
    }
}


gate_accumulator::gate_accumulator(const int ngen)
  : ng(ngen), len(0), Dacc(2*ngen+1,0)
{
  for (int a = -ngen; a <= ngen; ++a) Dacc[a + ngen] = a;
}

void gate_accumulator::push_back(const fold_derivative& d)
{
  if (ng == 0)
    {
      ng = d.ngen;
      Dacc.assign(2*ng+1,0);
      for (int a = -ng; a <= ng; ++a) Dacc[a + ng] = a;
    }
  if (d.ngen != ng) fail("Generator count changed along the path");

  // T <- T(d) u D(d)(T);  D <- D(d) o D.
  std::set<turn> T(d.turns);
  for (const turn& t : Tacc)
    T.insert(make_turn(d.apply(t.first),d.apply(t.second)));
  Tacc.swap(T);

  for (int a = -ng; a <= ng; ++a)
    if (a != 0) Dacc[a + ng] = d.apply(Dacc[a + ng]);
  ++len;
}


namespace {

// Shared analysis given the composed derivative and realised turns.
gate_analysis analyse_impl(const ttnumbering& N,
                           const std::vector<int>& D,
                           std::set<turn> T)
{
  const int ngen = N.nedges() + N.nprongs();
  if ((int)D.size() != 2*ngen + 1) fail("Derivative has the wrong size");
  auto apply = [&](int a) { return D[a + ngen]; };

  // D must send directions to directions, preserving peripherality.
  for (int a = -ngen; a <= ngen; ++a)
    {
      if (a == 0 || !is_direction(N,a)) continue;
      const int b = apply(a);
      if (!is_direction(N,b)) fail("Derivative of a direction is not a direction");
      if (N.is_side(a) != N.is_side(b)) fail("Derivative mixes main and peripheral directions");
    }

  // Close the turns under D x D.
  for (bool grew = true; grew;)
    {
      grew = false;
      std::vector<turn> cur(T.begin(),T.end());
      for (const turn& t : cur)
        {
          const turn u = make_turn(apply(t.first),apply(t.second));
          if (u.first == u.second) fail("Realised turn collapses under the derivative");
          if (T.insert(u).second) grew = true;
        }
    }

  // Directions at each Bestvina-Handel vertex, in cyclic order.
  std::map<int,std::vector<int> > dirs;   // vertex id -> directions
  std::map<int,std::pair<int,int> > vertex_loc;  // id -> (multigon, prong or -1)
  for (int m = 0; m < (int)N.prong_number.size(); ++m)
    {
      for (int p = 0; p < (int)N.prong_number[m].size(); ++p)
        {
          const int q = N.prong_number[m][p];
          const std::vector<int> d = N.directions_at_prong(q);
          if (N.prong[q].punctured)
            {
              dirs[q] = d;
              vertex_loc[q] = std::make_pair(m,q);
            }
          else
            {
              std::vector<int>& v = dirs[-(m+1)];
              v.insert(v.end(),d.begin(),d.end());
              vertex_loc[-(m+1)] = std::make_pair(m,-1);
            }
        }
    }

  // Eventual coincidence under D: D^K with K = (#directions)^2 is enough,
  // since two orbits that ever meet stay together and if they have not met
  // after that many steps they never will.
  int ndirs = 0;
  for (auto& kv : dirs) ndirs += kv.second.size();
  std::vector<int> DK(2*ngen+1);
  for (int a = -ngen; a <= ngen; ++a) DK[a + ngen] = a;
  for (int k = 0; k < ndirs*ndirs; ++k)
    for (int a = -ngen; a <= ngen; ++a)
      if (a != 0) DK[a + ngen] = apply(DK[a + ngen]);

  gate_analysis result;
  result.connected = true;
  result.shapes_ok = true;

  // Which vertex each direction belongs to, for checking turns.
  std::map<int,int> vertex_of;
  for (auto& kv : dirs) for (int a : kv.second) vertex_of[a] = kv.first;

  for (auto& kv : dirs)
    {
      const int vid = kv.first;
      const std::vector<int>& d = kv.second;
      gate_vertex_report rep;
      rep.multigon = vertex_loc[vid].first;
      rep.prong = vertex_loc[vid].second;
      rep.directions = d;
      {
        std::ostringstream os;
        const int q0 = (rep.prong >= 0 ? rep.prong : N.prong_number[rep.multigon][0]);
        if (rep.prong >= 0)
          os << "prong " << rep.prong << " of punctured "
             << N.prong[q0].nprongs << "-gon " << rep.multigon;
        else
          os << "unpunctured " << N.prong[q0].nprongs << "-gon " << rep.multigon;
        rep.name = os.str();
      }

      // Gates: group consecutive directions with equal D^K image.  (Gates
      // are consecutive in the cyclic order, but two separated blocks with
      // the same image would be a violation; use a map to be safe.)
      std::map<int,int> gate_of_image;
      std::vector<int> gate_of(d.size());
      for (std::size_t i = 0; i < d.size(); ++i)
        {
          const int img = DK[d[i] + ngen];
          auto it = gate_of_image.find(img);
          if (it == gate_of_image.end())
            {
              it = gate_of_image.insert(std::make_pair(img,(int)rep.gates.size())).first;
              rep.gates.push_back(std::vector<int>());
            }
          gate_of[i] = it->second;
          rep.gates[it->second].push_back(d[i]);
        }
      std::map<int,int> gate_of_dir;
      for (std::size_t i = 0; i < d.size(); ++i) gate_of_dir[d[i]] = gate_of[i];

      // Joins from realised turns at this vertex.
      const int k = rep.gates.size();
      dsu comp(k);
      std::set<std::pair<int,int> > joined_gates;
      for (const turn& t : T)
        {
          auto ia = vertex_of.find(t.first), ib = vertex_of.find(t.second);
          if (ia == vertex_of.end() || ib == vertex_of.end())
            fail("Realised turn involves a non-direction");
          if (ia->second != vid) continue;
          if (ib->second != vid) fail("Realised turn joins different vertices");
          const int ga = gate_of_dir[t.first], gb = gate_of_dir[t.second];
          if (ga == gb) fail("Realised turn inside a gate: map is not efficient");
          rep.joins.push_back(t);
          comp.unite(ga,gb);
          joined_gates.insert(std::make_pair(std::min(ga,gb),std::max(ga,gb)));
        }
      std::set<int> roots;
      for (int g = 0; g < k; ++g) roots.insert(comp.find(g));
      rep.components = roots.size();
      rep.connected = (rep.components == 1);

      // Shape (Bestvina-Handel Props. 3.3.3-3.3.4): at least two gates;
      // for k = 2 exactly one infinitesimal edge; for k > 2 the edges join
      // cyclically adjacent gates and form a k-gon or a k-gon minus one
      // side.  Gates are numbered in cyclic order, so adjacency is |ga-gb|
      // = 1 or {0,k-1}.
      rep.shape_ok = (k >= 2);
      if (k == 2) rep.shape_ok = rep.shape_ok && (joined_gates.size() == 1);
      else if (k > 2)
        {
          const int nj = joined_gates.size();
          if (nj != k && nj != k-1) rep.shape_ok = false;
          for (const auto& jg : joined_gates)
            {
              const int diff = jg.second - jg.first;
              if (!(diff == 1 || diff == k-1)) rep.shape_ok = false;
            }
        }
      if (!rep.connected) result.connected = false;
      if (!rep.shape_ok) result.shapes_ok = false;
      result.vertices.push_back(rep);
    }

  return result;
}

} // anonymous namespace


gate_analysis gate_accumulator::analyse(const ttnumbering& N) const
{
  if (ng != N.nedges() + N.nprongs()) fail("Numbering does not match the path");
  return analyse_impl(N,Dacc,Tacc);
}


gate_analysis analyse_gates(const ttnumbering& N, const jlt::freeauto<int>& AM)
{
  fold_derivative d(AM,N);
  return analyse_impl(N,d.D,d.turns);
}


std::ostream& gate_analysis::print(std::ostream& strm) const
{
  for (const gate_vertex_report& v : vertices)
    {
      strm << v.name << ": gates";
      for (const auto& g : v.gates)
        {
          strm << " {";
          for (std::size_t i = 0; i < g.size(); ++i) strm << (i ? " " : "") << g[i];
          strm << "}";
        }
      strm << "; joins";
      for (const turn& t : v.joins) strm << " " << t.first << "~" << t.second;
      strm << "; " << (v.connected ? "connected" : "NOT connected")
           << " (" << v.components << " component" << (v.components == 1 ? "" : "s") << ")"
           << (v.shape_ok ? "" : "; unexpected shape") << "\n";
    }
  strm << "gate test: " << (connected ? "connected at every vertex (pseudo-Anosov)"
                                      : "disconnected vertex found (reducible)")
       << (shapes_ok ? "" : "; some vertex has an unexpected infinitesimal-edge shape")
       << "\n";
  return strm;
}

} // namespace traintracks
