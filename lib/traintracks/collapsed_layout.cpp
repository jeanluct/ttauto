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

// The collapsed drawing of a train track.  See
// include/traintracks/collapsed_layout.hpp for what it guarantees.

#include <algorithm>
#include <cmath>
#include <vector>

#include "traintracks/collapsed_layout.hpp"

namespace traintracks {

namespace {

int multigon_of(const ttnumbering& num, const int q)
{
  return num.prong[q].multigon;
}

// The control point at each end runs along that end's prong direction,
// so edges sharing a prong are tangent there.  Its length grows with the
// horizontal span, so an edge reaching over intervening structure still
// arcs above it rather than cutting through.
void edge_controls(const vec2& a, const vec2& b,
                   const vec2& da, const vec2& db, vec2& c1, vec2& c2)
{
  const double h = 0.55*std::fabs(a.x-b.x) + 0.30;
  c1 = { a.x + h*da.x, a.y + h*da.y };
  c2 = { b.x + h*db.x, b.y + h*db.y };
}

} // namespace

vec2 cubic_point(const cubic& c, const double t)
{
  const double u = 1.0-t;
  const double w0 = u*u*u, w1 = 3*u*u*t, w2 = 3*u*t*t, w3 = t*t*t;
  return { w0*c.p0.x + w1*c.p1.x + w2*c.p2.x + w3*c.p3.x,
           w0*c.p0.y + w1*c.p1.y + w2*c.p2.y + w3*c.p3.y };
}

// Place the multigons from the embedding.
//
// Punctures sit on the axis at the position the boundary walk gives them.
// An unpunctured multigon sits over the middle of the puncture positions
// in its subtree, at a height that grows with the width of that span.
// The walk order is the planar order, so those spans nest instead of
// interleaving, and a wider structure therefore sits above whatever it
// contains.
collapsed_layout make_collapsed_layout(const ttnumbering& num,
                                       const tt_embedding& emb)
{
  const int M = (int)num.prong_number.size(), E = num.nedges();

  collapsed_layout L;
  L.mg.assign(M,{0.0,0.0});
  L.puncture_pos.assign(M,0);
  L.edge_mg.assign(E,std::make_pair(-1,-1));
  L.edge_pr.assign(E,std::make_pair(-1,-1));
  for (int e = 0; e < E; ++e)
    {
      L.edge_mg[e] = std::make_pair(multigon_of(num,num.edge_tail[e]),
                                    multigon_of(num,num.edge_head[e]));
      L.edge_pr[e] = std::make_pair(num.edge_tail[e],num.edge_head[e]);
    }

  for (int i = 0; i < emb.npunctures(); ++i)
    L.puncture_pos[multigon_of(num,emb.puncture_order[i])] = i+1;

  std::vector<std::vector<int> > adj(M);
  for (int e = 0; e < E; ++e)
    {
      adj[L.edge_mg[e].first].push_back(L.edge_mg[e].second);
      adj[L.edge_mg[e].second].push_back(L.edge_mg[e].first);
    }

  // The multigons and main edges form a tree, so rooting it and taking
  // the span of each subtree needs no search.
  const int root = multigon_of(num,emb.puncture_order[0]);
  std::vector<int> lo(M,M+1), hi(M,0), parent(M,-2), order;
  {
    std::vector<int> stack(1,root);
    parent[root] = -1;
    while (!stack.empty())
      {
        const int m = stack.back(); stack.pop_back();
        order.push_back(m);
        for (int k = 0; k < (int)adj[m].size(); ++k)
          if (parent[adj[m][k]] == -2)
            { parent[adj[m][k]] = m; stack.push_back(adj[m][k]); }
      }
  }
  for (int i = (int)order.size()-1; i >= 0; --i)
    {
      const int m = order[i];
      if (L.puncture_pos[m])
        {
          lo[m] = std::min(lo[m],L.puncture_pos[m]);
          hi[m] = std::max(hi[m],L.puncture_pos[m]);
        }
      if (parent[m] >= 0)
        {
          lo[parent[m]] = std::min(lo[parent[m]],lo[m]);
          hi[parent[m]] = std::max(hi[parent[m]],hi[m]);
        }
    }

  for (int m = 0; m < M; ++m)
    {
      if (L.puncture_pos[m])
        L.mg[m] = { (double)L.puncture_pos[m], 0.0 };
      else
        L.mg[m] = { 0.5*(lo[m]+hi[m]), 0.55*(hi[m]-lo[m]) };
    }

  // A multigon's prongs are equally spaced by angle around it, and the
  // edges are then attached to the prongs, with every edge at a prong
  // leaving along that prong's direction.  Edges at the same prong are
  // therefore tangent.
  //
  // Only the overall rotation of the star is free.  Pointing each prong
  // at its own neighbours instead lets the prongs bunch together, and the
  // multigon stops looking like a k-gon at all.
  L.prong_dir.assign(num.nprongs(),{0.0,1.0});

  // The sense in which the prong index runs round a multigon in the
  // plane.  It has to be the sense of the boundary walk that ordered the
  // punctures, or the prongs come out mirrored while the edges within a
  // prong do not, and no drawing avoids a crossing.
  const double sense = 1.0;

  // The rotation comes from the gaps between the edges.  Going round a
  // multigon anticlockwise, each edge leads to a block of consecutive
  // punctures (the edge towards puncture 1 to all the rest, cyclically),
  // so the gap between two consecutive edges faces a definite place:
  // the stretch of axis between punctures i and i+1, or straight up for
  // the gap between puncture n and puncture 1, where the walk is cut.
  // The rotation is the circular mean of what those gaps ask for.
  const int npunct = emb.npunctures();
  for (int m = 0; m < M; ++m)
    {
      // A puncture sits on the axis and the track lives above it, so its
      // one prong points straight up and the loop hangs below.
      if (L.puncture_pos[m]) continue;

      // The edges at m anticlockwise, with the angle each leaves at
      // relative to the star's rotation, and its block [first,last].
      struct end_at { double angle; int first, last; };
      std::vector<end_at> ends;
      const int k = (int)num.prong_number[m].size();
      for (int j = 0; j < k; ++j)
        {
          const int q = num.prong_number[m][j];
          for (const int letter : num.prong_letters[q])
            {
              const int e = num.edge_of(letter);
              const int o = (L.edge_mg[e].first == m ? L.edge_mg[e].second
                                                     : L.edge_mg[e].first);
              end_at x;
              x.angle = sense*2.0*M_PI*j/k;
              if (o == parent[m])
                { x.first = hi[m]+1; x.last = lo[m]-1; }
              else
                { x.first = lo[o]; x.last = hi[o]; }
              ends.push_back(x);
            }
        }

      double sx = 0.0, sy = 0.0;
      for (int i = 0; i < (int)ends.size(); ++i)
        {
          const end_at& u = ends[i];
          const end_at& w = ends[(i+1) % ends.size()];
          double off = w.angle;
          if (off < u.angle || (i+1 == (int)ends.size() && off == u.angle))
            off += 2.0*M_PI*sense;
          off = 0.5*(u.angle + off);
          const vec2 at = L.mg[m];
          double bearing;
          if (u.last == npunct || w.first == 1 || w.first > npunct)
            bearing = 0.5*M_PI;
          else
            bearing = std::atan2(-at.y,(u.last+0.5)-at.x);
          const double t = bearing - off;
          sx += std::cos(t); sy += std::sin(t);
        }
      const double theta = std::atan2(sy,sx);

      for (int j = 0; j < k; ++j)
        {
          const double a = theta + sense*2.0*M_PI*j/k;
          L.prong_dir[num.prong_number[m][j]] = { std::cos(a), std::sin(a) };
        }
    }

  L.arc.resize(E);
  for (int e = 0; e < E; ++e)
    {
      cubic& c = L.arc[e];
      c.p0 = L.mg[L.edge_mg[e].first];
      c.p3 = L.mg[L.edge_mg[e].second];
      edge_controls(c.p0,c.p3,L.prong_dir[L.edge_pr[e].first],
                    L.prong_dir[L.edge_pr[e].second],c.p1,c.p2);
    }

  return L;
}

} // namespace traintracks
