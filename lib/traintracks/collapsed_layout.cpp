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

// The longest control length along d from end A that keeps the control
// point within the x-range of the piece's two ends A and Z, widened by a
// margin.  A cubic lies within the hull of its control points, so with
// both controls inside that range it cannot run past Z and double back.
double max_control(const vec2& A, const vec2& d, const vec2& Z)
{
  const double margin = 0.1;
  if (std::fabs(d.x) < 1e-12) return 1e300;
  const double lo = std::min(A.x,Z.x) - margin, hi = std::max(A.x,Z.x) + margin;
  return ((d.x > 0 ? hi : lo) - A.x)/d.x;
}

// Edges leaving one prong share their tangent there, so which of them
// lies to the left of which is decided by their curvature at that point.
// Slots run anticlockwise, which looking out along the prong is leftward,
// so the curvature must increase along the slots.  For a cubic leaving
// end point A along d with control length h and next control B, the
// curvature there is (2/3) cross(d,B-A)/h^2, so it is put in order by
// shortening h when B lies to the left, or lengthening it when to the
// right.
void order_curvatures(const ttnumbering& num, collapsed_layout& L)
{
  for (int q = 0; q < num.nprongs(); ++q)
    {
      const std::vector<int>& letters = num.prong_letters[q];
      if (letters.size() < 2) continue;
      const vec2 d = L.prong_dir[q];
      const auto end_of = [&](const int letter, vec2*& A, vec2*& C, vec2*& B,
                              vec2** Z = nullptr)
        {
          const int e = num.edge_of(letter);
          vec2* z;
          if (letter > 0)
            { cubic& c = L.arc[e].front(); A = &c.p0; C = &c.p1; B = &c.p2; z = &c.p3; }
          else
            { cubic& c = L.arc[e].back(); A = &c.p3; C = &c.p2; B = &c.p1; z = &c.p0; }
          if (Z) *Z = z;
        };
      const auto lateral = [&](const int letter)
        {
          vec2 *A, *C, *B;
          end_of(letter,A,C,B);
          return d.x*(B->y-A->y) - d.y*(B->x-A->x);
        };
      const auto curvature = [&](const int letter)
        {
          vec2 *A, *C, *B;
          end_of(letter,A,C,B);
          const double h = std::hypot(C->x-A->x,C->y-A->y);
          return (2.0/3.0)*lateral(letter)/(h*h);
        };
      // Each end's control length may move between 0.2 and 3 times where
      // it started.  A pair out of order is moved apart from both sides,
      // the later end up and the earlier down, since one end alone often
      // cannot get there: lengthening only takes a curvature towards zero,
      // never past it.
      const auto length = [&](const int letter)
        {
          vec2 *A, *C, *B;
          end_of(letter,A,C,B);
          return std::hypot(C->x-A->x,C->y-A->y);
        };
      const auto scale_by = [&](const int letter, const double f)
        {
          vec2 *A, *C, *B;
          end_of(letter,A,C,B);
          C->x = A->x + f*(C->x-A->x);
          C->y = A->y + f*(C->y-A->y);
        };
      std::vector<double> h0(letters.size());
      for (int i = 0; i < (int)letters.size(); ++i) h0[i] = length(letters[i]);
      // Nudge the curvature of end i up (dir = +1) or down (dir = -1),
      // unless that would take its control length out of bounds.
      const auto nudge = [&](const int i, const int dir)
        {
          const double lat = lateral(letters[i]);
          if (std::fabs(lat) < 1e-12) return;
          // Shortening raises |curvature|, lengthening lowers it.
          const bool shorten = (lat > 0) == (dir > 0);
          const double f = shorten ? 0.85 : 1.2;
          const double h = length(letters[i]);
          if (h*f < 0.2*h0[i] || h*f > 3.0*h0[i]) return;
          vec2 *A, *C, *B, *Z;
          end_of(letters[i],A,C,B,&Z);
          if (f > 1 && h*f > max_control(*A,d,*Z)) return;
          scale_by(letters[i],f);
        };
      for (int pass = 0; pass < 60; ++pass)
        {
          bool sorted = true;
          for (int i = 1; i < (int)letters.size(); ++i)
            {
              const double prev = curvature(letters[i-1]);
              const double k = curvature(letters[i]);
              if (k > std::max(prev + 0.1,1.2*prev + 0.1)) continue;
              sorted = false;
              nudge(i,+1);
              nudge(i-1,-1);
            }
          if (sorted) break;
        }
    }
}

} // namespace

vec2 cubic_point(const cubic& c, const double t)
{
  const double u = 1.0-t;
  const double w0 = u*u*u, w1 = 3*u*u*t, w2 = 3*u*t*t, w3 = t*t*t;
  return { w0*c.p0.x + w1*c.p1.x + w2*c.p2.x + w3*c.p3.x,
           w0*c.p0.y + w1*c.p1.y + w2*c.p2.y + w3*c.p3.y };
}

// The drawing is built as an arc diagram over the axis.
//
// Rooted at puncture 1, where the boundary walk is cut, every other
// multigon cuts off with the edge to its parent the block of punctures
// lo..hi of its subtree.  Those blocks are nested or disjoint, and each
// gets a box over its stretch of axis, narrowed a little more at each
// depth so that a box lies strictly inside the box around it and
// siblings do not touch.  A multigon sits inside its own box, above the
// boxes of its children, and an edge crosses exactly one box boundary,
// the top of its child's box, straight above the child and vertically.
// Each edge is two cubics joined smoothly there, one inside the child's
// box and one in the parent's, so curves in different boxes cannot meet.
//
// The boxes are built bottom up, since a box has to be tall enough for
// what is drawn inside it: within the parent's box, the piece from a
// child passes over every sibling between it and the parent, and must
// clear both that sibling's box and that sibling's own piece.  So the
// pieces at a parent are built nearest first, each raised above the
// highest point actually reached, within its span, by what it passes
// over, and the parent's box is then closed above all of them.
collapsed_layout make_collapsed_layout(const ttnumbering& num,
                                       const tt_embedding& emb)
{
  const int M = (int)num.prong_number.size(), E = num.nedges();
  const double unit = 0.5;
  const vec2 up = { 0.0, 1.0 };
  const int nsample = 32;

  // The sense in which the prong index runs round a multigon in the
  // plane.  It has to be the sense of the boundary walk that ordered the
  // punctures, or the prongs come out mirrored while the edges within a
  // prong do not, and no drawing avoids a crossing.
  const double sense = 1.0;

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
  std::vector<int> lo(M,M+1), hi(M,0), parent(M,-2), depth(M,0), order;
  {
    std::vector<int> stack(1,root);
    parent[root] = -1;
    while (!stack.empty())
      {
        const int m = stack.back(); stack.pop_back();
        order.push_back(m);
        for (int k = 0; k < (int)adj[m].size(); ++k)
          if (parent[adj[m][k]] == -2)
            {
              parent[adj[m][k]] = m;
              depth[adj[m][k]] = depth[m]+1;
              stack.push_back(adj[m][k]);
            }
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

  std::vector<int> edge_to_parent(M,-1);
  for (int e = 0; e < E; ++e)
    {
      const int t = L.edge_mg[e].first, h = L.edge_mg[e].second;
      edge_to_parent[parent[t] == h ? t : h] = e;
    }

  std::vector<double> box_lo(M), box_hi(M), box_top(M,0.0);
  for (int m = 0; m < M; ++m)
    {
      const double shrink = 0.5*(1.0 - std::pow(0.7,depth[m]));
      box_lo[m] = lo[m] - 0.5 + shrink;
      box_hi[m] = hi[m] + 0.5 - shrink;
    }

  // A puncture sits on the axis and the track lives above it, so its one
  // prong points straight up and the loop hangs below.
  L.prong_dir.assign(num.nprongs(),up);
  L.arc.assign(E,std::vector<cubic>());
  std::vector<cubic> outer_of(M);

  for (int i = (int)order.size()-1; i >= 0; --i)
    {
      const int m = order[i];
      std::vector<int> kids;
      for (int o : adj[m]) if (parent[o] == m) kids.push_back(o);

      double kids_top = 0.0;
      for (int o : kids) kids_top = std::max(kids_top,box_top[o]);
      if (L.puncture_pos[m])
        L.mg[m] = { (double)L.puncture_pos[m], 0.0 };
      else
        L.mg[m] = { 0.5*(lo[m]+hi[m]), kids_top + 0.5*unit };

      // A multigon's prongs are equally spaced by angle around it, and
      // every edge at a prong leaves along that prong's direction, so
      // edges at the same prong are tangent.  Only the rotation of the
      // star is free, and it is the circular mean of where the edges go
      // next: straight up for the edge to the parent, and to the point
      // where it enters the child's box for the others.
      if (!L.puncture_pos[m])
        {
          const int k = (int)num.prong_number[m].size();
          double sx = 0.0, sy = 0.0;
          for (int j = 0; j < k; ++j)
            {
              const int q = num.prong_number[m][j];
              for (const int letter : num.prong_letters[q])
                {
                  const int e = num.edge_of(letter);
                  const int o = (L.edge_mg[e].first == m ? L.edge_mg[e].second
                                                         : L.edge_mg[e].first);
                  const double bearing = (o == parent[m]) ? 0.5*M_PI
                    : std::atan2(box_top[o]-L.mg[m].y,L.mg[o].x-L.mg[m].x);
                  const double t = bearing - sense*2.0*M_PI*j/k;
                  sx += std::cos(t); sy += std::sin(t);
                }
            }
          const double theta = std::atan2(sy,sx);
          for (int j = 0; j < k; ++j)
            {
              const double a = theta + sense*2.0*M_PI*j/k;
              L.prong_dir[num.prong_number[m][j]] = { std::cos(a), std::sin(a) };
            }
        }

      // The pieces from the children, nearest first.
      std::sort(kids.begin(),kids.end(),[&](const int u, const int w)
        { return std::fabs(L.mg[u].x-L.mg[m].x) < std::fabs(L.mg[w].x-L.mg[m].x); });
      double top = std::max(kids_top,L.mg[m].y);
      std::vector<int> done;
      for (int c : kids)
        {
          const int e = edge_to_parent[c];
          const bool tail_is_child = (L.edge_mg[e].first == c);
          const vec2 dc = L.prong_dir[tail_is_child ? L.edge_pr[e].first : L.edge_pr[e].second];
          const vec2 dp = L.prong_dir[tail_is_child ? L.edge_pr[e].second : L.edge_pr[e].first];
          const vec2 C = L.mg[c], P = L.mg[m], X = { C.x, box_top[c] };
          const double a = std::min(C.x,P.x), b = std::max(C.x,P.x);

          // The highest point of what this piece passes over.
          double clear = X.y;
          for (int o : done)
            {
              if (std::min(box_hi[o],b) <= std::max(box_lo[o],a)) continue;
              clear = std::max(clear,box_top[o]);
              for (int s = 0; s <= nsample; ++s)
                {
                  const vec2 q = cubic_point(outer_of[o],(double)s/nsample);
                  if (q.x > a && q.x < b) clear = std::max(clear,q.y);
                }
            }
          const double peak = clear + (clear > X.y ? 0.5*unit : 0.0);

          cubic inner;
          const double hi_ = 0.5*std::hypot(X.x-C.x,X.y-C.y);
          inner.p0 = C;
          const double hc = std::min(hi_,max_control(C,dc,X));
          inner.p1 = { C.x + hc*dc.x, C.y + hc*dc.y };
          inner.p2 = { X.x - hi_*up.x, X.y - hi_*up.y };
          inner.p3 = X;

          // Controls at height Y put the middle of the cubic at the peak;
          // with nothing to clear they stay short, so that the piece does
          // not rise above what it has to.
          const double ho = 0.5*std::hypot(P.x-X.x,P.y-X.y);
          const double Y = (peak - (X.y+P.y)/8.0)/0.75;
          double hx = std::max(0.25*ho,Y - X.y);
          double hp = (dp.y > 0.2) ? std::max(0.25*ho,(Y - P.y)/dp.y) : ho;

          // That only aims the middle.  Check the whole piece against what
          // it must stay out of or above -- the axis, every child's box,
          // and the pieces already built.  Near the parent, where sibling pieces
          // meet it by construction, only the boxes count.
          const auto dips = [&](const cubic& piece)
            {
              for (int s = 1; s < nsample; ++s)
                {
                  const vec2 q = cubic_point(piece,(double)s/nsample);
                  const bool near_p = std::hypot(q.x-P.x,q.y-P.y) < 0.15;
                  // The track lives in the upper half plane.
                  if (q.y < 0.0) return true;
                  if (q.x > box_lo[c] && q.x < box_hi[c] && q.y < box_top[c]-1e-9)
                    return true;
                  for (int o : kids)
                    if (o != c && q.x > box_lo[o] && q.x < box_hi[o] && q.y < box_top[o])
                      return true;
                  for (int o : done)
                    {
                      if (near_p) continue;
                      for (int r = 0; r < nsample; ++r)
                        {
                          const vec2 u = cubic_point(outer_of[o],(double)r/nsample);
                          const vec2 w = cubic_point(outer_of[o],(double)(r+1)/nsample);
                          if ((u.x-q.x)*(w.x-q.x) <= 0 && std::fabs(w.x-u.x) > 1e-12)
                            {
                              const double yo = u.y + (w.y-u.y)*(q.x-u.x)/(w.x-u.x);
                              if (q.y < yo + 0.1*unit) return true;
                            }
                        }
                    }
                }
              return false;
            };
          // Search the two control lengths over a range of scales, both
          // shorter and longer than the first guess, and keep the lowest
          // piece that clears.  Only lengthening, as before, sent pieces
          // hundreds of units up.  With nothing that clears, fall back on
          // lengthening until something does.
          const auto make = [&](const double a, const double b)
            {
              cubic q;
              q.p0 = X;
              q.p1 = { X.x, X.y + a };
              const double bb = std::min(b,max_control(P,dp,X));
              q.p2 = { P.x + bb*dp.x, P.y + bb*dp.y };
              q.p3 = P;
              return q;
            };
          const auto height = [&](const cubic& q)
            {
              double y = 0.0;
              for (int s = 0; s <= nsample; ++s)
                y = std::max(y,cubic_point(q,(double)s/nsample).y);
              return y;
            };
          static const double scale[] =
            { 0.2, 0.3, 0.45, 0.65, 1.0, 1.5, 2.2, 3.3, 5.0 };
          cubic outer = make(hx,hp);
          bool found = false;
          double best = 0.0;
          for (const double fa : scale)
            for (const double fb : scale)
              {
                const cubic q = make(fa*hx,fb*hp);
                if (dips(q)) continue;
                const double y = height(q);
                if (!found || y < best - 1e-9) { outer = q; best = y; found = true; }
              }
          for (int it = 0; !found && it < 30; ++it)
            {
              hx *= 1.2;
              if (dp.y > 0.2) hp *= 1.2;
              outer = make(hx,hp);
              found = !dips(outer);
            }
          outer_of[c] = outer;
          done.push_back(c);
          for (int s = 0; s <= nsample; ++s)
            top = std::max(top,cubic_point(outer,(double)s/nsample).y);

          if (tail_is_child)
            L.arc[e] = { inner, outer };
          else
            {
              std::swap(outer.p0,outer.p3); std::swap(outer.p1,outer.p2);
              std::swap(inner.p0,inner.p3); std::swap(inner.p1,inner.p2);
              L.arc[e] = { outer, inner };
            }
        }

      // Close the box above everything drawn in it.
      box_top[m] = top + (kids.empty() ? unit : 0.5*unit);
    }

  order_curvatures(num,L);

  return L;
}

} // namespace traintracks
