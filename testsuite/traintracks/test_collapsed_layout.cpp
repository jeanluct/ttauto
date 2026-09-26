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

// Planarity of the collapsed drawing (collapsed_layout.hpp).
//
// Arcs of a train track never cross; where two meet they are tangent.  In
// the collapsed drawing every edge at a multigon ends at that multigon's
// point, so nearly every pair of arcs shares an endpoint, and a checker
// that excuses such pairs wholesale tests almost nothing.  An earlier
// checker did exactly that and reported no crossings where there were
// hundreds.  This one excuses only the samples right at the shared point,
// and the self-check below proves it sees a crossing between two arcs
// leaving the same point.

#include <iostream>
#include <string>
#include <vector>
#include <jlt/vector.hpp>
#include "check.hpp"
#include "traintracks/build.hpp"
#include "traintracks/collapsed_layout.hpp"
#include "traintracks/traintrack.hpp"
#include "ttauto/ttfoldgraph.hpp"

using traintracks::cubic;
using traintracks::traintrack;
using traintracks::vec2;

static const int samples = 60;

// Sign of the turn p -> q -> r, with a dead band so that collinear or
// touching samples do not count as crossing.
static int orient(const vec2& p, const vec2& q, const vec2& r)
{
  const double v = (q.x-p.x)*(r.y-p.y) - (q.y-p.y)*(r.x-p.x);
  return (v > 1e-12) - (v < -1e-12);
}

static bool segments_cross(const vec2& a, const vec2& b,
                           const vec2& c, const vec2& d)
{
  return orient(a,b,c)*orient(a,b,d) < 0 && orient(c,d,a)*orient(c,d,b) < 0;
}

static bool same_point(const vec2& p, const vec2& q)
{
  return p.x == q.x && p.y == q.y;
}

// Do two arcs properly cross?  Samples within radius r of an endpoint
// the two arcs share are dropped: the arcs meet there by construction.
static bool arcs_cross(const cubic& A, const cubic& B, const double r = 0.02)
{
  std::vector<vec2> shared;
  if (same_point(A.p0,B.p0) || same_point(A.p0,B.p3)) shared.push_back(A.p0);
  if (same_point(A.p3,B.p0) || same_point(A.p3,B.p3)) shared.push_back(A.p3);
  auto near = [&](const vec2& p)
    {
      for (const vec2& s : shared)
        if ((p.x-s.x)*(p.x-s.x) + (p.y-s.y)*(p.y-s.y) < r*r) return true;
      return false;
    };

  std::vector<vec2> a(samples+1), b(samples+1);
  for (int i = 0; i <= samples; ++i)
    {
      a[i] = traintracks::cubic_point(A,(double)i/samples);
      b[i] = traintracks::cubic_point(B,(double)i/samples);
    }
  for (int i = 0; i < samples; ++i)
    {
      if (near(a[i]) || near(a[i+1])) continue;
      for (int j = 0; j < samples; ++j)
        {
          if (near(b[j]) || near(b[j+1])) continue;
          if (segments_cross(a[i],a[i+1],b[j],b[j+1])) return true;
        }
    }
  return false;
}

// Number of pairs of arcs of the track's drawing that cross.
static int crossings(const traintrack& tt)
{
  const traintracks::ttnumbering num = tt.numbering();
  const traintracks::collapsed_layout L =
    traintracks::make_collapsed_layout(num,traintracks::outer_embedding(num));
  int n = 0;
  for (int i = 0; i < (int)L.arc.size(); ++i)
    for (int j = i+1; j < (int)L.arc.size(); ++j)
      {
        bool x = false;
        for (const cubic& a : L.arc[i])
          for (const cubic& b : L.arc[j])
            if (!x && arcs_cross(a,b)) x = true;
        if (x) ++n;
      }
  return n;
}

int main()
{
  using std::cout;
  using std::endl;

  //
  // Self-check.  Two arcs leaving the origin, one straight up the
  // diagonal and one swinging out to the right and back over it, cross.
  // Two arcs leaving the origin tangent, then parting left and right, do
  // not.  Arcs with no endpoint in common are tested in full.
  //
  {
    const cubic diag  = { {0,0}, {1,1}, {2,2}, {3,3} };
    const cubic swing = { {0,0}, {3,0}, {3,0}, {0,3} };
    CHECK(arcs_cross(diag,swing));
    CHECK(arcs_cross(swing,diag));

    const cubic left  = { {0,0}, {0,1}, {-1,1}, {-1,0} };
    const cubic right = { {0,0}, {0,1}, { 1,1}, { 1,0} };
    CHECK(!arcs_cross(left,right));

    const cubic high = { {-1,1}, {0,1}, {0,1}, {1,1} };
    const cubic low  = { {-1,2}, {0,0}, {0,0}, {1,2} };
    CHECK(arcs_cross(high,low));
  }

  //
  // The two tracks of Figure 4 of the ttauto paper, (a) and (b):
  // vertices 2 and 0 of the automaton of build_traintrack_list(4)[1].
  //
  const char* fig4[] = {
    "1111 1122 1112 1311 2311 1111 3311 1111",
    "1111 1311 2311 1111 3312 1111 3322 1111"
  };
  for (const char* c : fig4)
    {
      const int x = crossings(traintrack(c));
      CHECK_MSG(x == 0, c << ": " << x << " crossings");
    }
  cout << "Figure 4 tracks: no crossings" << endl;

  //
  // Every vertex of every automaton for n = 3..6.  At n = 7, 98 of 3272
  // tracks still cross (2026-09-26): a prong carrying several children
  // can be rotated to point up, and those edges then have to climb over
  // everything to come back down.
  //
  int ntracks = 0, nbad = 0, ncross = 0;
  for (int n = 3; n <= 6; ++n)
    {
      jlt::vector<traintrack> ttv = traintracks::build_traintrack_list(n);
      for (int s = 0; s < (int)ttv.size(); ++s)
        {
          const ttauto::ttfoldgraph<traintrack> ttg(ttv[s]);
          for (int v = 0; v < ttg.vertices(); ++v)
            {
              const int x = crossings(ttg.traintrack(v));
              ++ntracks;
              if (x) { ++nbad; ncross += x; }
            }
        }
    }
  cout << ntracks << " tracks for n = 3..6: " << nbad << " with crossings, "
       << ncross << " crossing pairs in all" << endl;
  CHECK(ntracks == 428);
  CHECK(nbad == 0);

  cout << "\ntest_collapsed_layout: all checks passed" << endl;
  return 0;
}
