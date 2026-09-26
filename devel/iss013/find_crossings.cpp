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

// Every automaton vertex for n punctures whose collapsed drawing
// (traintracks::make_collapsed_layout) has a proper crossing.  One line
// per track, fields separated by '|':
//
//   stratum | vertex | height of drawing | x,y;x,y;... crossings | coding
//
// Stratum and vertex are 1-based.  The crossing test is the one in
// testsuite/traintracks/test_collapsed_layout.cpp, without its speedups.
// Build and run from the repository root; see devel/iss013/README.md.
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>
#include "traintracks/build.hpp"
#include "traintracks/collapsed_layout.hpp"
#include "traintracks/traintrack.hpp"
#include "ttauto/ttfoldgraph.hpp"
using namespace traintracks;
static int orient(const vec2& p, const vec2& q, const vec2& r)
{ const double v=(q.x-p.x)*(r.y-p.y)-(q.y-p.y)*(r.x-p.x); return (v>1e-12)-(v<-1e-12); }
static bool sx(const vec2&a,const vec2&b,const vec2&c,const vec2&d)
{ return orient(a,b,c)*orient(a,b,d)<0 && orient(c,d,a)*orient(c,d,b)<0; }
static bool eq(const vec2&p,const vec2&q){return p.x==q.x&&p.y==q.y;}
static bool cross(const cubic& A,const cubic& B, vec2& at)
{
  std::vector<vec2> sh;
  if (eq(A.p0,B.p0)||eq(A.p0,B.p3)) sh.push_back(A.p0);
  if (eq(A.p3,B.p0)||eq(A.p3,B.p3)) sh.push_back(A.p3);
  auto nr=[&](const vec2&p){for(auto&s:sh) if((p.x-s.x)*(p.x-s.x)+(p.y-s.y)*(p.y-s.y)<4e-4) return true; return false;};
  const int S=60; std::vector<vec2> a(S+1),b(S+1);
  for(int i=0;i<=S;++i){a[i]=cubic_point(A,(double)i/S);b[i]=cubic_point(B,(double)i/S);}
  for(int i=0;i<S;++i){ if(nr(a[i])||nr(a[i+1])) continue;
    for(int j=0;j<S;++j){ if(nr(b[j])||nr(b[j+1])) continue;
      if(sx(a[i],a[i+1],b[j],b[j+1])){ at={0.5*(a[i].x+a[i+1].x),0.5*(a[i].y+a[i+1].y)}; return true;}}}
  return false;
}
int main(int argc, char** argv)
{
  const int n = (argc > 1 ? std::atoi(argv[1]) : 8);
  jlt::vector<traintrack> ttv=build_traintrack_list(n);
  for (int s=0;s<(int)ttv.size();++s){
    const ttauto::ttfoldgraph<traintrack> g(ttv[s]);
    for (int v=0;v<g.vertices();++v){
      const traintrack& tt=g.traintrack(v);
      const ttnumbering num=tt.numbering();
      const collapsed_layout L=make_collapsed_layout(num,outer_embedding(num));
      std::vector<vec2> pts; double ymax=0;
      for (auto& arc : L.arc) for (auto& c : arc)
        for (int i=0;i<=20;++i) ymax=std::max(ymax,cubic_point(c,i/20.0).y);
      for (int i=0;i<(int)L.arc.size();++i) for (int j=i+1;j<(int)L.arc.size();++j){
        bool got=false;
        for (auto& A:L.arc[i]) for (auto& B:L.arc[j]){ vec2 at; if(!got && cross(A,B,at)){pts.push_back(at);got=true;} }
      }
      if (pts.empty()) continue;
      std::ostringstream o; tt.print_coding(o);
      std::cout<<s+1<<"|"<<v+1<<"|"<<ymax<<"|";
      for (auto&p:pts) std::cout<<p.x<<","<<p.y<<";";
      std::cout<<"|"<<o.str()<<"\n";
    }}
}
