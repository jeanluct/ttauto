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

#include <iostream>
#include <cstdlib>
#include "traintracks/coding.hpp"
#include "traintracks/util.hpp"

namespace traintracks {
namespace coding {

namespace {

inline void append_block(traintrack::intVec& code,
                         const int prong, const int nprongs,
                         const int label, const int edge,
                         const int nedges)
{
  code.push_back(prong);
  code.push_back(nprongs);
  code.push_back(label);
  code.push_back(edge);
  code.push_back(nedges);
}

inline int multigon_index(const traintrack& tt, const multigon* mm)
{
  for (int m = 0; m < tt.multigons(); ++m)
    {
      if (&tt.Multigon(m) == mm) return m;
    }
  std::cerr << "Error in coding::multigon_index(): multigon not found.\n";
  std::exit(1);
}

void recursive_coding(const traintrack& tt, const multigon& mm,
                      const int pin, const int ein,
                      traintrack::intVec& code, const int dir)
{
  int p = pin, e = ein;

  if (mm.edges() == 1)
    {
      append_block(code,0,1,mm.label(),0,1);
      return;
    }

  do
    {
      const int prong = traintracks::mod(dir*(p-pin),mm.prongs());
      const int nprongs = mm.prongs();
      const int label = mm.label();
      const int nedges = mm.edges(p);
      const int edge = (dir == 1 ? e : nedges-1-e);

      append_block(code,prong,nprongs,label,edge,nedges);

      int pout, eout;
      multigon *ed = mm.Edge(p,e)->target_multigon(&mm,pout,eout);

      if (!(p == pin && e == ein))
        {
          recursive_coding(tt,*ed,pout,eout,code,dir);
        }

      mm.cycle_edges(p,e,dir);
    }
  while (!(p == pin && e == ein));
}

void recursive_canonical_coding(const traintrack& tt, const multigon& mm,
                                const int pin, const int ein,
                                traintrack::intVec& code, const int dir,
                                jlt::vector<int>& first_pin)
{
  const int mi = multigon_index(tt,&mm);
  if (first_pin[mi] < 0) first_pin[mi] = pin;
  const int p0 = first_pin[mi];

  int p = pin, e = ein;

  if (mm.edges() == 1)
    {
      append_block(code,0,1,mm.label(),0,1);
      return;
    }

  do
    {
      const int prong = traintracks::mod(dir*(p-p0),mm.prongs());
      const int nprongs = mm.prongs();
      const int label = mm.label();
      const int nedges = mm.edges(p);
      const int edge = (dir == 1 ? e : nedges-1-e);

      append_block(code,prong,nprongs,label,edge,nedges);

      int pout, eout;
      multigon *ed = mm.Edge(p,e)->target_multigon(&mm,pout,eout);

      if (!(p == pin && e == ein))
        {
          recursive_canonical_coding(tt,*ed,pout,eout,code,dir,first_pin);
        }

      mm.cycle_edges(p,e,dir);
    }
  while (!(p == pin && e == ein));
}

} // namespace

traintrack::intVec coding_from_monogon(const traintrack& tt,
                                       const int mono,
                                       const int dir)
{
  if (tt.Multigon(mono).edges() > 1)
    {
      std::cerr << "Not an uncusped monogon in traintracks::coding::coding_from_monogon.\n";
      std::exit(1);
    }

  traintrack::intVec code;
  append_block(code,0,1,tt.Multigon(mono).label(),0,1);

  int pmono, pemono;
  multigon *egmono =
    tt.Multigon(mono).Edge(0,0)->target_multigon(&tt.Multigon(mono),pmono,pemono);

  recursive_coding(tt,*egmono,pmono,pemono,code,dir);

  return code;
}

traintrack::intVec coding(const traintrack& tt, const int dir)
{
  int mono = 0;
  traintrack::intVec codemin = coding_from_monogon(tt,mono,dir);

  while (tt.Multigon(++mono).edges() <= 1)
    {
      traintrack::intVec code = coding_from_monogon(tt,mono,dir);
      if (code < codemin) codemin = code;
    }

  return codemin;
}

traintrack::intVec canonical_coding_from_monogon(const traintrack& tt,
                                                 const int mono,
                                                 const int dir)
{
  if (tt.Multigon(mono).edges() > 1)
    {
      std::cerr << "Not an uncusped monogon in traintracks::coding::canonical_coding_from_monogon.\n";
      std::exit(1);
    }

  traintrack::intVec code;
  append_block(code,0,1,tt.Multigon(mono).label(),0,1);

  int pmono, pemono;
  multigon *egmono =
    tt.Multigon(mono).Edge(0,0)->target_multigon(&tt.Multigon(mono),pmono,pemono);

  jlt::vector<int> first_pin(tt.multigons(),-1);
  recursive_canonical_coding(tt,*egmono,pmono,pemono,code,dir,first_pin);

  return code;
}

traintrack::intVec canonical_coding(const traintrack& tt, const int dir)
{
  int mono = 0;
  traintrack::intVec codemin = canonical_coding_from_monogon(tt,mono,dir);

  while (tt.Multigon(++mono).edges() <= 1)
    {
      traintrack::intVec code = canonical_coding_from_monogon(tt,mono,dir);
      if (code < codemin) codemin = code;
    }

  return codemin;
}

} // namespace coding
} // namespace traintracks
