// <LICENSE
//   ttauto: a C++ library for building train track automata
//
//   https://github.com/jeanluct/ttauto
//
//   Copyright (C) 2010-2014  Jean-Luc Thiffeault   <jeanluc@math.wisc.edu>
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
#include <list>
#include <jlt/matrix.hpp>
#ifdef DOUBLEHUMP_TRACK
#include "doublehump.hpp"
#else
#include "traintrack.hpp"
#include "ttfoldgraph.hpp"
#include "ttauto.hpp"
#endif
#include "traintracks_util.hpp"
#include "folding_path.hpp"
#include "ttauto.hpp"
#include "badwords.hpp"


int main()
{
  using std::cout;
  using std::endl;
  using ttauto::traintrack;
  using ttauto::folding_path;

#ifdef DOUBLEHUMP_TRACK
  typedef doublehump TrTr;
  doublehump tt(1,3,1, 1,1);
#if 1
  // Append a 3-prong to the 1st free prong of b branch, then another 3-prong.
  tt.b.append(1,3);//->append(1,3);//->append(1,3)->append(1,3)->append(1,3);
  // Append a 3-prong to the 2nd free prong of b branch, then another 3-prong.
  //tt.b.append(2,3);//->append(1,3);
#endif

  tt.check();
#else
  typedef traintrack TrTr;
#if 1
  const int n = 6;
  traintrack tt(n,3);
  int maxplen = 5;
#else
  const int n = 5; // Ham&Song
  traintrack tt(n,3,3);
  int maxplen = 6;
#endif
#endif

  typedef std::list<folding_path<TrTr> >::const_iterator path_it;

  cout << "Train track has " << tt.punctures() << " punctures and ";
  cout << tt.edges() << " edges\n";

  ttauto::ttfoldgraph<TrTr> ttg(tt);

  cout << "TT folding graph has " << ttg.vertices() << " vertices\n";

  jlt::matrix<std::list<folding_path<TrTr> > > pbad(badwords(ttg,maxplen));

  cout << "Bad words:\n";

  for (int plen = 0; plen < maxplen; ++ plen)
    {
      for (int v = 0; v < ttg.vertices(); ++v)
	{
	  for (path_it i = pbad(v,plen).begin(); i != pbad(v,plen).end(); ++i)
	    {
	      cout << "Vertex " << v << " " << "\tfolding path " << *i;
	      cout << "\t(vertices " << i->vertices() << ")\n";
	    }
	}
    }

  return 0;
}
