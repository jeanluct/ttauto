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

// Set traintrack::label_multiprongs = true in traintrack.hpp.

#include <iostream>
#include <jlt/vector.hpp>
#include "traintrack.hpp"
#include "ttfoldgraph.hpp"
#include "ttauto.hpp"

using namespace traintracks;

int main()
{
  using std::cout;
  using std::endl;

  typedef ttfoldgraph<traintrack>			ttgraph;
  typedef std::list<ttgraph>::iterator			ttgit;
  typedef std::list<ttgraph>::const_iterator		cttgit;
  typedef jlt::vector<traintrack>			ttVec;
  typedef ttVec::const_iterator				ttVeccit;

  const int n = 5;		// Number of punctures.
  int trk = 3;			// Initial train track to search.
  double dilmax = 3;		// Max dilatation to look for.

  ttVec ttv = ttbuild_list(n);

  //
  // 5 punctures:
  //
  // Track #3 has 2 x 3-prongs.
  // Give one of them (multigon 6) the label "1" instead of default 0:
  const int label = 1;
  ttv[trk].set_label(6,label);
  ttv[trk].print(cout) << endl;

  // Find that the minimum dilatation is 2.015, as when the 3-prongs
  // are identical, but there are far fewer pA's here.  To see this,
  // recompile and run with set_label(6,0) above, rather than (6,1).
  //
  // Find 8 pA's with dilatation <= 3 for identical 3-prongs, but only
  // two for distinct 3-prongs.  This means that the six others
  // permute the singularities.

  cout << "Train track has " << ttv[trk].punctures() << " punctures and ";
  cout << ttv[trk].edges() << " edges\n";

  // Make a list of train track graphs.  The main graph is the first
  // element.
  std::list<ttgraph> ttg(subgraphs(ttgraph(ttv[trk])));

  cout << "\nFolding subgraphs from initial train track: \n";
  print_subgraphs(ttg);

  for (cttgit i = ttg.begin(); i != ttg.end(); ++i)
    {
      int fg = std::distance(cttgit(ttg.begin()),i);
      cout << "\n\nFOLDING GRAPH " << fg;
      cout << " with " << i->vertices();
      cout << (i->vertices() == 1 ? " vertex\n\n" : " vertices\n\n");

      ttauto<traintrack> tta(*i);
      tta.max_dilatation(dilmax).check_norms();
      tta.search();
    }
}
