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

// What labelling a multigon does to the automaton and to the
// pseudo-Anosovs it finds.  This is a paired comparison: the same search
// is run with and without a label on one of the two trigons, and the
// point is the difference between them, not either count on its own.

#include <iostream>
#include <jlt/vector.hpp>
#include "traintracks/traintrack.hpp"
#include "traintracks/build.hpp"
#include "ttauto/ttfoldgraph.hpp"
#include "ttauto/ttauto.hpp"


int main()
{
  using std::cout;
  using std::endl;
  using traintracks::traintrack;

  using ttgraph = ttauto::ttfoldgraph<traintrack>;
  using ttVec = jlt::vector<traintrack>;

  const int n = 5;		// Number of punctures.
  int trk = 3;			// Initial train track to search.
  double dilmax = 3;		// Max dilatation to look for.

  ttVec ttv = traintracks::build_traintrack_list(n);

  //
  // 5 punctures:
  //
  // Track #3 has 2 x 3-prongs.
  // Give one of them (multigon 6) the label "1" instead of default 0:
  const int label = 1;
  ttv[trk].set_label(6,label);
  ttv[trk].print(cout) << endl;

  // Labelling one trigon eliminates some of the pseudo-Anosovs: those
  // whose folding path permutes the two trigons no longer close up, so
  // the classes that survive are a subset of the ones found without the
  // label.  To see the other side of the comparison, recompile and run
  // with set_label(6,0) above, rather than (6,1).
  //
  // The minimum dilatation is the same either way, 2.01536 at path
  // length 4, so it is the control here rather than the finding.
  //
  // Labelling also doubles the automaton, 9 vertices to 18: each
  // unlabelled track corresponds to two labelled ones, according to
  // which of its two trigons carries the label.
  //
  // Measured 2026-09-25, with dilmax = 3, badword_length 0 and the
  // path-length cap below: 10 classes without the label, 6 with, so 4 of
  // the 10 permute the trigons.  Quote the configuration whenever you
  // quote the counts -- they move with the window, the cap and the
  // pruning, which is how the previous version of this comment came to
  // claim 8 and 2.  Those numbers could not be reproduced under any
  // setting tried; see devel/iss016/labelled_automaton.md.

  cout << "Train track has " << ttv[trk].punctures() << " punctures and ";
  cout << ttv[trk].edges() << " edges\n";

  // Make a list of train track graphs.  The main graph is the first
  // element.
  std::list<ttgraph> ttg(ttauto::subgraphs(ttgraph(ttv[trk])));

  cout << "\nFolding subgraphs from initial train track: \n";
  ttauto::print_subgraphs(ttg);

  for (auto i = ttg.begin(); i != ttg.end(); ++i)
    {
      int fg = std::distance(ttg.begin(),i);
      cout << "\n\nFOLDING GRAPH " << fg;
      cout << " with " << i->vertices();
      cout << (i->vertices() == 1 ? " vertex\n\n" : " vertices\n\n");

      ttauto::ttauto<traintrack> tta(*i);
      tta.max_dilatation(dilmax).check_norms();
      // check_norms() derives its own path-length bound from the
      // dilatation window, and at dilmax = 3 on this stratum that is 734
      // folds, which is why this example used to run for about a hundred
      // minutes.  The classes found are unchanged at every cap from 10 to
      // 24, so 12 is used here; set it after check_norms(), which would
      // otherwise overwrite it.  This is a measured cap, not a proved
      // one: tta.path_length_exceeded() is nonzero at 12, so the search
      // is complete for the window only up to that length.
      tta.max_pathlength(12);
      tta.search();
    }
}
