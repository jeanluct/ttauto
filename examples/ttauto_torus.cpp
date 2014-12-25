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
#include <vector>
#include "traintrack.hpp"
#include "ttfoldgraph.hpp"
#include "ttauto.hpp"


int main()
{
  using std::cout;
  using std::endl;
  using std::vector;
  using ttauto::traintrack;

  // Number of punctures.
  const int n = 3;

  // Build list of initial train tracks.
  vector<traintrack> ttv = ttauto::ttbuild_list(n);

  // Minimum dilatations for each stratum (see Ham & Song 2007).
  vector<double> dilmax(ttv.size());
  dilmax[0] = 50;

  // Loop over strata and find pA's with dilatation <= dilmax.
  for (int trk = 0; trk < (int)ttv.size(); ++trk)
    {
      cout << "\n===================================================";
      cout << "\nSTRATUM " << trk+1 << ":  ";
      ttv[trk].print_singularity_data(cout) << endl;

      // Create the train track graph (automaton).
      ttauto::ttfoldgraph<traintrack> ttg(ttv[trk]);

      cout << "\nSearching automaton graph with " << ttg.vertices();
      cout << " vertices for pseudo-Anosovs...\n";
      ttauto::ttauto<traintrack> tta(ttg);
      tta.max_dilatation(dilmax[trk]).badword_length(0).check_norms();
      tta.search();

      cout << "\n\nLowest dilatation found on stratum is root of ";
      cout << tta.pA_list().begin()->second.polynomial() << endl;
    }
}
