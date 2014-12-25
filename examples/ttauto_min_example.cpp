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

//
// ttauto_min_example
//
//  Reproduce the results of Song, Ko & Los (2002) and Ham & Song
//  (2007) on the minimum dilation of braids on the discs with 3 to 5
//  punctures.  The dilatation is checked for each stratum on a given
//  punctured disc.  For 6 and 7 punctures, attempt to directly test
//  the results of Lanneau & Thiffeault (2009), but for several strata
//  this is likely to take too long on most computers.
//
//  This example program doesn't use more powerful techniques such as
//  decomposing the automata into invariant subgraphs or looking for
//  "bad words".
//
//  References
//
//  W. T. Song, K. H. Ko, and J. E. Los, "Entropies of braids,"
//  J. Knot Th. Ramifications 11 (2002), 647-666.
//
//  J.-Y. Ham and W. T. Song, "The minimum dilatation of
//  pseudo-Anosov 5-braids," Experiment. Math. 16 (2007), 167-179.
//
//  E. Lanneau and J.-L. Thiffeault, "On the minimum dilatation of
//  braids on punctured discs," preprint (2009).
//


int main()
{
  using std::cout;
  using std::cerr;
  using std::endl;
  using std::vector;
  using traintracks::traintrack;

  // Number of punctures.
  const int n = 5;

  // Build list of initial train tracks.
  vector<traintrack> ttv = traintracks::ttbuild_list(n);

  // Minimum dilatations for each stratum (see Ham & Song 2007).
  vector<double> dilmax(ttv.size());
  if (n == 3)
    {
      // Classical
      dilmax[0] = 2.62;
    }
  else if (n == 4)
    {
      // Song, Ko & Los (2002)
      dilmax[0] = 2.62; dilmax[1] = 2.30;
    }
  else if (n == 5)
    {
      // Ham & Song (2007)
      dilmax[0] = 1.73; dilmax[1] = 1.73; dilmax[2] = 2.16; dilmax[3] = 2.02;
    }
  else if (n == 6)
    {
      // Lanneau & Thiffeault (2009)
      cerr << "\n**** This will run a while! ****\n";
      dilmax[0] = 1.89; dilmax[1] = 1.84; dilmax[2] = 1.89; dilmax[3] = 2.09;
      dilmax[4] = 2.09; dilmax[5] = 1.89; dilmax[6] = 2.18;
    }
  else if (n == 7)
    {
      // Lanneau & Thiffeault (2009)
      cerr << "\n**** This will run a while! ****\n";
      dilmax[0] = 1.56; dilmax[1] = 1.47; dilmax[2] = 1.47; dilmax[3] = 1.56;
      dilmax[4] = 2.05; dilmax[5] = 1.62; dilmax[6] = 2.48; dilmax[7] = 1.81;
      dilmax[8] = 1.76; dilmax[9] = 1.62; dilmax[10] = 2.05; dilmax[11] = 2.22;
    }
  else
    {
      cerr << "\n**** Too many punctures! ****\n\n";
      exit(-1);
    }

  // Loop over strata and find pA's with dilatation <= dilmax.
  for (int trk = 0; trk < (int)ttv.size(); ++trk)
    {
      cout << "\n===================================================";
      cout << "\nSTRATUM " << trk+1 << ":  ";
      ttv[trk].print_singularity_data(cout) << endl;

      // Create the train track graph (automaton).
      traintracks::ttfoldgraph<traintrack> ttg(ttv[trk]);

      cout << "\nSearching automaton graph with " << ttg.vertices();
      cout << " vertices for pseudo-Anosovs...\n";
      traintracks::ttauto<traintrack> tta(ttg);
      tta.max_dilatation(dilmax[trk]).badword_length(0).check_norms();
      tta.search();

      cout << "\n\nLowest dilatation found on stratum is root of ";
      cout << tta.pA_list().begin()->second.polynomial() << endl;
    }
}
