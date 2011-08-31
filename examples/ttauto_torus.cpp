#include <iostream>
#include <vector>
#include "traintrack.hpp"
#include "ttfoldgraph.hpp"
#include "ttauto.hpp"

using namespace traintracks;
using namespace std;

int main()
{
  // Number of punctures.
  const int n = 3;

  // Build list of initial train tracks.
  vector<traintrack> ttv = ttbuild_list(n);

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
      ttfoldgraph<traintrack> ttg(ttv[trk]);

      cout << "\nSearching automaton graph with " << ttg.vertices();
      cout << " vertices for pseudo-Anosovs...\n";
      ttauto<traintrack> tta(ttg);
      tta.max_dilatation(dilmax[trk]).badword_length(0).check_norms();
      tta.search();

      cout << "\n\nLowest dilatation found on stratum is root of ";
      cout << tta.pA_list().begin()->second.polynomial() << endl;
    }
}
