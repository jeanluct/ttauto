// <LICENSE
// LICENSE>

#include <iostream>
#include "traintrack.hpp"
#include "ttfoldgraph.hpp"

using namespace traintracks;

int main()
{
  for (int n = 3; n <= 9; ++n)
    {
      // Train track with n monogons and an (n-2)-gon on the boundary.
      traintrack tt(n);

      // Make train track folding automaton graph.
      ttfoldgraph<traintrack> ttg(tt);

      // Number of vertices in full automaton.
      std::cout << n << "\t" << ttg.vertices();
      // Number of vertices in main automaton (that supports pA's).
      std::cout << "\t" << subgraphs(ttg).front().vertices() << std::endl;
    }
}
