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

// ttbraid: the braid of a closed folding path.
//
//   ttbraid <punctures> <stratum> <vertex> <branch> <branch> ...
//
// Strata and vertices are 1-based, as everywhere the program talks to the
// user, and branches are 0-based as in folding_path.  With no branches the
// program lists the strata of that many punctures.
//
// A closed folding path defines a homeomorphism of the punctured disc: do
// the folds, then identify the final track with the initial one by the
// canonical numbering.  Reading that off as a braid needs the punctures to
// have positions along the real axis, which is what a proper embedding
// gives; see traintracks/embedding.hpp and traintracks/braid.hpp.
//
// The braid printed is checked against the path it came from: its growth,
// computed from the action on Dynnikov coordinates, must equal the Perron
// root of the path's transition matrix.  That check only means anything
// when the path is pseudo-Anosov, so it is reported and not enforced.

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include "traintracks/braid.hpp"
#include "traintracks/build.hpp"
#include "traintracks/embedding.hpp"
#include "traintracks/traintrack.hpp"
#include "ttauto/folding_path.hpp"
#include "ttauto/path_braid.hpp"
#include "ttauto/ttfoldgraph.hpp"

using traintracks::traintrack;

int main(int argc, char* argv[])
{
  using std::cout;
  using std::cerr;

  if (argc < 3)
    {
      cerr << "usage: " << argv[0]
           << " <punctures> <stratum> [<vertex> <branch> ...]\n";
      return 1;
    }

  const int n = std::atoi(argv[1]);
  if (n < 3)
    {
      cerr << "Need at least three punctures.\n";
      return 1;
    }

  jlt::vector<traintrack> ttv = traintracks::build_traintrack_list(n);
  const int stratum = std::atoi(argv[2]);
  if (stratum < 1 || stratum > (int)ttv.size())
    {
      cerr << "Stratum must be between 1 and " << ttv.size() << " for "
           << n << " punctures.\n";
      for (int s = 0; s < (int)ttv.size(); ++s)
        {
          cerr << "  " << s+1 << "  ";
          ttv[s].print_singularity_data(cerr) << "\n";
        }
      return 1;
    }

  ttauto::ttfoldgraph<traintrack> ttg(ttv[stratum-1]);
  cout << "Stratum " << stratum << ": ";
  ttv[stratum-1].print_singularity_data(cout);
  cout << "\nAutomaton has " << ttg.vertices() << " vertices.\n";

  if (argc < 5) return 0;

  const int v0 = std::atoi(argv[3]) - 1;
  if (v0 < 0 || v0 >= ttg.vertices())
    {
      cerr << "Vertex must be between 1 and " << ttg.vertices() << ".\n";
      return 1;
    }

  ttauto::folding_path<traintrack> p(ttg,v0);
  for (int i = 4; i < argc; ++i)
    {
      const int br = std::atoi(argv[i]);
      if (br < 0 || br >= ttg.foldings(p.final_vertex()))
        {
          cerr << "Branch " << br << " does not exist at vertex "
               << p.final_vertex()+1 << ", which has "
               << ttg.foldings(p.final_vertex()) << ".\n";
          return 1;
        }
      p.push_back(br);
    }

  cout << "\nPath visits vertices";
  for (int i = 0; i < (int)p.vertices().size(); ++i)
    cout << " " << p.vertices()[i]+1;
  cout << "\n";

  if (!p.closed())
    {
      cerr << "\nThe path is not closed, so it defines no braid.\n";
      return 1;
    }

  const double lambda = ttauto::detail::perron_root(p.transition_matrix());
  cout << "Growth of the path                 " << lambda << "\n";
  cout << "Bestvina-Handel gates connected    "
       << (p.gates().connected ? "yes" : "no") << "\n";

  // The proper embedding of the initial track, which is where the
  // positions of the punctures come from.
  const traintracks::tt_embedding emb
    = traintracks::outer_embedding(ttg.traintrack(v0).numbering());
  cout << "Punctures of the initial track, left to right, by prong number:";
  for (int i = 0; i < emb.npunctures(); ++i)
    cout << " " << emb.puncture_order[i];
  cout << "\n";

  bool verified = false;
  const traintracks::braidword b = ttauto::folding_path_braid(p,&verified);

  cout << "\nBraid on " << b.strings() << " strings, length " << b.length()
       << ", exponent sum " << b.exponent_sum() << ":\n  ";
  b.print(cout) << "\n";
  cout << "Mathematica form: ";
  b.printMathematicaForm(cout) << "\n";
  cout << "Permutation:";
  const std::vector<int> perm = b.permutation();
  for (int i = 0; i < (int)perm.size(); ++i) cout << " " << perm[i];
  cout << "\nGrowth of the braid                " << b.growth() << "\n";
  cout << "Agrees with the path               "
       << (verified ? "yes" : "no") << "\n";
  if (!verified)
    cout << "  (A path that is not pseudo-Anosov is not expected to agree."
         << "  For one\n   that is, see the limitation in"
         << " doc/CODE_STRUCTURE.md.)\n";

  return 0;
}
