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
#include <cassert>
#include <list>
#include <jlt/stlio.hpp>
#include "traintrack.hpp"
#include "ttfoldgraph.hpp"
#include "folding_path.hpp"
#include "freeword.hpp"
#include "traintracks_util.hpp"


int main()
{
  using std::cout;
  using std::endl;
  using jlt::operator<<;
  using namespace traintracks;

  // Core free_word/free_auto sanity checks.
  {
    free_word<int> w({1,2,-2,3,4,-4,-3});
    w.reduce();
    assert((w == free_word<int>({1}))); // 2,-2 and 3,4,-4,-3 cancel

    free_word<int> wi = w.inverse();
    assert((wi == free_word<int>({-1}))); // inverse of {1}

    free_auto<int> id(3);
    free_auto<int> a(3);
    a[1] = {2};
    a[2] = {-1,3};
    a[3] = {3};

    // Right identity under composition as currently implemented.
    free_auto<int> aid = a * id;
    assert((aid[1] == a[1]));
    assert((aid[2] == a[2]));
    assert((aid[3] == a[3]));
  }

#if 0
  free_word<int> w({-5,-1,1,2,2,-2,1,1,2,1,-2,2,-1});
  cout << w << endl << endl;
  cout << w.inverse() << endl << endl;
  cout << w*w << endl << endl;
  cout << w*-11 << endl;
  cout << -11*w << endl;
  cout << w.reduce() << endl;

  // Define automorphism: 2 main edges, 3 infinitesimal edges.
  // {1,2,3,a,b} = {1,2,3,4,5}
  // Defaults to the identity automorphism.
  cout << "\nIdentity automorphism:\n" << free_auto<int>(5);

  free_auto<int> T1(5);
  T1[1] = {2};
  T1[2] = {1};
  T1[3] = {3};
  T1[4] = {-4};
  T1[5] = {4,-2,5};
  cout << "\nT1:\n" << T1;

  // Compose the two automorphisms.
  cout << endl << "T1*id\n" << T1*free_auto<int>(5) << endl;

  free_auto<int> T2(5);
  T2[1] = {1};
  T2[2] = {3};
  T2[3] = {2};
  T2[4] = {4,-2,5};
  T2[5] = {-5};
  cout << "\nT2:\n" << T2;
  // Compose the two automorphisms.
  cout << endl << "T1*T2:\n" << T1*T2 << endl;
#endif
  typedef ttfoldgraph<traintrack>			ttgraph;
  typedef jlt::vector<traintrack>			ttVec;

  int n = 3;
  ttVec ttv = ttbuild_list(n);
  int trk = 0;

  cout << "\nTrain track has " << ttv[trk].punctures() << " punctures and ";
  cout << ttv[trk].edges() << " edges\n";
  cout << "\nConstructing automaton...\n";

  ttgraph ttg(ttv[trk]);
  cout << "\nFolding automaton has " << ttg.vertices();
  cout << (ttg.vertices() > 1 ? " vertices\n" : " vertex\n");

  // Make a folding path through the automaton.
  // This is the example in devel/iss03
  folding_path<traintrack> p(ttg,0);
  p.push_back(1); p.push_back(0);

  cout << "\nTransition matrix (transposed):\n";
  jlt::mathmatrix<int> TMp = p.transition_matrix().transpose();
  TMp.printMatrixForm(cout) << endl;

  cout << "\nTrain track map:\n";
  free_auto<int> AMp = p.traintrack_map();
  cout << AMp << endl;

  // Main-edge transition matrix should match map-derived one.
  jlt::mathmatrix<int> TMfromAMp = transition_matrix_from_map_transposed(ttv[trk],AMp);
  assert(TMp == TMfromAMp);

  //  return 0;

  /* Something's weird in the rest... */
  /* Doesn't reproduce the correct automaton? */
  n = 4;
  ttVec ttv2 = ttbuild_list(n);
  trk = 1;

  cout << "\nTrain track has " << ttv2[trk].punctures() << " punctures and ";
  cout << ttv2[trk].edges() << " edges\n";
  cout << "\nConstructing automaton...\n";

  ttgraph ttg2(ttv2[trk]);
  cout << "\nFolding automaton has " << ttg2.vertices();
  cout << (ttg2.vertices() > 1 ? " vertices\n" : " vertex\n");

  // Make a folding path through the automaton.
  folding_path<traintrack> p2(ttg2,0);
  p2.push_back(0); p2.push_back(1); p2.push_back(0);

  cout << "\nTransition matrix (transposed):\n";
  jlt::mathmatrix<int> TMp2 = p2.transition_matrix().transpose();
  TMp2.printMatrixForm(cout) << endl;

  cout << "\nTrain track map:\n";
  free_auto<int> AMp2 = p2.traintrack_map();
  cout << AMp2 << endl;

  jlt::mathmatrix<int> TMfromAMp2 = transition_matrix_from_map_transposed(ttv2[trk],AMp2);
  assert(TMp2 == TMfromAMp2);
}
