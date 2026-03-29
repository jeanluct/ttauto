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
#include "freeauto.hpp"
#include "traintracks_util.hpp"
#include "ttmap_labels.hpp"


int main()
{
  using std::cout;
  using std::endl;
  using jlt::operator<<;
  using namespace traintracks;

  // Core freeword/freeauto sanity checks.
  {
    freeword<int> w({1,2,-2,3,4,-4,-3});
    w.reduce();
    assert((w == freeword<int>({1}))); // 2,-2 and 3,4,-4,-3 cancel

    freeword<int> wi = w.inverse();
    assert((wi == freeword<int>({-1}))); // inverse of {1}

    freeauto<int> id(3);
    freeauto<int> a(3);
    a[1] = {2};
    a[2] = {-1,3};
    a[3] = {3};

    // Right identity under composition as currently implemented.
    freeauto<int> aid = a * id;
    assert((aid[1] == a[1]));
    assert((aid[2] == a[2]));
    assert((aid[3] == a[3]));

    // Decode helper sanity checks (transposed convention).
    {
      jlt::mathmatrix<int> I(jlt::identity_matrix<int>(3));
      permplus1_decode d = decode_transposed_permplus1(I);
      assert(d.is_perm);
      assert(d.row2 == -1 && d.col2 == -1);
      assert(d.perm[0] == 0 && d.perm[1] == 1 && d.perm[2] == 2);
      assert(d.perm_inv[0] == 0 && d.perm_inv[1] == 1 && d.perm_inv[2] == 2);
    }

    {
      // A transposed permutation+1 example:
      // identity permutation with +1 at (row2,col2)=(0,1).
      jlt::mathmatrix<int> M(3,3,0);
      M(0,0) = 1;
      M(1,1) = 1;
      M(2,2) = 1;
      M(0,1) = 1;

      permplus1_decode d = decode_transposed_permplus1(M);
      assert(!d.is_perm);
      assert(d.row2 == 0);
      assert(d.col2 == 1);
      assert(d.perm[0] == 0 && d.perm[1] == 1 && d.perm[2] == 2);
      assert(d.perm_inv[0] == 0 && d.perm_inv[1] == 1 && d.perm_inv[2] == 2);
    }
  }

#if 0
  freeword<int> w({-5,-1,1,2,2,-2,1,1,2,1,-2,2,-1});
  cout << w << endl << endl;
  cout << w.inverse() << endl << endl;
  cout << w*w << endl << endl;
  cout << w*-11 << endl;
  cout << -11*w << endl;
  cout << w.reduce() << endl;

  // Define automorphism: 2 main edges, 3 infinitesimal edges.
  // {1,2,3,a,b} = {1,2,3,4,5}
  // Defaults to the identity automorphism.
  cout << "\nIdentity automorphism:\n" << freeauto<int>(5);

  freeauto<int> T1(5);
  T1[1] = {2};
  T1[2] = {1};
  T1[3] = {3};
  T1[4] = {-4};
  T1[5] = {4,-2,5};
  cout << "\nT1:\n" << T1;

  // Compose the two automorphisms.
  cout << endl << "T1*id\n" << T1*freeauto<int>(5) << endl;

  freeauto<int> T2(5);
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

  // The fold->infinitesimal index map should stay within range.
  for (int f = 0; f < ttv[trk].foldings(); ++f)
    {
      int infix = ttv[trk].fold_infinitesimal_index(f);
      assert(infix >= 0);
      assert(infix < ttv[trk].total_prongs());

      multigon* mmc = 0;
      int pc = -1, ec = -1;
      ttv[trk].fold_cusp_location(f,mmc,pc,ec);
      assert(mmc != 0);
      assert(pc >= 0);
      assert(ec >= 0);
    }

  // Make a folding path through the automaton.
  // This is the example in devel/iss03
  folding_path<traintrack> p(ttg,0);
  p.push_back(1); p.push_back(0);

  // Diagnostic scaffold for issue #3:
  // list one-step fold maps from the initial vertex, so a hand-worked
  // example can be matched against exact fold indices and composition order.
  cout << "\nOne-step fold maps from initial vertex:\n";
  for (int f = 0; f < ttg.foldings(0); ++f)
    {
      cout << "f=" << f << " (" << (f % 2 ? "clockwise" : "counterclockwise")
           << ")\n";
      cout << ttg.traintrack_map(0,f) << endl;
    }

  // For each one-step fold from each automaton vertex, verify matrix/map
  // consistency in both non-transposed and transposed conventions.
  for (int v = 0; v < ttg.vertices(); ++v)
    {
      for (int f = 0; f < ttg.foldings(v); ++f)
        {
          const traintrack& ttvf = ttg.traintrack(v);
          const ttmap_labeler labels(ttvf.edges(),ttvf.total_prongs());
          const freeauto<int>& AMstep = ttg.traintrack_map(v,f);
          const mathmatrix_permplus1& PM = ttg.transition_matrix(v,f);
          permplus1_decode dec = decode_fold_map_structure(ttvf,f);
          assert(dec.is_perm == PM.is_perm());

          jlt::mathmatrix<int> TMstep(ttg.transition_matrix(v,f).full());
          jlt::mathmatrix<int> TMfromMap =
            transition_matrix_from_map(ttg.traintrack(v),ttg.traintrack_map(v,f));
          assert(TMstep == TMfromMap);

          jlt::mathmatrix<int> TMstepT(TMstep);
          TMstepT.transpose();
          jlt::mathmatrix<int> TMfromMapT =
            transition_matrix_from_map_transposed(ttg.traintrack(v),ttg.traintrack_map(v,f));
          assert(TMstepT == TMfromMapT);

          // In one-step fold maps, the selected infinitesimal generator should
          // appear once with the fixed global orientation in non-permutation folds.
          int infg = std::abs(ttvf.fold_infinitesimal_generator(f,labels.nmain));
          int ninf = 0, nneg = 0, npos = 0;
          for (int g = 1; g <= labels.nmain; ++g)
            {
              for (auto img : AMstep.get_action(g))
                {
                  if (img == -infg) { ++ninf; ++nneg; }
                  else if (img == infg) { ++ninf; ++npos; }
                }
            }
          if (PM.is_perm())
            {
              assert(ninf == 0);
            }
          else
            {
              assert(ninf == 1);
              assert(npos == 1);
              assert(nneg == 0);
            }
        }
    }

  // Hand-checked mapping example alignment (issue #3):
  // a=1, b=2, infinitesimal generator at folded cusp has index 5.
  // Step 1: f=1 (clockwise): a->-a, b->a 5 b.
  const freeauto<int>& AMf1 = ttg.traintrack_map(0,1);
  assert((AMf1[1] == freeword<int>({-1})));
  assert((AMf1[2] == freeword<int>({1,5,2})));

  // Step 2 candidate: f=0 from the target of step 1.
  // Keep this as diagnostic until we fully lock geometric labeling/orientation.
  const int v_after_f1 = ttg.target_vertex(0,1);
  const freeauto<int>& AMf0_after_f1 = ttg.traintrack_map(v_after_f1,0);
  assert((AMf0_after_f1[1] == freeword<int>({1,5,2})));
  assert((AMf0_after_f1[2] == freeword<int>({-2})));

  cout << "\nTransition matrix (transposed):\n";
  jlt::mathmatrix<int> TMp = p.transition_matrix().transpose();
  TMp.printMatrixForm(cout) << endl;

  cout << "\nTrain track map:\n";
  freeauto<int> AMp = p.traintrack_map();
  cout << AMp << endl;

  // Composition order check:
  // path [1,0] should compose as AM(1) * AM(0), matching folding_path.
  freeauto<int> AMcheck(ttg.traintrack_map(0,1));
  AMcheck *= ttg.traintrack_map(ttg.target_vertex(0,1),0);
  assert(AMp[1] == AMcheck[1]);
  assert(AMp[2] == AMcheck[2]);

  // Hand-checked composed map for step sequence [1,0]
  // under the current infinitesimal orientation convention.
  assert((AMp[1] == freeword<int>({-2,-5,-1})));
  assert((AMp[2] == freeword<int>({1,5,2,5,-2})));

  // Main-edge transition matrix should match map-derived one.
  jlt::mathmatrix<int> TMfromAMp = transition_matrix_from_map_transposed(ttv[trk],AMp);
  assert(TMp == TMfromAMp);

  //  return 0;

  // Secondary scenario: broader consistency checks over a larger automaton.
  n = 4;
  ttVec ttv2 = ttbuild_list(n);
  trk = 1;

  cout << "\nTrain track has " << ttv2[trk].punctures() << " punctures and ";
  cout << ttv2[trk].edges() << " edges\n";
  cout << "\nConstructing automaton...\n";

  ttgraph ttg2(ttv2[trk]);
  cout << "\nFolding automaton has " << ttg2.vertices();
  cout << (ttg2.vertices() > 1 ? " vertices\n" : " vertex\n");

  for (int f = 0; f < ttv2[trk].foldings(); ++f)
    {
      int infix = ttv2[trk].fold_infinitesimal_index(f);
      assert(infix >= 0);
      assert(infix < ttv2[trk].total_prongs());

      multigon* mmc = 0;
      int pc = -1, ec = -1;
      ttv2[trk].fold_cusp_location(f,mmc,pc,ec);
      assert(mmc != 0);
      assert(pc >= 0);
      assert(ec >= 0);
    }

  for (int v = 0; v < ttg2.vertices(); ++v)
    {
      for (int f = 0; f < ttg2.foldings(v); ++f)
        {
          const traintrack& ttvf = ttg2.traintrack(v);
          const ttmap_labeler labels(ttvf.edges(),ttvf.total_prongs());
          const freeauto<int>& AMstep = ttg2.traintrack_map(v,f);
          const mathmatrix_permplus1& PM = ttg2.transition_matrix(v,f);
          permplus1_decode dec = decode_fold_map_structure(ttvf,f);
          assert(dec.is_perm == PM.is_perm());

          jlt::mathmatrix<int> TMstep(ttg2.transition_matrix(v,f).full());
          jlt::mathmatrix<int> TMfromMap =
            transition_matrix_from_map(ttg2.traintrack(v),ttg2.traintrack_map(v,f));
          assert(TMstep == TMfromMap);

          jlt::mathmatrix<int> TMstepT(TMstep);
          TMstepT.transpose();
          jlt::mathmatrix<int> TMfromMapT =
            transition_matrix_from_map_transposed(ttg2.traintrack(v),ttg2.traintrack_map(v,f));
          assert(TMstepT == TMfromMapT);

          int infg = std::abs(ttvf.fold_infinitesimal_generator(f,labels.nmain));
          int ninf = 0, nneg = 0, npos = 0;
          for (int g = 1; g <= labels.nmain; ++g)
            {
              for (auto img : AMstep.get_action(g))
                {
                  if (img == -infg) { ++ninf; ++nneg; }
                  else if (img == infg) { ++ninf; ++npos; }
                }
            }
          if (PM.is_perm())
            {
              assert(ninf == 0);
            }
          else
            {
              assert(ninf == 1);
              assert(npos == 1);
              assert(nneg == 0);
            }
        }
    }

  // Make a folding path through the automaton.
  folding_path<traintrack> p2(ttg2,0);
  p2.push_back(0); p2.push_back(1); p2.push_back(0);

  cout << "\nTransition matrix (transposed):\n";
  jlt::mathmatrix<int> TMp2 = p2.transition_matrix().transpose();
  TMp2.printMatrixForm(cout) << endl;

  cout << "\nTrain track map:\n";
  freeauto<int> AMp2 = p2.traintrack_map();
  cout << AMp2 << endl;

  jlt::mathmatrix<int> TMfromAMp2 = transition_matrix_from_map_transposed(ttv2[trk],AMp2);
  assert(TMp2 == TMfromAMp2);
}
