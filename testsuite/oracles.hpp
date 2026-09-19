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

// Oracles for the testsuite: slow, independent implementations kept to
// cross-check the library.

#ifndef TTAUTO_TESTSUITE_ORACLES_HPP
#define TTAUTO_TESTSUITE_ORACLES_HPP

#include <jlt/mathmatrix.hpp>
#include "traintracks/traintrack.hpp"

namespace oracles {

// One-fold transition matrix TM(new, old) obtained the original way: put
// unit weight on each old edge in turn, fold a copy, read the new weights.
// Costs n folds; the library now reads the matrix off the fold record.
inline jlt::mathmatrix<int>
unit_weight_transition_matrix(const traintracks::traintrack& tt0, const int f)
{
  const int n = tt0.edges();
  jlt::mathmatrix<int> TM(n,n,0);
  for (int i = 0; i < n; ++i)
    {
      traintracks::traintrack tt(tt0);
      traintracks::traintrack::dblVec wv(n,0.0);
      wv[i] = 1;
      tt.weights(wv.begin());
      tt.fold(f);
      wv = tt.weights();
      for (int j = 0; j < n; ++j) TM(j,i) = (int)wv[j];
    }
  return TM;
}

} // namespace oracles

#endif // TTAUTO_TESTSUITE_ORACLES_HPP
