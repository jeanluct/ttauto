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

#include <cassert>
#include <jlt/mathmatrix.hpp>
#include "traintracks/mathmatrix_permplus1.hpp"


static void check_case(const jlt::mathmatrix<int>& Mpm)
{
  // Verify that a matrix accepted by mathmatrix_permplus1 can be losslessly
  // represented and expanded back to dense form.
  const int n = Mpm.dim();
  traintracks::mathmatrix_permplus1 pm(Mpm);

  assert(pm.full() == Mpm);
  assert(pm.order() > 0);

  jlt::mathmatrix<int> M(jlt::identity_matrix<int>(n));
  for (int i = 0; i < n; ++i)
    {
      for (int j = 0; j < n; ++j)
        {
          M(i,j) = i - 3 + n*j + i*j;
        }
    }

  // Verify mixed products agree with dense multiplication on both sides.
  assert(Mpm * M == pm * M);
  assert(M * Mpm == M * pm);
}


int main()
{
  // Case 1: plain permutation matrix.
  const int n = 6;

  jlt::mathmatrix<int> Mperm(n,n);
  Mperm(0,4) = 1;
  Mperm(1,0) = 1;
  Mperm(2,1) = 1;
  Mperm(3,2) = 1;
  Mperm(4,3) = 1;
  Mperm(5,5) = 1;

  check_case(Mperm);

  // Case 2: permutation-plus-one matrix.
  jlt::mathmatrix<int> Mplus1(Mperm);
  Mplus1(5,4) = 1;

  check_case(Mplus1);

  return 0;
}
