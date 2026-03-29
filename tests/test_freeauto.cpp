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
#include "freeauto.hpp"
#include "traintracks_util.hpp"


int main()
{
  using namespace traintracks;

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

  return 0;
}
