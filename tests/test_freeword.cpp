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
#include <jlt/stlio.hpp>
#include <list>
#include "freeword.hpp"


int main()
{
  using std::cout;
  using std::endl;
  using jlt::operator<<;
  using namespace ttauto;

  free_word<int> w = {-1,1,2,2,-2,1,1,2,1,-2,2,-1};
  cout << w << endl << endl;
  cout << w*-11 << endl;
  cout << -11*w << endl;
  cout << w.reduce() << endl;
}
