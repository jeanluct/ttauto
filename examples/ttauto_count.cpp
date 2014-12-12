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
