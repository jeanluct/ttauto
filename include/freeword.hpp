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

#ifndef FREEWORD_HPP
#define FREEWORD_HPP

namespace ttauto {

class free_elem
{
protected:
  int inv;

public:
  free_elem(int inv_ = 0) : inv(inv_) {}
};

class main_edge : public free_elem
{
  // Main edge is numbered by edge#.
  int eno;

public:
  main_edge(int eno_ = 0, int inv_ = 0) : free_elem(inv_), eno(eno_) {}
};

class inf_edge : public free_elem
{
  // Infinitesimal edge is numbered by multigon# and prong#.
  int mno, pno;

public:
  inf_edge(int mno_ = 0, int pno_ = 0, int inv_ = 0)
    : free_elem(inv_), mno(mno_), pno(pno_) {}
};

} // namespace ttauto

#endif // FREEWORD_HPP
