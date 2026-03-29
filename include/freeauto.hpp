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

#ifndef FREEAUTO_HPP
#define FREEAUTO_HPP

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>

#include "freeword.hpp"

namespace traintracks {

template<class T>
class freeauto : std::vector<freeword<T> >
{
public:
  freeauto(const int ngen_) : std::vector<freeword<T> >(ngen_+1)
  {
    // Initialize to the identity automorphism.
    for (size_t i = 1; i <= numgens(); ++i)
      std::vector<freeword<T> >::operator[](i).push_back(i);
  }

  // Construct from base class.
  freeauto(const std::vector<freeword<T> >& a) : std::vector<freeword<T> >(a)
  {
  }

  // Number of generators.  Note that the vector size has to be one larger.
  size_t numgens() const { return (this->size()-1); }

  // Action of the automorphism on a single generator a.
  // Use this for assignment.
  freeword<T>& operator[](const T a)
  {
    if (a == 0)
      {
	std::cerr << "Cannot specify action on identity.\n";
	std::exit(1);
      }
    if (a < 0)
      {
	std::cerr << "Assign to positive generators only.\n";
	std::exit(1);
      }
    // Only store the action on positive generators.
    return std::vector<freeword<T> >::operator[](a);
  }

  // Const version of operator[].
  // Must still return a reference, so that pointer operations like
  // begin() are still valid.
  const freeword<T>& operator[](const T a) const
  {
    if (a == 0)
      {
	std::cerr << "Cannot specify action on identity.\n";
	std::exit(1);
      }
    if (a < 0)
      {
	std::cerr << "Assign to positive generators only.\n";
	std::exit(1);
      }
    // Only store the action on positive generators.
    return std::vector<freeword<T> >::operator[](a);
  }

  // get_action is like operator[], but returns by value.
  // It can thus return a freeword for a<0.
  const freeword<T> get_action(const T a) const
  {
    if (a == 0)
      {
	std::cerr << "Cannot specify action on identity.\n";
	std::exit(1);
      }
    if (a > 0)
      return std::vector<freeword<T> >::operator[](a);
    else // For negative a, call freeword<T>::inverse().
      return (std::vector<freeword<T> >::operator[](-a)).inverse();
  }

  // Compose automorphisms: this = this followed by a.
  freeauto<T>& operator*=(const freeauto<T>& a)
  {
    if (numgens() != a.numgens())
      {
	std::cerr << "Can only compose freeauto objects with ";
	std::cerr << "the same number of generators.\n";
	std::exit(1);
      }
    freeauto<T> res(std::vector<freeword<T> >(numgens()+1));
    for (size_t i = 1; i <= numgens(); ++i)
      {
	for (auto j : this->operator[](i))
	  {
	    res[i] = res[i]*a.get_action(j);
	  }
      }
    return (*this = res);
  }
};

template<class T>
inline freeauto<T> operator*(const freeauto<T>& a, const freeauto<T>& b)
{
  return (freeauto<T>(a) *= b);
}

template<class T>
std::ostream& operator<<(std::ostream& strm, const freeauto<T>& a)
{
  constexpr int wid = 5;

  for (size_t i = 1; i <= a.numgens(); ++i)
    {
      strm << std::setw(wid);
      strm << i << " -> ";
      for (auto j : a[i]) strm << std::setw(wid) << j;
      strm << std::endl;
    }
  return strm;
}

} // namespace traintracks

#endif // FREEAUTO_HPP
