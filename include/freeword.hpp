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

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <initializer_list>
#include <iostream>
#include <list>

namespace traintracks {

// Remove from container adjacent elements satisfying a comparison function.
// Similar to unique (http://www.cplusplus.com/reference/algorithm/unique/)
template <class ForwardIterator, class BinaryPredicate>
ForwardIterator adjacent_remove_if(ForwardIterator first,
				   ForwardIterator last,
				   BinaryPredicate pred)
{
  if (first == last) return last;

  ForwardIterator result = first;
  while (++first != last)
    {
      if (!pred(*result,*first))
	{
	  // The elements are not equal, so copy to result and move on.
	  *(++result) = std::move(*first);
	}
      else
	{
	  // The elements are equal.
	  // Point the the next element (after the elements to be cancelled).
	  // If we're at the end of the container, we're done.
	  if (++first == last) return result;
	  // Copy the next element after the cancelled ones.
	  *result = std::move(*first);
	}
    }
  return ++result;
}


// A default inverse function, appropriate for ints and doubles.
template<class T>
inline T inverse(const T& i) { return -i; }


// Forward declarations.
template<class T> class freeword;
template<class T>
freeword<T> operator*(const freeword<T>&, const T&);
template<class T>
freeword<T> operator*(const T&, const freeword<T>&);
template<class T>
freeword<T> operator*(const freeword<T>&, const freeword<T>&);


template<class T>
class freeword : public std::list<T>
{
  int ngen;

public:
  freeword(int ngen_ = 1) : ngen(ngen_) {}

  // Forward initializer list to base class.
  freeword(std::initializer_list<T> l, int ngen_ = 0) : std::list<T>(l)
  {
    /* This is specialized to ints for now. */
    int ngen_found = 0;
    for (auto i : *this)
      {
	ngen_found = std::max(std::abs(i),ngen_found);
      }
    if (ngen_ && ngen_found > ngen_)
      {
	std::cerr << "ngen is too small.\n";
	std::exit(1);
      }
    ngen = std::max(ngen_found,ngen_);
  }

  /* Make a copy constructor and operator= that check gen. */

  // Reduce a word by cancelling adjacent inverses.
  freeword<T>& reduce();

  // Return the inverse of the free word.
  freeword<T> inverse() const
  {
    freeword<T> iw(ngen);
    for (auto i = this->rbegin(); i != this->rend(); ++i) iw.push_back(-(*i));
    return iw;
  }

  // Multiplication just appends and prepends elements to the word.
  friend freeword<T> operator*<>(const freeword<T>&, const T&);
  friend freeword<T> operator*<>(const T&, const freeword<T>&);
  friend freeword<T> operator*<>(const freeword<T>&, const freeword<T>&);
};

// Reduce a word by cancelling adjacent inverses.
template<class T>
inline freeword<T>& freeword<T>::reduce()
{
  struct IsInv {
    bool operator() (const T& a, const T& b)
    { return (a == traintracks::inverse(b)); }
    // Make sure to specify inverse() function, not method.
  } isinv;

  while (true)
    {
      auto i = adjacent_remove_if(this->begin(),this->end(),isinv);
      if (i == this->end()) break;
      this->erase(i,this->end());
    }

  return *this;
}

template<class T>
inline freeword<T> operator*(const freeword<T>& ww, const T& ee)
{
  freeword<T> ww2 = ww;
  ww2.push_back(ee);
  return ww2;
}

template<class T>
inline freeword<T> operator*(const T& ee, const freeword<T>& ww)
{
  freeword<T> ww2 = ww;
  ww2.push_front(ee);
  return ww2;
}

template<class T>
inline freeword<T> operator*(const freeword<T>& w1, const freeword<T>& w2)
{
  freeword<T> w12 = w1;
  w12.insert(w12.end(),w2.begin(),w2.end());
  return w12;
}


} // namespace traintracks

#endif // FREEWORD_HPP
