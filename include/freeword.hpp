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

#include <iostream>

namespace ttauto {

// Remove equal adjacent elements from container.
// This can be used to cancel inverses by using the right comparison function.
// Similar to unique (http://www.cplusplus.com/reference/algorithm/unique/)
template <class ForwardIterator, class BinaryPredicate>
ForwardIterator adjacent_remove(ForwardIterator first, ForwardIterator last,
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
template<class T> class free_word;
template<class T> free_word<T> operator*(const free_word<T>&, const T&);
template<class T> free_word<T> operator*(const T&, const free_word<T>&);


template<class T>
class free_word
{
  std::list<T> word;

public:
  // Perfect-forward arguments to base constructor.
  template<typename... Args>
  free_word(Args&&... _args) : word(std::forward<Args>(_args)...) {}

  // Forward initializer list to base class.
  free_word(std::initializer_list<T> l) : word(l) {}

  // Reduce a word by cancelling adjacent inverses.
  free_word<T>& reduce();

  // Multiplication just appends and prepends elements to the word.
  friend free_word<T> operator*<>(const free_word<T>&, const T&);
  friend free_word<T> operator*<>(const T&, const free_word<T>&);

  std::ostream& print(std::ostream& strm) const
  {
    copy(word.begin(), word.end(), std::ostream_iterator<T>(strm, "\t"));

    return strm;
  }
};

// Reduce a word by cancelling adjacent inverses.
template<class T>
inline free_word<T>& free_word<T>::reduce()
{
  // This is the "equality" comparison class that allows cancellation
  // of generators.
  struct IsInv {
    bool operator() (const T& a, const T& b) { return (a == inverse(b)); }
  } isinv;

  // Keep cancelling until nothing changes.
  while (true)
    {
      auto i = adjacent_remove(word.begin(),word.end(),isinv);
      if (i == word.end()) break;
      word.erase(i,word.end());
    }

  return *this;
}

template<class T>
inline free_word<T> operator*(const free_word<T>& ww, const T& ee)
{
  free_word<T> ww2(ww.word);
  ww2.word.push_back(ee);
  return ww2;
}

template<class T>
inline free_word<T> operator*(const T& ee, const free_word<T>& ww)
{
  free_word<T> ww2(ww.word);
  ww2.word.push_front(ee);
  return ww2;
}




//template<class T> bool operator==(const T&, const T&);
//template<class T> T inverse(const T&);


class free_elem
{
protected:
  int inv;

public:
  free_elem(int inv_ = 1) : inv(inv_)
  {
    if (!(inv_ == 1 || inv_ == -1))
      {
	std::cerr << "Error: bad value " << inv_ <<" for inv.\n";
	std::exit(-1);
      }
  }

  friend bool operator==(const free_elem&, const free_elem&);
  friend free_elem inverse<>(const free_elem&);

  virtual std::ostream& print(std::ostream& strm) const
  {
    strm << "(" << "*" << (inv > 0 ? "+" : "-" ) << ")";
    return strm;
  }
};

class main_edge : public free_elem
{
  // Main edge is numbered by edge#.
  int eno;

public:
  main_edge(int eno_ = 0, int inv_ = 1) : free_elem(inv_), eno(eno_) {}

  friend bool operator==(const main_edge&, const main_edge&);
  friend main_edge inverse<>(const main_edge&);

  virtual std::ostream& print(std::ostream& strm) const
  {
    strm << "(" << eno << (inv > 0 ? "+" : "-" ) << ")";
    return strm;
  }
};

class inf_edge : public free_elem
{
  // Infinitesimal edge is numbered by multigon# and prong#.
  int mno, pno;

public:
  inf_edge(int mno_ = 0, int pno_ = 0, int inv_ = 1)
    : free_elem(inv_), mno(mno_), pno(pno_) {}

  friend bool operator==(const inf_edge&, const inf_edge&);
  friend inf_edge inverse<>(const inf_edge&);

  virtual std::ostream& print(std::ostream& strm) const
  {
    strm << "(" << mno << "," << pno << (inv > 0 ? "+" : "-" ) << ")";
    return strm;
  }
};


template<>
inline main_edge inverse(const main_edge& a)
{
  return main_edge(a.eno,-a.inv);
}

template<>
inline inf_edge inverse(const inf_edge& a)
{
  return inf_edge(a.mno,a.pno,-a.inv);
}

// A main edge is never equal to an infinitesimal one.
bool operator==(const main_edge&, const inf_edge&) { return false; }
bool operator==(const inf_edge&, const main_edge&) { return false; }

// Main edge comparison function.
bool operator==(const main_edge& a, const main_edge& b)
{
  return (a.eno == b.eno && a.inv == b.inv);
}

// Infinitesimal edge comparison function.
bool operator==(const inf_edge& a, const inf_edge& b)
{
  return (a.mno == b.mno && a.pno == b.pno && a.inv == b.inv);
}

#if 0
std::ostream& operator<<(std::ostream& strm, const free_elem& a)
{
  strm << "(" << "*" << (a.inv > 0 ? "+" : "-" ) << ")";
  return strm;
}

std::ostream& operator<<(std::ostream& strm, const main_edge& a)
{
  strm << "(" << a.eno << (a.inv > 0 ? "+" : "-" ) << ")";
  return strm;
}

std::ostream& operator<<(std::ostream& strm, const inf_edge& a)
{
  strm << "(" << a.mno << "," << a.pno << (a.inv > 0 ? "+" : "-" ) << ")";
  return strm;
}
#else
template<class T>
std::ostream& operator<<(std::ostream& strm, const free_word<T>& w)
{
  return w.print(strm);
}

std::ostream& operator<<(std::ostream& strm, const free_elem& a)
{
  return a.print(strm);
}

std::ostream& operator<<(std::ostream& strm, const main_edge& a)
{
  return a.print(strm);
}

std::ostream& operator<<(std::ostream& strm, const inf_edge& a)
{
  return a.print(strm);
}
#endif


} // namespace ttauto

#endif // FREEWORD_HPP
