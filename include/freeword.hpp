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

template <class ForwardIterator>
ForwardIterator reduce(ForwardIterator first, ForwardIterator last)
{
  if (first == last) return last;

  ForwardIterator result = first;
  while (++first != last)
    {
      if (!(*result == *first))
	{
	  *(++result) = *first;
	}
      else
	{
	  *result = *(++first);
	  if (first == last) { --result; break; }
	}
    }
  return ++result;
}

class free_elem
{
protected:
  int inv;

public:
  free_elem(int inv_ = 0) : inv(inv_) {}

  free_elem inverse() { return free_elem(-inv); }
};

// Forward declarations.
template<class T> class free_word;
template<class T> free_word<T> operator*(const free_word<T>&, const T&);
template<class T> free_word<T> operator*(const T&, const free_word<T>&);

#if 0
template<class T>
class free_word
{
  std::list<T> w;

public:
  free_word(std::list<T> ww = std::list<T>()) : w(ww) {}

  free_word(std::initializer_list<T> l) : w(l) {}

  free_word<T> reduce()
  {
    for (auto j = w.begin(), i = j++; j != w.end; ++i, ++j)
      {
	if (*i == j->inverse())
	  {
	    std::remove(i,j);
	  }
      }
  return *this;
  }

  friend free_word<T> operator*<>(const free_word<T>&, const T&);
  friend free_word<T> operator*<>(const T&, const free_word<T>&);
};

template<class T>
inline free_word<T> operator*(const free_word<T>& ww, const T& ee)
{
  std::list<T> w = ww.w;
  w.push_back(ee);
  return free_word<T>(w);
}

template<class T>
inline free_word<T> operator*(const T& ee, const free_word<T>& ww)
{
  std::list<T> w = ww.w;
  w.push_front(ee);
  return free_word<T>(w);
}
#else
template<class T>
class free_word : public std::list<T>
{
public:
  // Forward initializer list.
  free_word(std::initializer_list<T> l) : std::list<T>(l) {}

  free_word<T> reduce()
  {
    for (auto j = this->begin(), i = j++; j != this->end(); ++i, ++j)
      {
	if (*i == j->inverse())
	  {
	    std::remove(i,j);
	  }
      }
  return *this;
  }

  friend free_word<T> operator*<>(const free_word<T>&, const T&);
  friend free_word<T> operator*<>(const T&, const free_word<T>&);
};

template<class T>
inline free_word<T> operator*(const free_word<T>& ww, const T& ee)
{
  free_word<T> ww2 = ww;
  ww2.push_back(ee);
  return ww2;
}

template<class T>
inline free_word<T> operator*(const T& ee, const free_word<T>& ww)
{
  free_word<T> ww2 = ww;
  ww2.push_front(ee);
  return ww2;
}
#endif


class main_edge : public free_elem
{
  // Main edge is numbered by edge#.
  int eno;

public:
  main_edge(int eno_ = 0, int inv_ = 0) : free_elem(inv_), eno(eno_) {}

  /*
  main_edge inv()
  {
    return main_edge(eno,-inv);
  }
  */
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
