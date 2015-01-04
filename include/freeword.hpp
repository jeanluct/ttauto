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
template<class T>
free_word<T> operator*(const free_word<T>&, const T&);
template<class T>
free_word<T> operator*(const T&, const free_word<T>&);
template<class T>
free_word<T> operator*(const free_word<T>&, const free_word<T>&);


template<class T>
class free_word : public std::list<T>
{
  int ngen;

public:
  free_word(int ngen_ = 1) : ngen(ngen_) {}

  // Forward initializer list to base class.
  free_word(std::initializer_list<T> l, int ngen_ = 0) : std::list<T>(l)
  {
    /* This is specialized to ints for now. */
    int ngen = 0;
    for (auto i : *this)
      {
	ngen = std::max(std::abs(i),ngen);
      }
    if (ngen_ && ngen > ngen_)
      {
	std::cerr << "ngen is too small.\n";
	std::exit(-1);
      }
    ngen = std::max(ngen,ngen_);
  }

  /* Make a copy constructor and operator= that check gen. */

  // Reduce a word by cancelling adjacent inverses.
  free_word<T>& reduce();

  // Return the inverse of the free word.
  free_word<T> inverse() const
  {
    free_word<T> iw(ngen);
    for (auto i = this->rbegin(); i != this->rend(); ++i) iw.push_back(-(*i));
    return iw;
  }

  // Multiplication just appends and prepends elements to the word.
  friend free_word<T> operator*<>(const free_word<T>&, const T&);
  friend free_word<T> operator*<>(const T&, const free_word<T>&);
  friend free_word<T> operator*<>(const free_word<T>&, const free_word<T>&);
};

// Reduce a word by cancelling adjacent inverses.
template<class T>
inline free_word<T>& free_word<T>::reduce()
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

template<class T>
inline free_word<T> operator*(const free_word<T>& w1, const free_word<T>& w2)
{
  free_word<T> w12 = w1;
  w12.insert(w12.end(),w2.begin(),w2.end());
  return w12;
}


template<class T>
class free_auto : std::vector<free_word<T> >
{
public:
  free_auto(const int ngen_) : std::vector<free_word<T> >(ngen_+1)
  {
    // Initialize to the identity automorphism.
    for (size_t i = 1; i <= numgens(); ++i)
      std::vector<free_word<T> >::operator[](i).push_back(i);
  }

  // Construct from base class.
  free_auto(const std::vector<free_word<T> >& a) : std::vector<free_word<T> >(a)
  {
  }

  // Number of generators.  Note that the vector size has to be one larger.
  size_t numgens() const { return (this->size()-1); }

  // Action of the automorphism on a single generator a.
  // Use this for assignment.
  free_word<T>& operator[](const T a)
  {
    if (a == 0)
      {
	std::cerr << "Cannot specify action on identity.\n";
	std::exit(-1);
      }
    if (a < 0)
      {
	std::cerr << "Assign to positive generators only.\n";
	std::exit(-1);
      }
    // Only store the action on positive generators.
    return std::vector<free_word<T> >::operator[](a);
  }

  // Const version of operator[].
  // Must still return a reference, so that pointer operations like
  // begin() are still valid.
  const free_word<T>& operator[](const T a) const
  {
    if (a == 0)
      {
	std::cerr << "Cannot specify action on identity.\n";
	std::exit(-1);
      }
    if (a < 0)
      {
	std::cerr << "Assign to positive generators only.\n";
	std::exit(-1);
      }
    // Only store the action on positive generators.
    return std::vector<free_word<T> >::operator[](a);
  }

  // get_action is like operator[], but returns by value.
  // It can thus return a free_word for a<0.
  const free_word<T> get_action(const T a) const
  {
    if (a == 0)
      {
	std::cerr << "Cannot specify action on identity.\n";
	std::exit(-1);
      }
    if (a > 0)
      return std::vector<free_word<T> >::operator[](a);
    else // For negative a, call free_word<T>::inverse().
      return (std::vector<free_word<T> >::operator[](-a)).inverse();
  }

  free_auto<T>& operator*=(const free_auto<T>& a)
  {
    if (numgens() != a.numgens())
      {
	std::cerr << "Can only compose free_auto objects with ";
	std::cerr << "the same number of generators.\n";
	std::exit(-1);
      }
    free_auto<T> res(std::vector<free_word<T> >(numgens()+1));
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
inline free_auto<T> operator*(const free_auto<T>& a, const free_auto<T>& b)
{
  return (free_auto<T>(a) *= b);
}

template<class T>
std::ostream& operator<<(std::ostream& strm, const free_auto<T>& a)
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

#endif // FREEWORD_HPP
