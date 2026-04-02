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

#ifndef TRAINTRACKS_UTIL_HPP
#define TRAINTRACKS_UTIL_HPP

#include <jlt/vector.hpp>

namespace traintracks {


// If two vectors are cyclically equivalent, return a vector p0v of
// offsets between them such that v1[v] == v2[(v+p0v[i]) % size()].
// Return empty p0v if they are not cyclically equivalent.
template<class Vec>
jlt::vector<int> cyclic_shift(const Vec& v1, const Vec& v2)
{
  if (v1.size() != v2.size()) return jlt::vector<int>();

  jlt::vector<int> p0v;

  const int k = v1.size();

  // Shift the starting point of the second path to test for 'cyclic
  // equality'.
  // Try all possibilities, and return a list of matching offsets.
  for (int p0 = 0; p0 < k; ++p0)
    {
      int p = 0;
      for (; p < k; ++p)
	{
	  int p1 = (p+p0) % k;
	  if (!(v1[p] == v2[p1])) break;
	}
      // We made it through without breaking out, so must tbe equal.
      if (p == k) p0v.push_back(p0);
    }

  return p0v;
}


template<class Vec>
inline bool cyclic_equal(const Vec& v1, const Vec& v2)
{
  return !(cyclic_shift(v1,v2).empty());
}


// Small function to mod and make the result nonnegative.
template<class Int, class Int2>
inline Int2 mod(const Int m, const Int2 n)
{
  // Be careful with unsigned types for Int!
  Int mm = m % n;
  // mm has a value between -(n-1) and (n-1).  Make it nonnegative.
  return (Int2)(mm >= 0 ? mm : mm + n);
}

} // namespace traintracks

#endif // TRAINTRACKS_UTIL_HPP
