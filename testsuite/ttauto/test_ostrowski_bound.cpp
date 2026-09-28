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

// Checks ttauto::ostrowski_schneider_bound, the lower bound behind the
// Ostrowski-Schneider prune test (issue #23), against its definition:
// the smallest weighted mean sum_i s_i x_i / sum_i x_i of the sums s over
// the box delta <= x_i <= 1.  The weighted mean is a ratio of linear
// functions, so its minimum over the box is at a vertex, and brute force
// over the 2^m vertices gives the exact value.  Also checks the
// properties the prune test relies on: the bound does not change when
// the sums are permuted, does not decrease when delta or a sum grows, and
// lies between the smallest sum and the mean.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
#include "check.hpp"
#include "traintracks/traintrack.hpp"
#include "ttauto/ttauto.hpp"

typedef ttauto::ttauto<traintracks::traintrack> tt;

static double bound(std::vector<int> s, const double delta)
{
  return tt::ostrowski_schneider_bound(s,delta);
}

static double brute(const std::vector<int>& s, const double delta)
{
  const int m = (int)s.size();
  double best = 0;
  for (int mask = 0; mask < (1 << m); ++mask)
    {
      double num = 0, den = 0;
      for (int i = 0; i < m; ++i)
        {
          const double x = (mask & (1 << i)) ? 1 : delta;
          num += s[i]*x; den += x;
        }
      if (mask == 0 || num/den < best) best = num/den;
    }
  return best;
}

int main()
{
  const double eps = 1e-12;
  std::mt19937 gen(20260928);
  std::uniform_int_distribution<int> size(1,8), sum(1,200);
  std::uniform_real_distribution<double> unit(0,1);
  int cases = 0;
  for (int trial = 0; trial < 20000; ++trial)
    {
      const int m = size(gen);
      std::vector<int> s(m);
      for (int& v : s) v = sum(gen);
      // Weight ratios near 0 (as in the search, Lambda^-(d-1)), and up to 1.
      const double delta = trial % 2 ? unit(gen) : std::pow(unit(gen),8);

      const double b = bound(s,delta);
      CHECK_MSG(std::abs(b - brute(s,delta)) <= eps*brute(s,delta),
                "m=" << m << " delta=" << delta);

      const int smin = *std::min_element(s.begin(),s.end());
      double mean = 0;
      for (int v : s) mean += v;
      mean /= m;
      CHECK(b >= smin*(1 - eps) && b <= mean*(1 + eps));

      std::vector<int> t = s;
      std::shuffle(t.begin(),t.end(),gen);
      CHECK(std::abs(bound(t,delta) - b) <= eps*b);

      CHECK(bound(s,std::min(1.0,delta*1.5 + 1e-3)) >= b*(1 - eps));
      t = s;
      ++t[trial % m];
      CHECK(bound(t,delta) >= b*(1 - eps));
      ++cases;
    }

  // delta = 1 gives the mean; delta = 0 gives the smallest sum.
  CHECK(std::abs(bound({1,2,3,6},1) - 3) < eps);
  CHECK(std::abs(bound({4,2,9},0) - 2) < eps);

  // The sums are sorted in place, and the caller's buffer reused.
  std::vector<int> buf = {5,1,4};
  tt::ostrowski_schneider_bound(buf,0.5);
  CHECK((buf == std::vector<int>{1,4,5}));

  std::cout << "test_ostrowski_bound: " << cases
            << " random cases agree with brute force" << std::endl;
  return 0;
}
