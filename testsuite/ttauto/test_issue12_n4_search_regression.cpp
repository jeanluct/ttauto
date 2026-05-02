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

#include <cmath>
#include <iostream>
#include <list>
#include <algorithm>
#include "traintracks/build.hpp"
#include "traintracks/traintrack.hpp"
#include "ttauto/ttauto.hpp"
#include "ttauto/ttfoldgraph.hpp"

#define REQUIRE(cond) \
  do { \
    if (!(cond)) { \
      std::cerr << "Requirement failed: " << #cond << " at " \
                << __FILE__ << ":" << __LINE__ << "\n"; \
      return 1; \
    } \
  } while (0)

int main()
{
  using traintracks::traintrack;
  using ttauto::ttfoldgraph;

  auto ttv = traintracks::build_traintrack_list(4);
  REQUIRE(static_cast<int>(ttv.size()) >= 1);

  ttfoldgraph<traintrack> full(ttv[0]);
  REQUIRE(full.vertices() == 4);

  std::list<ttfoldgraph<traintrack> > sub = ttauto::subgraphs(full);
  REQUIRE(sub.size() == 2);

  auto it = std::max_element(sub.begin(),sub.end(),
                             [](const ttfoldgraph<traintrack>& a,
                                const ttfoldgraph<traintrack>& b)
                             { return a.vertices() < b.vertices(); });
  REQUIRE(it != sub.end());
  REQUIRE(it->vertices() == 3);

  ttauto::ttauto<traintrack> search(*it);
  search.max_pathlength(5);
  search.badword_length(0);
  search.search();

  bool found_2618 = false;
  for (auto pit = search.pA_list().begin(); pit != search.pA_list().end(); ++pit)
    {
      if (std::abs(pit->second.dilatation() - 2.61803398875) < 1e-4)
        {
          found_2618 = true;
          break;
        }
    }

  REQUIRE(found_2618);
  return 0;
}
