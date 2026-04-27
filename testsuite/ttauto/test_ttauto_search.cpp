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

#include <algorithm>
#include <cassert>
#include <sstream>
#include "traintracks/build.hpp"
#include "traintracks/traintrack.hpp"
#include "ttauto/pAclass.hpp"
#include "ttauto/ttauto.hpp"
#include "ttauto/ttfoldgraph.hpp"


int main()
{
  using traintracks::build_traintrack_list;
  using traintracks::traintrack;
  using ttauto::ttauto;
  using ttauto::ttfoldgraph;

  // Deterministic small fixture: same strata family as examples/ttauto
  // but with a tiny graph and bounded path length for fast CTest runs.
  const int n = 5;
  const int trk = 0;
  auto ttv = build_traintrack_list(n);
  assert(!ttv.empty());
  assert(static_cast<int>(ttv.size()) > trk);

  ttfoldgraph<traintrack> ttg(ttv[trk]);
  assert(ttg.vertices() > 0);

  ttauto<traintrack> tta(ttg);
  tta.max_pathlength(8);
  tta.max_dilatation(3.0);
  tta.badword_length(0);
  tta.search();

  // We avoid checking exact totals because tiny internal changes can produce
  // slightly different counts while still being correct. Instead, we verify
  // that every result has the basic properties we expect.
  const auto& pal = tta.pA_list();
  assert(!pal.empty());
  assert(pal.size() < 1000);

  for (auto it = pal.begin(); it != pal.end(); ++it)
    {
      [[maybe_unused]] const auto& cls = it->second;
      assert(cls.number_of_paths() > 0);
      assert(cls.shortest() > 0);
      assert(cls.longest() >= cls.shortest());
      assert(cls.dilatation() > 1.0);
    }

  // Ensure reporting helpers stay callable and non-empty on valid results.
  std::ostringstream out;
  tta.print_pA_list(out);
  assert(!out.str().empty());

  out.str("");
  out.clear();
  tta.print_pA_list_MathematicaForm(out);
  assert(!out.str().empty());

  return 0;
}
