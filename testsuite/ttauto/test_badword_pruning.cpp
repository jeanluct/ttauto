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

// Bad-word pruning: what it buys and what it costs.
//
// A bad word is the square h*h of a closed subpath h whose transition
// matrix has the same pattern of zeros as its square (badwords.hpp).
// Traversing h twice instead of once can only raise the dilatation, so
// the minimiser never contains one and the prune is safe when hunting a
// minimum -- but the repeated path is usually a pseudo-Anosov class in
// its own right, so the prune is NOT safe when enumerating every class.
// See ttauto::badword_length and doc/ttauto.tex.
//
// Both halves are checked here, because the second is the reason the
// prune is off by default, and prose alone has not kept people from
// assuming it is free.
//
// Before 2026-09 this block could not be tested at all: pruning was
// silently gated on check_norms as well, and every caller disabled it.

#include <cmath>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <jlt/vector.hpp>
#include "traintracks/build.hpp"
#include "traintracks/traintrack.hpp"
#include "ttauto/ttfoldgraph.hpp"
#include "ttauto/ttauto.hpp"
#include "check.hpp"

using traintracks::traintrack;
using ttgraph = ttauto::ttfoldgraph<traintrack>;

// One length-bounded search, with bad words on or off.
struct outcome
{
  int classes = 0;
  double lowest = -1;
  long long omitted = 0;
  // Keyed by characteristic polynomial, which is what pAclass itself is
  // keyed on.  Comparing dilatations instead would be wrong: the same
  // class reached by a different path differs in the last few digits.
  std::set<std::string> keys;
};

static outcome search_with(const ttgraph& ttg, const int bwl, const int len)
{
  ttauto::ttauto<traintrack> tta(ttg);
  tta.badword_length(bwl).max_pathlength(len);

  // The search and the bad-word table both report progress on cout.
  std::ostringstream sink;
  std::streambuf* const saved = std::cout.rdbuf(sink.rdbuf());
  tta.search();
  std::cout.rdbuf(saved);

  outcome o;
  o.classes = (int)tta.pA_list().size();
  o.omitted = tta.badwords_omitted();
  if (!tta.pA_list().empty())
    o.lowest = tta.pA_list().begin()->second.dilatation();
  for (auto it = tta.pA_list().begin(); it != tta.pA_list().end(); ++it)
    {
      std::ostringstream poly;
      poly << it->first;
      o.keys.insert(poly.str());
    }
  return o;
}

int main()
{
  using std::cout;
  using std::endl;

  // Four punctures, second stratum: three vertices, so a path of length
  // 10 has room to go round a loop several times.  Small enough that both
  // searches together take a few hundredths of a second.
  jlt::vector<traintrack> ttv = traintracks::build_traintrack_list(4);
  CHECK((int)ttv.size() > 1);
  const ttgraph ttg(ttv[1]);
  CHECK(ttg.vertices() == 3);

  const int len = 10;
  const outcome off = search_with(ttg,0,len);
  const outcome on  = search_with(ttg,1,len);

  cout << "n=4 stratum 2, " << ttg.vertices() << " vertices, length "
       << len << ":\n";
  cout << "  without bad words: " << off.classes << " classes, min "
       << off.lowest << "\n";
  cout << "  with bad words:    " << on.classes << " classes, min "
       << on.lowest << ", " << on.omitted << " branches cut\n";

  // The prune has to actually fire, or nothing below means anything.
  CHECK(off.omitted == 0);
  CHECK(on.omitted > 0);

  // What it buys: the minimum is untouched.  This is the claim the paper
  // makes (a repeated loop can only raise the dilatation) and it had
  // never been executed.
  CHECK(off.classes > 0);
  CHECK(std::fabs(on.lowest - off.lowest) < 1e-10);
  CHECK(std::fabs(off.lowest - 2.29663) < 1e-5);

  // What it costs: classes go missing, and not marginally.  Pruning can
  // only remove branches, so what survives must be a subset of what an
  // unpruned search finds; anything else would mean the prune is
  // reaching paths it should not.
  CHECK(on.classes < off.classes);
  CHECK((int)on.keys.size() == on.classes);
  CHECK((int)off.keys.size() == off.classes);
  for (auto it = on.keys.begin(); it != on.keys.end(); ++it)
    CHECK_MSG(off.keys.count(*it) == 1, *it);

  // Pin the numbers, so a change in either direction is noticed.
  CHECK(off.classes == 62);
  CHECK(on.classes == 39);

  // Off is the default: a freshly built search prunes nothing.
  {
    ttauto::ttauto<traintrack> fresh(ttg);
    fresh.max_pathlength(len);
    std::ostringstream sink;
    std::streambuf* const saved = std::cout.rdbuf(sink.rdbuf());
    fresh.search();
    std::cout.rdbuf(saved);
    CHECK(fresh.badwords_omitted() == 0);
    CHECK((int)fresh.pA_list().size() == off.classes);
  }

  // check_norms() must not turn pruning on behind the caller's back, as
  // it used to by rebuilding the table from a default of 2.
  {
    ttauto::ttauto<traintrack> norms(ttg);
    norms.max_dilatation(3.0).check_norms().max_pathlength(len);
    std::ostringstream sink;
    std::streambuf* const saved = std::cout.rdbuf(sink.rdbuf());
    norms.search();
    std::cout.rdbuf(saved);
    CHECK(norms.badwords_omitted() == 0);
  }

  cout << "\ntest_badword_pruning: all checks passed" << endl;
  return 0;
}
