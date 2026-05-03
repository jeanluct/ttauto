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

#include <cassert>
#include <jlt/vector.hpp>
#include "traintracks/build.hpp"
#include "traintracks/traintrack.hpp"


int main()
{
  using traintracks::traintrack;

  // Build a nontrivial train track and check initial consistency.
  const int N = 9;
  jlt::vector<int> Kv(3);
  Kv[0] = 4;
  Kv[1] = 5;
  Kv[2] = 3;

  traintrack tt(N,Kv);
  if (traintrack::label_multiprongs)
    tt.set_label(0,8);
  tt.check();

  traintrack tt_copy(tt);
  assert(tt == tt_copy);

  // Repeated normalisation should preserve consistency.
  for (int i = 0; i < 10; ++i)
    {
      tt_copy.normalise();
      tt_copy.check();
    }

  traintrack tt_from_code(tt.coding());
  assert(tt == tt_from_code);

  // Mutation stress: repeatedly apply one valid fold and re-check invariants.
  traintrack tt_mut(tt_from_code);
  for (int step = 0; step < 20; ++step)
    {
      bool folded = false;
      for (int f = 0; f < tt_mut.foldings(); ++f)
        {
          traintrack candidate(tt_mut);
          if (!candidate.fold(f)) continue;
          candidate.check();
          candidate.normalise();
          candidate.check();
          tt_mut = candidate;
          folded = true;
          break;
        }
      if (!folded) break;
    }

  jlt::vector<double> w(tt.edges());
  for (int i = 0; i < tt.edges(); ++i) w[i] = i + 1;
  tt.weights(w.begin());

  // Weights traversal should remain aligned with edge count.
  assert(tt.weights().size() == static_cast<jlt::vector<double>::size_type>(tt.edges()));

  // A single canonical fold (if available) should produce a valid track.
  if (tt.foldings() > 0)
    {
      [[maybe_unused]] const bool ok = tt.fold(0);
      assert(ok);
      tt.check();
    }

  jlt::vector<traintrack> ttv = traintracks::build_traintrack_list(6);

  // Enumerated representatives should be non-empty and valid.
  assert(!ttv.empty());
  ttv[0].check();

  return 0;
}
