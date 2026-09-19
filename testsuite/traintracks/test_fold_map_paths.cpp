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

// One-fold train-track maps are continuous edge paths in the canonical
// numbering, abelianise to the transition matrix, and compose along paths
// (including through automorphic vertices) into continuous paths that map
// tails to tails and heads to heads.

#include <algorithm>
#include <iostream>
#include <random>
#include <vector>
#include <jlt/freeauto.hpp>
#include "check.hpp"
#include "traintracks/build.hpp"
#include "traintracks/map.hpp"
#include "traintracks/traintrack.hpp"
#include "ttauto/folding_path.hpp"
#include "ttauto/ttfoldgraph.hpp"

using jlt::freeauto;
using jlt::freeword;
using traintracks::traintrack;
using traintracks::ttnumbering;
using ttauto::folding_path;
using ttauto::ttfoldgraph;

// Every image word alternates main and side letters, begins and ends with
// a main letter, and consecutive letters share an endpoint in numbering N.
static void check_words_continuous(const freeauto<int>& AM, const ttnumbering& N,
                                   const int nmain, const char* what)
{
  for (int g = 1; g <= nmain; ++g)
    {
      const freeword<int> w = AM.get_action(g);
      std::vector<int> v(w.begin(),w.end());
      CHECK_MSG(!v.empty(), what);
      CHECK_MSG(N.is_main(v.front()) && N.is_main(v.back()), what);
      for (std::size_t i = 0; i < v.size(); ++i)
        {
          const bool main = N.is_main(v[i]);
          CHECK_MSG(main || N.is_side(v[i]), what);
          CHECK_MSG(main == (i % 2 == 0), what);
          if (i + 1 < v.size())
            CHECK_MSG(N.head_of(v[i]) == N.tail_of(v[i+1]), what);
        }
    }
  // Sides map to single sides.
  for (int q = 0; q < N.nprongs(); ++q)
    {
      const freeword<int> w = AM.get_action(N.side_letter(q));
      CHECK_MSG(w.size() == 1 && N.is_side(*w.begin()) && *w.begin() > 0, what);
    }
}

// Image of edge e starts at the image of its tail prong and ends at the
// image of its head prong (prong images read off the side letters).
static void check_endpoints_map(const freeauto<int>& AM, const ttnumbering& Nfrom,
                                const ttnumbering& Nto, const char* what)
{
  const int nmain = Nfrom.nedges();
  auto prong_image = [&](int q) {
    return Nto.side_of(*AM.get_action(Nfrom.side_letter(q)).begin());
  };
  for (int e = 0; e < nmain; ++e)
    {
      const freeword<int> w = AM.get_action(e+1);
      CHECK_MSG(Nto.tail_of(*w.begin()) == prong_image(Nfrom.edge_tail[e]), what);
      CHECK_MSG(Nto.head_of(*w.rbegin()) == prong_image(Nfrom.edge_head[e]), what);
    }
  // Prong images form a bijection preserving multigon type and side order.
  std::vector<int> img(Nfrom.nprongs());
  for (int q = 0; q < Nfrom.nprongs(); ++q) img[q] = prong_image(q);
  std::vector<int> sorted(img); std::sort(sorted.begin(),sorted.end());
  for (int q = 0; q < Nfrom.nprongs(); ++q) CHECK_MSG(sorted[q] == q, what);
  for (int q = 0; q < Nfrom.nprongs(); ++q)
    {
      CHECK_MSG(Nto.prong[img[q]].nprongs == Nfrom.prong[q].nprongs, what);
      CHECK_MSG(Nto.prong[img[q]].punctured == Nfrom.prong[q].punctured, what);
      CHECK_MSG(Nto.side_to(img[q]) == img[Nfrom.side_to(q)], what);
    }
}

static void check_graph(const int n, const int trk, std::mt19937& rng,
                        const int npaths, const int maxlen, int& nbranches)
{
  jlt::vector<traintrack> ttv = traintracks::build_traintrack_list(n);
  ttfoldgraph<traintrack> ttg(ttv[trk]);
  const int nmain = ttg.edges();

  // Every branch.
  for (int v = 0; v < ttg.vertices(); ++v)
    {
      const ttnumbering Nv = ttg.traintrack(v).numbering();
      for (int b = 0; b < ttg.foldings(v); ++b)
        {
          const freeauto<int> AM = ttg.traintrack_map(v,b);
          const ttnumbering Nt = ttg.traintrack(ttg.target_vertex(v,b)).numbering();
          check_words_continuous(AM,Nt,nmain,"one-step word");
          check_endpoints_map(AM,Nv,Nt,"one-step endpoints");
          CHECK_MSG(traintracks::transition_matrix_from_map(ttg.traintrack(v),AM)
                    == ttg.transition_matrix(v,b).full(), "one-step matrix");
          // Exactly one side letter, in exactly one image, for a real fold.
          int nsideletters = 0;
          for (int g = 1; g <= nmain; ++g)
            for (auto x : AM.get_action(g)) if (Nt.is_side(x)) ++nsideletters;
          CHECK_MSG(nsideletters == 1, "one side letter per fold");
          ++nbranches;
        }
    }

  // Random paths: composed words are continuous, map endpoints correctly
  // and abelianise to the path matrix; closed paths iterate continuously.
  std::uniform_int_distribution<int> pick_v(0,ttg.vertices()-1);
  std::uniform_int_distribution<int> pick_len(1,maxlen);
  for (int i = 0; i < npaths; ++i)
    {
      const int v0 = pick_v(rng);
      folding_path<traintrack> p(ttg,v0);
      const int len = pick_len(rng);
      for (int j = 0; j < len; ++j)
        {
          std::uniform_int_distribution<int> pick_b(0,ttg.foldings(p.final_vertex())-1);
          p.push_back(pick_b(rng));
        }
      const freeauto<int> AM = p.traintrack_map();
      const ttnumbering N0 = ttg.traintrack(v0).numbering();
      const ttnumbering N1 = ttg.traintrack(p.final_vertex()).numbering();
      check_words_continuous(AM,N1,nmain,"path word");
      check_endpoints_map(AM,N0,N1,"path endpoints");
      CHECK_MSG(traintracks::transition_matrix_from_map(ttg.traintrack(v0),AM)
                == p.transition_matrix(), "path matrix");
      if (p.closed())
        {
          freeauto<int> P(AM);
          for (int k = 0; k < 2; ++k)
            {
              P = P * AM;
              check_words_continuous(P,N0,nmain,"iterated closed path");
            }
        }
    }
}

int main()
{
  std::mt19937 rng(20260919);
  int nbranches = 0;
  check_graph(3,0,rng,40,6,nbranches);
  check_graph(4,1,rng,40,8,nbranches);
  // n=5: several strata, including tracks with cyclic symmetry.
  for (int trk = 0; trk < 4; ++trk) check_graph(5,trk,rng,25,8,nbranches);
  CHECK(nbranches > 0);
  std::cout << "test_fold_map_paths: OK (" << nbranches << " branches checked)\n";
  return 0;
}
