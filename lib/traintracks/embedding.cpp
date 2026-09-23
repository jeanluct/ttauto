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

// The boundary walk of a train track and the order it puts the punctures
// in.  See include/traintracks/embedding.hpp for what this is for.

#include <cstdlib>
#include <iostream>
#include <vector>

#include "traintracks/embedding.hpp"

namespace traintracks {

std::vector<int> all_directions_at_prong(const ttnumbering& num, const int q)
{
  std::vector<int> d;
  d.push_back(-num.side_letter(num.side_from_prev(q)));
  d.insert(d.end(),num.prong_letters[q].begin(),num.prong_letters[q].end());
  d.push_back(num.side_letter(q));
  return d;
}

namespace {

// Where each dart sits in the cyclic order at its tail prong, indexed by
// letter through dart_index().
struct dart_position
{
  int prong;
  int slot;
};

inline int dart_index(const int letter, const int nletters)
{
  return (letter > 0 ? letter-1 : nletters - letter - 1);
}

} // namespace

tt_embedding outer_embedding(const ttnumbering& num, const int cut_dart)
{
  const int nprongs = num.nprongs(), nedges = num.nedges();
  const int nletters = nedges + nprongs;

  // The rotation system: the cyclic order of the darts leaving each prong.
  std::vector<std::vector<int> > dirs(nprongs);
  std::vector<dart_position> where(2*nletters,dart_position{-1,-1});
  for (int q = 0; q < nprongs; ++q)
    {
      dirs[q] = all_directions_at_prong(num,q);
      for (int i = 0; i < (int)dirs[q].size(); ++i)
        where[dart_index(dirs[q][i],nletters)] = dart_position{q,i};
    }

  // Where to cut.  Dart 0 means the canonical choice, the loop of the root
  // monogon of the coding, which is prong 0.
  const int d0 = (cut_dart == 0 ? num.side_letter(0) : cut_dart);
  if (!num.is_main(d0) && !(num.is_side(d0) && d0 > 0))
    {
      std::cerr << "Dart " << d0 << " is not on the boundary region"
                << " in traintracks::outer_embedding.\n";
      std::exit(1);
    }

  tt_embedding emb;
  emb.cut_dart = d0;
  emb.position_of.assign(nprongs,0);

  // Walk the boundary region.  The face following the dart d continues
  // along the dart after -d in the cyclic order at the head of d; with the
  // coding's cyclic order that convention traverses a complementary polygon
  // along its negative side letters, and the boundary region along the
  // positive ones and along both darts of every main edge.
  int d = d0;
  do
    {
      emb.walk.push_back(d);

      if (num.is_side(d) && d > 0)
        {
          const int q = num.side_of(d);
          if (num.prong[q].punctured)
            {
              if (num.prong[q].nprongs != 1)
                {
                  std::cerr << "Punctured multigon with " << num.prong[q].nprongs
                            << " prongs in traintracks::outer_embedding:"
                            << " only monogon punctures are supported.\n";
                  std::exit(1);
                }
              emb.position_of[q] = (int)emb.puncture_order.size() + 1;
              emb.puncture_order.push_back(q);
            }
        }

      const dart_position& w = where[dart_index(-d,nletters)];
      const std::vector<int>& list = dirs[w.prong];
      d = list[(w.slot+1) % (int)list.size()];
    }
  while (d != d0);

  // The boundary region runs twice along every main edge, since the
  // multigons and main edges form a tree and so every main edge is a
  // bridge, and once along every side.  Any other length means the
  // rotation system is not the one this assumes.
  if ((int)emb.walk.size() != 2*nedges + nprongs)
    {
      std::cerr << "Boundary walk has length " << emb.walk.size()
                << ", expected " << 2*nedges + nprongs
                << ", in traintracks::outer_embedding.\n";
      std::exit(1);
    }

  // A corner of the boundary region is an exterior cusp when both its
  // letters are main: the two smooth corners at each prong are the ones
  // where a main letter meets a side.
  const int L = (int)emb.walk.size();
  for (int i = 0; i < L; ++i)
    {
      if (num.is_main(emb.walk[i]) && num.is_main(emb.walk[(i+1) % L]))
        emb.cusp_after.push_back(i);
    }

  return emb;
}

int transported_cut_dart(const ttnumbering& before, const fold_map_data& fm,
                         const int cut_dart)
{
  if (before.is_main(cut_dart))
    {
      const int e = before.edge_of(cut_dart);
      const int im = fm.edge_image[e];
      return (cut_dart > 0 ? im : -im);
    }
  if (before.is_side(cut_dart) && cut_dart > 0)
    return fm.nmain + 1 + fm.prong_image[before.side_of(cut_dart)];

  std::cerr << "Dart " << cut_dart << " is not on the boundary region"
            << " in traintracks::transported_cut_dart.\n";
  std::exit(1);
}

std::ostream& tt_embedding::print(std::ostream& strm) const
{
  strm << "cut at dart " << cut_dart
       << "\npunctures (prong numbers, position 1 first):";
  for (int i = 0; i < (int)puncture_order.size(); ++i)
    strm << " " << puncture_order[i];
  strm << "\nboundary walk:";
  for (int i = 0; i < (int)walk.size(); ++i) strm << " " << walk[i];
  strm << "\nexterior cusps after walk positions:";
  for (int i = 0; i < (int)cusp_after.size(); ++i) strm << " " << cusp_after[i];
  return strm << "\n";
}

} // namespace traintracks
