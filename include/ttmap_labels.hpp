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

#ifndef TTMAP_LABELS_HPP
#define TTMAP_LABELS_HPP

#include <cstdlib>
#include <iostream>

namespace traintracks {

// Unified generator labeling for train-track maps.
//
// Main edges are numbered 1..nmain.
// Infinitesimal edges are numbered nmain+1..nmain+ninf.
// (Use "peripheral" specifically for edges around monogons.)
// Orientation is encoded by sign (negative means inverse orientation).

struct ttmap_labeler
{
  // Number of main edges.
  const int nmain;
  // Number of infinitesimal edges.
  const int ninf;

  ttmap_labeler(const int nmain_, const int ninf_)
    : nmain(nmain_), ninf(ninf_)
  {
    if (nmain_ < 0 || ninf_ < 0)
      {
        std::cerr << "Negative generator counts in traintracks::ttmap_labeler::ttmap_labeler.\n";
        std::exit(1);
      }
  }

  // Total generator count (main + infinitesimal).
  int num_generators() const { return nmain + ninf; }

  // True when g is a nonzero generator index in range.
  bool is_valid_generator(const int g) const
  {
    return (g != 0 && std::abs(g) <= num_generators());
  }

  // True when g indexes a main-edge generator (up to sign).
  bool is_main_generator(const int g) const
  {
    return (g != 0 && std::abs(g) <= nmain);
  }

  // Convert signed main generator label to zero-based edge index.
  int main_generator_index(const int g) const
  {
    if (!is_main_generator(g))
      {
        std::cerr << "Not a main generator in traintracks::ttmap_labeler::main_generator_index.\n";
        std::exit(1);
      }
    return std::abs(g) - 1;
  }

  // Convert zero-based main-edge index to positive generator label.
  int main_gen(const int edge_index) const
  {
    if (edge_index < 0 || edge_index >= nmain)
      {
        std::cerr << "Main edge index out of range in traintracks::ttmap_labeler::main_gen.\n";
        std::exit(1);
      }
    return edge_index + 1;
  }

  // Convert zero-based infinitesimal index to positive generator label.
  int infinitesimal_gen(const int infinitesimal_index) const
  {
    if (infinitesimal_index < 0 || infinitesimal_index >= ninf)
      {
        std::cerr << "Infinitesimal index out of range in traintracks::ttmap_labeler::infinitesimal_gen.\n";
        std::exit(1);
      }
    return nmain + infinitesimal_index + 1;
  }
};

} // namespace traintracks

#endif // TTMAP_LABELS_HPP
