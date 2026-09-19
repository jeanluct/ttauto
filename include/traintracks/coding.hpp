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

#ifndef TRAINTRACKS_CODING_HPP
#define TRAINTRACKS_CODING_HPP

#include <iosfwd>
#include <jlt/vector.hpp>

namespace traintracks {

class traintrack;
class mathmatrix_permplus1;
class multigon;

namespace detail {

// Fixed-size block encoding one directed edge-step in canonical coding.
struct coding_block
{
  using coding_vec = jlt::vector<int>;
  static const int length = 5;

  int prong;
  int nprongs;
  int label;
  int edge;
  int nedges;

  coding_block(int p = 0, int np = 1, int lb = 0, int e = 0, int ne = 1)
    : prong(p), nprongs(np), label(lb), edge(e), nedges(ne) {}

  coding_block(coding_vec::const_iterator& ci)
    : prong(*ci++), nprongs(*ci++), label(*ci++), edge(*ci++), nedges(*ci++) {}

  void append_to(coding_vec& v) const
  {
    v.push_back(prong); v.push_back(nprongs);
    v.push_back(label);
    v.push_back(edge); v.push_back(nedges);
  }

  bool operator==(const coding_block& b) const
  {
    return (prong == b.prong && nprongs == b.nprongs &&
            label == b.label &&
            edge == b.edge && nedges == b.nedges);
  }

  bool operator!=(const coding_block& b) const { return !operator==(b); }
};

class coding_engine
{
public:
  using coding_vec = jlt::vector<int>;

  // Return canonical coding in orientation dir without mutating tt.
  static coding_vec coding(const traintrack& tt, int dir);

  // Rotate tt so monogon 0 is the lexicographic coding minimiser.
  static coding_vec minimise_coding(traintrack& tt);

  // Compute branch permutation induced by cyclic coding symmetry.
  static mathmatrix_permplus1 cyclic_symmetry(traintrack& tt);

  // Write canonical coding blocks to stream.
  static std::ostream& print_coding(const traintrack& tt,
                                    std::ostream& strm,
                                    int dir);

private:
  // Build coding sequence rooted at a chosen uncusped monogon.
  static coding_vec coding_from_monogon(const traintrack& tt,
                                        int mono,
                                        int dir = 1);

  // DFS emitter used by coding_from_monogon().
  static void recursive_coding(const multigon& mm,
                               int pin,
                               int ein,
                               coding_vec& code,
                               int dir);
};

} // namespace detail
} // namespace traintracks

#endif // TRAINTRACKS_CODING_HPP
