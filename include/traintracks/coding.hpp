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
#include <vector>
#include <jlt/vector.hpp>

namespace traintracks {

class traintrack;
class mathmatrix_permplus1;
class multigon;
class edge;

// Canonical numbering of the prongs and edges of a normalised train track.
//
// Both numberings come from the same depth-first walk that defines the
// coding and the order of weights(): start at monogon 0 (or the monogon
// passed to coding_engine::numbering), cross its edge, and at each
// multigon visit the remaining edges in cycle_edges order, descending into
// every multigon that is not an uncusped monogon.
//
// - Edge e (0-based) is the e-th edge met by the walk, so edge numbers
//   agree with the positions used by weights() and the transition
//   matrices.  Its canonical orientation runs from the multigon that first
//   reaches it (tail) to the multigon it leads to (head).
// - When a multigon is entered at prong pin, its prongs pin, pin+1, ...
//   (mod k, in cycle_edges order) receive consecutive numbers.  Monogon 0
//   has prong number 0.
// - Side q of a multigon is the infinitesimal (or, for a punctured
//   multigon, peripheral) side from prong q to side_to(q), the next prong
//   in cycle_edges order.  For a monogon it is the loop around the
//   puncture.
//
// Letters of a train-track map (jlt::freeauto<int>) are signed: main edge
// e is +-(e+1), side q is +-(nedges()+1+q), a negative sign meaning the
// reverse orientation.  tail_of()/head_of() give the prong numbers at
// which a signed letter starts and ends.
//
// "Label" is deliberately not used here: it is the puncture label of a
// multigon (multigon::label()).
struct ttnumbering
{
  struct prong_info
  {
    int multigon;     // index in the track's multigon list
    int prong;        // prong index within that multigon
    int nprongs;      // number of prongs of that multigon
    bool punctured;
  };

  int start_monogon = 0;
  std::vector<prong_info> prong;                 // by prong number
  std::vector<std::vector<int> > prong_number;   // [multigon][prong]
  std::vector<int> edge_tail;                    // by edge number
  std::vector<int> edge_head;                    // by edge number
  // Identity of the edge objects, by edge number.  Valid only as long as
  // the track is not modified; used to follow edges through a fold.
  std::vector<const edge*> edge_ptr;

  int nprongs() const { return prong.size(); }
  int nedges() const { return edge_tail.size(); }

  // Prong reached from prong q along side q.
  int side_to(const int q) const;

  // Letter classification and endpoints (see above for the encoding).
  bool is_main(const int letter) const;
  bool is_side(const int letter) const;
  int main_letter(const int e) const { return e + 1; }
  int side_letter(const int q) const { return nedges() + 1 + q; }
  int edge_of(const int letter) const;   // 0-based edge number of a main letter
  int side_of(const int letter) const;   // prong number q of a side letter
  int tail_of(const int letter) const;
  int head_of(const int letter) const { return tail_of(-letter); }

  bool operator==(const ttnumbering& o) const;
  bool operator!=(const ttnumbering& o) const { return !operator==(o); }

  std::ostream& print(std::ostream& strm) const;
};

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

  // Canonical prong/edge numbering rooted at uncusped monogon mono.
  // Walks the raw structure; does not require the track to be normalised
  // (the public traintrack::numbering() does, and uses monogon 0).
  static ttnumbering numbering(const traintrack& tt, int mono);

private:
  // DFS used by numbering(): number the prongs of mm on entry, record the
  // head of the entry edge, then number the outgoing edges in
  // cycle_edges order.
  static void recursive_numbering(const traintrack& tt,
                                  const multigon& mm,
                                  int pin,
                                  int ein,
                                  int entry_edge,
                                  ttnumbering& num);

  // Give consecutive numbers to the prongs of multigon mi, starting at pin.
  static void number_prongs(const traintrack& tt, int mi, int pin,
                            ttnumbering& num);

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
