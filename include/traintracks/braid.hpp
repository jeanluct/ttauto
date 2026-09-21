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

#ifndef TRAINTRACKS_BRAID_HPP
#define TRAINTRACKS_BRAID_HPP

#include <iosfwd>
#include <vector>
#include "traintracks/embedding.hpp"
#include "traintracks/fold_map.hpp"

namespace traintracks {

// A word in the braid group B_n, as signed generator indices: +i is the
// right-handed half twist sigma_i of the punctures at positions i and i+1,
// -i its inverse.  Read left to right as a motion in time.
class braidword
{
public:
  braidword(const int n_ = 0) : n(n_) {}
  braidword(const int n_, const std::vector<int>& w_) : n(n_), w(w_) {}

  int strings() const { return n; }
  int length() const { return (int)w.size(); }
  const std::vector<int>& word() const { return w; }

  // Sum of the exponents, which is invariant under conjugation.
  int exponent_sum() const;

  // Where each puncture ends up: entry i-1 is the final position of the
  // puncture that started at position i.
  std::vector<int> permutation() const;

  braidword& operator*=(const braidword& b);
  braidword inverse() const;

  // Free reduction of adjacent inverse pairs.  Nothing stronger: this is
  // not a normal form.
  braidword& reduce();

  // Growth of the braid, from the action on Dynnikov coordinates: iterate
  // the braid on a loop and measure how fast it grows.  For a
  // pseudo-Anosov braid this is the dilatation, to about ten digits.  When
  // the growth is not exponential the loop grows polynomially instead and
  // the iteration only approaches 1 slowly from above, so read a value
  // near 1 as "not pseudo-Anosov" and no more.  This shares nothing with
  // the rest of the library, which is what makes it worth checking a braid
  // read off a folding path against.
  double growth() const;

  // Block A = [i, i+a-1] moves right past block B = [i+a, i+a+b-1], every
  // crossing of the given sign.  a == 1 is the braid
  // beta_{i,j} = sigma_i ... sigma_{j-1} of Lanneau and Thiffeault,
  // Prop. 4, with j = i+b; b == 1 is the same motion in reverse.  The
  // general case, with both blocks longer than one puncture, does occur.
  static braidword block_swap(const int n, const int i, const int a,
                              const int b, const int sign);

  // Rotation by one position: the puncture at position k moves to k-1, and
  // the one at position 1 wraps round to n.  sign gives the handedness.
  // delta(n,1)^n is the full twist, which generates the centre.
  static braidword delta(const int n, const int sign);

  std::ostream& print(std::ostream& strm) const;
  std::ostream& printMathematicaForm(std::ostream& strm) const;

private:
  int n;
  std::vector<int> w;
};

inline braidword operator*(braidword a, const braidword& b) { return a *= b; }

// The motion of the punctures that one fold asks for.
//
// A fold slides the end of one edge along another to the target multigon
// and round one corner of it.  When the target is an unpunctured polygon
// the new route crosses only infinitesimal edges and the track stays
// properly embedded, so nothing moves.  When it is a punctured monogon the
// edge goes once around the puncture, crossing the segment that joins it to
// the boundary, and the punctures have to be moved to undo that
// (Lanneau and Thiffeault, Prop. 4).
//
// What moves is read off the track before the fold.  In the boundary walk
// the punctures S of the subtree the fold moves are consecutive, since the
// multigons and main edges form a tree, and the target puncture T sits
// right next to them: after S when the fold runs clockwise, before S when
// it runs anticlockwise, because the fold is only legal when the edge it
// folds onto is the outermost one at the target.  So the motion is always
// that S and T change places.
struct fold_swap
{
  int first = 0;        // position of the leftmost puncture that moves
  int nleft = 0;        // size of the left block
  int nright = 0;       // size of the right block
  int dir = 0;          // fold direction, copied from fold_map_data
  int rotate = 0;       // see below
  bool trivial = true;

  // The dart to cut the folded track at, so that the punctures keep the
  // positions this swap leaves them in.  Feed it to outer_embedding.
  int next_cut_dart = 0;

  // The two blocks are next to each other around the disc, but the cut
  // that turns that cyclic order into positions 1..n can fall between
  // them.  rotate is how many places the cut has to move so that it does
  // not, and the braid then starts with that many copies of delta.  The
  // other fields are in the rotated frame.
};

// What one fold does to the punctures, read off the track before it.
// `before` must be the embedding of the track the fold was applied to.
// Exits if the punctures the fold moves are not consecutive, which would
// mean the embedding and the fold record disagree.
fold_swap fold_block_swap(const ttnumbering& before_num,
                          const fold_map_data& fm,
                          const tt_embedding& before);

// k copies of the rotation that moves the cut one place, in the handedness
// the fold braids use.  Also what the final identification of a closed
// folding path needs.
braidword rotation_braid(const int n, const int k);

// The same, as a braid word: rotate copies of delta, then the swap.
braidword fold_braid(const ttnumbering& before_num, const fold_map_data& fm,
                     const tt_embedding& before);

} // namespace traintracks

#endif // TRAINTRACKS_BRAID_HPP
