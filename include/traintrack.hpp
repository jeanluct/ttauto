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

#ifndef TRAINTRACK_HPP
#define TRAINTRACK_HPP

#include <memory>
#include <jlt/vector.hpp>
#include <jlt/mathmatrix.hpp>
#include "edge.hpp"
#include "multigon.hpp"
#include "traintracks_util.hpp"
#include "mathmatrix_permplus1.hpp"


namespace ttauto {

class traintrack
{
  static const bool label_multiprongs = true;

public:
  typedef multigon::edgep				edgep;
#if __cplusplus > 199711L && !defined(TTAUTO_NO_SHARED_PTR)
  // There is only one owner for each multigon, so use std::unique_ptr (C++11).
  typedef std::unique_ptr<multigon>			mgonp;
#else
  typedef multigon*					mgonp;
#endif
  typedef jlt::vector<mgonp>				mgpVec;
  typedef multigon::intVec				intVec;
  typedef jlt::vector<double>				dblVec;
  typedef mgpVec::const_iterator			cmit;
  typedef mgpVec::iterator				mit;

  static const int debug = 0;
  static const bool exploit_symmetries = true;

private:

  mgpVec mgv;

public:
  //
  // Constructors
  //

  traintrack() {}

  // Copy constructor.
  traintrack(const traintrack& tt) { copy(*this,tt); }

  // Make a train track from its coding.
  traintrack(const intVec& code);

  // Make a train track from a string of its coding.
  // Only works when the numbers are < 10.
  // Train track must be one-indexed: e.g.:
  //   1111 1311 2311 1111 3312 1111 3322 1111
  // as output by print_coding().
  traintrack(const char* codes);

  // Assignment operator.
  traintrack& operator=(const traintrack& tt)
  {
    mgv.clear(); copy(*this,tt); return *this;
  }


  //
  // Constructors for specific train tracks (defined in ttbuild.cpp)
  //

  // Track with N monogons around punctures and a (N-2)-gon on the
  // boundary (all cusps).
  traintrack(const int N);

  // Track with N monogons around punctures, an unpunctured K-gon,
  // and an (N-K)-gon on the boundary.
  traintrack(const int N, const int K);

  // Track with N monogons around punctures, an unpunctured K1-gon, an
  // unpunctured K2-gon, and an (N-(K1+K2)+2)-gon on the boundary.
  traintrack(const int N, const int K1, const int K2);

  // Track with N monogons around punctures, unpunctured K1-, K2-,
  // and K3-gons, and an (N-(K1+K2+K3)+4)-gon on the boundary.
  traintrack(const int N, const int K1, const int K2, const int K3);

  // Track with N monogons around punctures and unpunctured K-gons
  // specified by a vector Kv.
  traintrack(const int N, const intVec& Kv);

#if 0
  // Track with all multigons specified by Kv, including monogons.
  traintrack(const intVec& Kv);
  /* Cannot do this since it clashes with creating a track from coding. */
  /* Introduce coding class. */
#endif

  //
  // Public methods
  //

  const multigon& Multigon(const int m) const;

  int edges() const { return mgv.size()-1; } /* Count explicitly */

  int multigons() const { return mgv.size(); }

  int monogons() const;

  int punctures() const;

  // Count the total number of exterior cusps in graph.
  // Exterior cusps are the same as boundary prongs.
  int cusps() const;

  // Isotopy of train tracks.
  bool operator==(const traintrack& tt) const;

  // Return the coding, minimised over uncusped monogons, but
  // leave the monogons untouched.
  intVec coding(const int dir = 1) const;

  // Return true if the train track is reflection-symmetric.
  bool is_reflection_symmetric() const { return (coding(1) == coding(-1)); }

  // Return positive (order-1) if the train track is cyclically-symmetric.
  int is_cyclically_symmetric() { return (cyclic_symmetry().order()-1); }

  // Return permutation for cyclically-symmetric train track.
  // The symmetry is Delta^order, where Delta = sigma_1 ... sigma_(n-1).
  mathmatrix_permplus1 cyclic_symmetry();

  void set_label(const int m, const int lb);

  // Give unique label to each punctured multigon; label 0 if unpunctured.
  void pure_braid();

  //
  // Put into normal form.
  //
  // Order multigons according to increasing prongness,
  // then ascending total number of branches,
  // and finally a lexicographical compare of their maximised edge sequence.
  //
  // Then put the uncusped monogon that minimises the train track
  // coding (lexicographically) at the start.
  //
  void normalise();

  // Check that everything is hooked properly.
  void check() const;

  // The maximum number of possible foldings.
  // (There could be fewer due to adjacent cusps.)
  int foldings() const { return 2*cusps(); }

  // Fold at a cusp until the next edge.
  // This combines two edges.
  // Cusps have to be given an ordering.
  //  f even = fold counterclockwise
  //  f odd  = fold clockwise
  bool fold(const int f);

  // Fold and find transition matrix.
  jlt::mathmatrix<int> fold_transition_matrix(const int f)
  {
    jlt::mathmatrix<int> M(ttauto::fold_transition_matrix(*this,f));
    fold(f);
    return M;
  }

  // Return vector of edge weights.
  dblVec weights(const int mono = 0) const;

  // Set edge weights from iterator.
  dblVec::const_iterator weights(dblVec::const_iterator wi);

  // Print some information about the traintrack.
  std::ostream& print(std::ostream& strm = std::cout) const;

  // Print singularity data of the train track.
  std::ostream& print_singularity_data(std::ostream& strm) const;

  // Print in a format that can be used by Mathematica.
  std::ostream& printMathematicaForm(std::ostream& strm = std::cout) const;

  // Print coding.
  std::ostream& print_coding(std::ostream& strm = std::cout,
			     const int dir = 1) const;

private:
  //
  // Helper methods for building tracks (defined in ttbuild.cpp)
  //

  // Track with N monogons around punctures and a (N-2)-gon on the
  // boundary (all cusps).
  void ttbuild_all_monogons(const int N);

  // Track with N monogons around punctures, an unpunctured K-gon,
  // and an (N-K)-gon on the boundary.
  void ttbuild_monogons_and_multigon(const int N, const int K);

  // Track with N monogons around punctures, an unpunctured K1-gon, an
  // unpunctured K2-gon, and an (N-(K1+K2)+2)-gon on the boundary.
  void ttbuild_monogons_and_two_multigons(const int N,
					  const int K1, const int K2);

  // Track with N monogons around punctures, unpunctured K1-, K2-,
  // and K3-gons, and an (N-(K1+K2+K3)+4)-gon on the boundary.
  void ttbuild_monogons_and_three_multigons(const int N, const int K1,
					    const int K2, const int K3);

  // Track with N monogons around punctures and unpunctured K-gons
  // specified by a vector Kv.
  void ttbuild_monogons_and_multigons(const int N,
				      const jlt::vector<int>& Kv);

  // Track with N monogons around punctures and unpunctured K-gons
  // specified by a vector Kv.
  //
  //   (All the multigons are anchored at a monogon, so not all tracks
  //    can be made this way.)
  void ttbuild_monogoncusps(const int N, const intVec& Kv);

  // Add a monogon to all the free edge endings.
  void capoff();

  // Replace an uncusped monogon at position L in mgv by a K-pronged multigon.
  void monogon_to_multigon(const int L, const int K);

  //
  // Other helper functions and private members
  //

  // non-const reference to multigon is private.
  multigon& Multigon(const int m);

  // The coding of a train track is a sequence of coding_blocks.
  struct coding_block;

  void copy(traintrack& ttnew, const traintrack& ttexist);

  // Find the index of a multigon edge in the vector, given a pointer.
  int multigon_index(const multigon* mm) const;

  // Compute an index of for a multigon prong.
  int multigon_prong_index(const int mi, const int pi) const;

  // Do two normalised tracks have the same multigons?
  bool same_multigons(const traintrack& tt) const;

  void recursive_build(edgep& ee, intVec::const_iterator& cd);

  void recursive_coding(const multigon& mm, const int pp, const int ee,
			intVec& code, const int dir) const;

  void recursive_get_weights(const multigon& mm, const int pp, const int ee,
			     dblVec& wv) const;

  void recursive_set_weights(const multigon& mm, const int pp, const int ee,
			     dblVec::const_iterator& wi);

  bool recursive_find_cusp(multigon& mm, const int pp, const int ee,
			   int& fcusp, multigon*& mmc, int& pc, int& ec) const;

  // Minimise coding over uncusped monogons, such that the min is
  // obtained from the first position.
  intVec minimise_coding();

  // Return a coding vector for the train track, which is used to test
  // for equality.  The coding vector starts at a given monogon.
  intVec coding_from_monogon(const int mono, const int dir = 1) const;

  // Sort ascending using the strict order relation for multigons.
  void sort();

  // Swap the position of two multigons in the vector.
  void swap(const int m1, const int m2);

  // Fold cusp c of prong p of multigon m in direction dir.
  //   dir = 1 clockwise, dir = 1 anticlockwise.
  bool fold(multigon& mm, const int p, const int c, const int dir);

  // Print a coding block.
  friend std::ostream& operator<<(std::ostream& strm, const coding_block& b);
};


//
// Helper function for building tracks (defined in ttbuild.cpp)
//

jlt::vector<traintrack> ttbuild_list2(const int N);

jlt::vector<traintrack> ttbuild_list(const int N, const int N2 = 0);

// Comparison function for sorting strata.
bool compare_strata(const jlt::vector<int>& first,
		    const jlt::vector<int>& second);

// How many prongs on the boundary for singularity data sdata?
int boundary_prongs(jlt::vector<int> sdata);

// How many prongs on the boundary for N punctures and multigons given
// by the vector K?
int boundary_prongs(const int N, const jlt::vector<int> K);


//
// Class member definitions
//

// The coding of a train track is a sequence of coding_blocks.
struct traintrack::coding_block
{
  static const int length = 5;

  int prong;
  int nprongs;
  int label;
  int edge;
  int nedges;

  // Default is uncusped monogon with label 0.
  coding_block(const int p = 0, const int np = 1, const int lb = 0,
	       const int e = 0, const int ne = 1)
    : prong(p), nprongs(np), label(lb), edge(e), nedges(ne) {}

  coding_block(intVec::const_iterator& ci)
    : prong(*ci++), nprongs(*ci++), label(*ci++), edge(*ci++), nedges(*ci++) {}

  void append_to(intVec& v) const
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


//
// Inline method definitions
//

inline multigon& traintrack::Multigon(const int m)
{
  if (debug && (m < 0 || m >= (int)mgv.size()))
    {
      std::cerr << "Nonexistent multigon " << m;
      std::cerr << " in traintrack::Multigon\n";
      std::exit(1);
    }
  return *mgv[m];
}

inline const multigon& traintrack::Multigon(const int m) const
{
  if (debug && (m < 0 || m >= (int)mgv.size()))
    {
      std::cerr << "Nonexistent multigon " << m;
      std::cerr << " in traintrack::Multigon\n";
      std::exit(1);
    }
  return *mgv[m];
}

inline int traintrack::monogons() const
{
  int nm = 0;
  for (int m = 0; m < multigons(); ++m)
    {
      if (Multigon(m).prongs() == 1) ++nm;
    }
  return nm;
}

inline int traintrack::punctures() const
{
  int np = 0;
  for (int m = 0; m < multigons(); ++m)
    {
      if (Multigon(m).punctured()) ++np;
    }
  return np;
}

inline int traintrack::cusps() const
{
  int c = 0;

  for (cmit it = mgv.begin(); it != mgv.end(); ++it)
    {
      c += (*it)->cusps();
    }
  return c;
}

// Isotopy of train tracks.
// Use coding to decide equality.
inline bool traintrack::operator==(const traintrack& tt) const
{
  // A minimum requirement is to have the same multigon structure,
  // with the same number of edges hooked to each prongs.
  if (!same_multigons(tt)) return false;

  return (coding() == tt.coding());
}

// Compare two train tracks on the basis of multigons only.
// Both tracks must be normalised.
inline bool traintrack::same_multigons(const traintrack& tt) const
{
  if (edges() != tt.edges() || multigons() != tt.multigons()) return false;

  for (int m = 0; m < multigons(); ++m)
    {
      if (*mgv[m] != *tt.mgv[m]) return false;
    }

  return true;
}

// Put into normal form.
inline void traintrack::normalise()
{
  for (int m = 0; m < multigons(); ++m)
    {
      Multigon(m).normalise();
    }
  sort();
  minimise_coding();
}


//
// Friend functions
//

// Swap the position of two multigons in the vector.
inline void traintrack::swap(const int m1, const int m2)
{
  if (debug && (m1 < 0 || m1 >= multigons() || m2 < 0 || m2 >= multigons()))
    {
      std::cerr << "Nonexistent multigon in traintrack::swap.\n";
      std::exit(1);
    }
  ttauto::swap(*mgv[m1],*mgv[m2]);
}

inline std::ostream&
operator<<(std::ostream& strm, const traintrack::coding_block& b)
{
  // Print prong/label/edge with 1 offset rather than 0 offset.
  if (traintrack::label_multiprongs)
    strm << b.prong+1 << b.nprongs << b.label+1 << b.edge+1 << b.nedges;
  else
    strm << b.prong+1 << b.nprongs << b.edge+1 << b.nedges;
  return strm;
}

} // namespace ttauto

#endif // TRAINTRACK_HPP
