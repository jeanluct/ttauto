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

#ifndef MULTIGON_HPP
#define MULTIGON_HPP

#include <iostream>
#include <cstdlib>
#include <memory>
#include <jlt/vector.hpp>
#include "edge.hpp"
#include "traintracks_util.hpp"

#ifndef TTAUTO_NO_SHARED_PTR
#if __cplusplus > 199711L
#include <memory>
namespace ttauto
{
  using std::shared_ptr;
}
#else
#include <boost/shared_ptr.hpp>
namespace ttauto
{
  using boost::shared_ptr;
}
#endif
#endif


namespace ttauto {

class edge;
class traintrack;

class multigon : public std::enable_shared_from_this<multigon>
{
  static const int debug = 0;

public:
#ifdef TTAUTO_NO_SHARED_PTR
  typedef edge*				edgep;
#else
  typedef shared_ptr<edge>		edgep;
#endif
  typedef jlt::vector<edgep>		egpVec;
  typedef jlt::vector<egpVec>		egpVecVec;
  typedef jlt::vector<int>		intVec;

private:
  int k;			// Number of prongs.
  egpVecVec egv;		// Each prong points to one or more edges.
  int lab;			// Label set.
  bool punct;			// Punctured?

public:

  // Create a k-pronged multigon.
  // Allocate one branch per prong.
  multigon(const int k_ = 1, const int l_ = 0) :
    k(k_), egv(k_,egpVec(1)) , lab(l_), punct(false)
  {
    if (k <= 0)
      {
	std::cerr << "Illegal number of prongs in multigon.\n";
	std::exit(1);
      }

      // Monogons and bigons must be punctured.
    if (k == 1 || k == 2) punct = true;
  }

  // How to copy multigons?  The problem is that if they already point
  // elsewhere, we'll end up with multiple pointers to the same
  // location.  And then when we copy the edge, how do we choose
  // which pointer to update?  Maybe this is where "smart pointers"
  // would be useful.
  //
  // The other option is to code like my previous attempt and
  // recursively copy everything.  In that case only one object would
  // create everything.

  // Copy constructor.
  // Be careful: pointers to edges are duplicated!
  multigon(const multigon& mm) :
    k(mm.k), egv(mm.egv), lab(mm.lab), punct(mm.punct)
  {
    if (debug) std::cerr << "multigon::copy constructor\n";
  }

  // Be careful: pointers to edges are duplicated!
  multigon& operator=(const multigon& mm)
  {
    if (debug) std::cerr << "multigon::operator=\n";

    k = mm.k;
    egv = mm.egv;
    lab = mm.lab;
    punct = mm.punct;

    return *this;
  }

  edgep Edge(const int p, const int e);

  const edgep Edge(const int p, const int e) const;

  // Cycle through edges, incrementing the prong as needed.
  void cycle_edges(int& p, int &e, const int dir = 1) const;

  // The number of prongs of the multigon.
  int prongs() const { return k; }

  // Is the multigon punctured?
  bool punctured() const { return punct; }

  // Puncture the multigon.
  void puncture() { punct = true; }

  // The label of the multigon.
  int label() const { return lab; }

  // Set the label of the multigon.
  void set_label(const int lb) { lab = lb; }

  // Return the total number of exterior cusps in the multigon.
  int cusps() const { return (edges()-k); }

  // Return the number of edges hooked to the pth prong.
  int edges(const int p) const { return egv[p].size(); }

  // Return the total number of edges hooked to the multigon.
  int edges() const;

  double weight(const int p, const int b) const;

  void weight(const int p, const int b, const double w) const;

  bool check() const;

  // Attach a new edge at prong p, branch e.
  void attach_edge(const int p = 0, const int e = 0);

  // Attach an edge eg at prong p, branch e.
  void attach_edge(edgep eg, const int p = 0, const int e = 0);

  // Insert a new edge at prong p, branch e (existing branches are
  // moved out of the way).
  void insert_edge(const int p = 0, const int e = 0);

  // Insert an edge eg at prong p, branch e (existing branches are
  // moved out of the way).
  void insert_edge(edgep eg, const int p = 0, const int e = 0);

  // Is there an edge already attached at prong p, edge e?
  bool prong_edge_is_unattached(const int p, const int e) const;

  // Move all prongs clockwise by i (anticlockwise for i < 0).
  void cycle_prongs(const int i = 1);

  // A vector with the number of edges sequentially from prong p0.
  intVec edge_sequence(const int p0 = 0) const;

  // True if multigons have the same edge structure.
  bool operator==(const multigon& mm) const;

  bool operator!=(const multigon& mm) const { return !(mm == *this); }

  bool operator<(const multigon& mm) const;

  bool operator>(const multigon& mm) const { return (mm < *this); }

  bool operator<=(const multigon& mm) const { return !(mm < *this); }

  bool operator>=(const multigon& mm) const { return !(mm > *this); }

  // Cycle the prongs to maximise the edge sequence norm.
  intVec normalise();

  // Print some information about the multigon.
  std::ostream& print(std::ostream& strm = std::cout) const;

private:
  // Cycle through edge pointers, find the end pointing to pm_old, and
  // update to point to this.
  /* void update_edge_prong_pointers(const multigon *pm_old);*/
  void update_edge_prong_pointers(const std::shared_ptr<const multigon> pm_old);

  bool check_range(const int p, const int e) const
  {
    if (debug && (p < 0 || p >= prongs()))
      {
	std::cerr << "Prong out of range";
	std::cerr << " in multigon::check_range.\n";
	std::exit(1);
      }
    if (debug && (e < 0 || e >= edges(p)))
      {
	std::cerr << "Edge out of range";
	std::cerr << " in multigon::check_range.\n";
	std::exit(1);
      }
    return true;
  }

  void erase_edge_pointer(const int p, const int e);

  void point_to_edge(edgep eg, const int p, const int b);

  friend class edge;
  friend void swap(std::shared_ptr<multigon> m1, std::shared_ptr<multigon> m2);
  // friend void swap(multigon& m1, multigon& m2);
};

// Swap two multigons.
void swap(std::shared_ptr<multigon> m1, std::shared_ptr<multigon> m2);
//void swap(multigon& m1, multigon& m2);


//
// Inline method definitions
//

inline multigon::edgep multigon::Edge(const int p, const int e)
{
  if (debug) check_range(p,e);

  return egv[p][e];
}

inline const multigon::edgep multigon::Edge(const int p, const int e) const
{
  if (debug) check_range(p,e);

  return egv[p][e];
}

inline bool multigon::prong_edge_is_unattached(const int p, const int e) const
{
  if (debug && (p < 0 || p >= prongs()))
    {
      std::cerr << "Prong out of range";
      std::cerr << " in multigon::check_range.\n";
      std::exit(1);
    }
  if (debug && (e < 0 || e >= edges(p)))
    {
      std::cerr << "Edge out of range";
      std::cerr << " in multigon::check_range.\n";
      std::exit(1);
    }

  // If the edge is out of range but valid (positive), it is unattached.
  if (e >= edges(p)) return true;

  return (egv[p][e] == 0);
}

// Cycle through edges, incrementing the prong as needed.
inline void multigon::cycle_edges(int& p, int &e, const int dir) const
{
  if (dir == 1)
    {
      // Cycle clockwise.
      if (++e == edges(p)) { e = 0; p = ttauto::mod(p+1,prongs()); }
    }
  else
    {
      // Cycle anticlockwise.
      // Important to update the prong first!
      if (--e == -1) { p = ttauto::mod(p-1,prongs()); e = edges(p)-1; }
    }
}

inline int multigon::edges() const
{
  int c = 0;
  for (int p = 0; p < k; ++p)
    {
      c += edges(p);
    }
  return c;
}

} // namespace ttauto

#endif // MULTIGON_HPP
