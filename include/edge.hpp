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

#ifndef EDGE_HPP
#define EDGE_HPP

#include <iostream>
#include <cstdlib>
#include <jlt/vector.hpp>
#include "multigon.hpp"

namespace ttauto {

class multigon;
class traintrack;

class edge
{
  static const int debug = 0;
  static const int nends = 2;	// An edge has two ends.

public:
  typedef jlt::vector<multigon*>	mgpVec;
  typedef multigon::intVec		intVec;

private:
  double wt;			// Weight on edge.
  mgpVec mg;			// The multigon at each end of the edge.
  intVec pr;			// The mother prong of each multigon.
  intVec pre;			// Is this the 0th, 1st, etc. edge

public:
  edge(const double wt_ = 0) : wt(wt_), mg(nends), pr(nends), pre(nends) {}

  // Copy constructor.
  edge(const edge& eg) { copy(*this,eg); }

  edge& operator=(const edge& eg) { copy(*this,eg); return *this; }

#if 0
  /* Unused.  See definition to see why it's maybe a bad idea. */
  void attach_to_multigon(multigon* mm, const int p = 0, const int e = 0);
#endif

  // Is the edge ending en unattached?
  bool is_unattached(const int en) const { return (mg[en] == 0); }

  // True if either edge ending is unattached.
  bool is_unattached() const { return (mg[0] == 0 || mg[1] == 0); }

  // Get the weight on the edge.
  double weight() const { return wt; }

  // Set the weight on the edge.
  void weight(double wt_) { wt = wt_; }

  // Check for consistency of pointers.
  bool check() const;

  // Which edge ending is a multigon attached to?
  int which_ending(const multigon* mm) const;

  // Detach from a given multigon.
  // class traintrack's fold method needs this.
  void detach_from_multigon(const multigon* mm);

  // Return the pointer, prong, and edge to a multigon that the edge
  // is hooked to from mm.  Needed by class traintrack.
  multigon* target_multigon(const multigon* mm, int& t_pr, int& t_pre) const;

private:
  // Be careful: pointers to edges are duplicated!
  void copy(edge& enew, const edge& eold)
  {
    if (debug) std::cerr << "edge::copy\n";

    enew.wt = eold.wt;
    enew.mg = eold.mg;
    enew.pr = eold.pr;
    enew.pre = eold.pre;
  }

  // Are the two edge endings attached to the same multigon?
  //
  // This is one of the few methods that explicitly assumes a branch
  // has 2 ends (doesn't respect edge::nends).
  bool attached_to_same_multigon() const
  { return (mg[0] == mg[1]); }

  void swap_endings()
  {
    // This is one of the few methods that explicitly assumes a branch
    // has 2 ends (doesn't respect edge::nends).
    std::swap(mg[0],mg[1]);
    std::swap(pr[0],pr[1]);
    std::swap(pre[0],pre[1]);
  }

  // Detach from both multigons completely.  Keeps only the weight.
  void detach_from_multigons()
  {
    for (int en = 0; en < nends; ++en)
      {
	// Let mother multigon know we've detached.
	mg[en]->erase_edge_pointer(pr[en],pre[en]);
	// Zero everything.
	mg[en] = 0;
	pr[en] = 0;
	pre[en] = 0;
      }
  }

  void point_to_multigon(multigon* mm, const int p, const int e);

  friend class multigon;
  // friend class traintrack; // Only for detach_from_multigon and
			   // target_multigon.
  friend void swap(multigon* m1, multigon* m2);
};


//
// Method definitions
//

#if 0
/* Unused. */
inline void edge::attach_to_multigon(multigon* mm, const int p, const int e)
{
  point_to_multigon(mm,p,e);

  /* If we decide we eventually need this method, we can't do this:
   edgep is now a shared_ptr, so creating from the this pointer is a
   bad idea.  If we do need tto do this, we have to derive edge from
   std::enable_shared_from_this<edge>.  See
   http://en.cppreference.com/w/cpp/memory/enable_shared_from_this */
  multigon::edgep ep(this);
  /* Use instead: */
  // multigon::edgep ep(shared_from_this());

  // Make sure the multigon knows the edge is attached.
  mm->point_to_edge(ep,p,e);
}
#endif

inline void edge::point_to_multigon(multigon* mm, const int p, const int e)
{
  // Attach to the first unattached ending.
  for (int en = 0; en < nends; ++en)
    {
      if (is_unattached(en))
	{
	  mg[en] = mm;
	  pr[en] = p;
	  pre[en] = e;
	  return;
	}
    }
  std::cerr << "No unattached edge ending in edge::point_to_multigon.\n";
  std::exit(1);
}

inline int edge::which_ending(const multigon* mm) const
{
  for (int en = 0; en < nends; ++en)
    {
      if (mg[en] == mm) return en;
    }
  std::cerr << "Could not find pointer back to prong";
  std::cerr << " in edge::which_ending.\n";
  std::exit(1);
}

inline void edge::detach_from_multigon(const multigon* mm)
{
  int en = which_ending(mm);
  // Let mother multigon know we've detached.
  mg[en]->erase_edge_pointer(pr[en],pre[en]);
  // Zero everything on that ending.
  mg[en] = 0;
  pr[en] = 0;
  pre[en] = 0;
}

// Return the pointer, prong, and edge to a multigon that the edge
// is hooked to from mm.
inline multigon* edge::target_multigon(const multigon* mm,
				       int& t_pr, int& t_pre) const
{
  // Target ending is the one mm is not on.
  int t_en = 1 - which_ending(mm);
  t_pr = pr[t_en];
  t_pre = pre[t_en];
    return mg[t_en];
}

inline bool edge::check() const
{
  for (int en = 0; en < nends; ++en)
    {
      if (is_unattached(en))
	{
	  std::cerr << "Unhooked edge in edge::check.\n";
	  std::exit(1);
	}
#ifdef TTAUTO_NO_SHARED_PTR
      if (mg[en]->egv[pr[en]][pre[en]] != this)
#else
      if (mg[en]->egv[pr[en]][pre[en]].get() != this)
#endif
	{
	  std::cerr << "Inconsistent pointer in edge::check.\n";
	  std::exit(1);
	}
    }
  if (attached_to_same_multigon())
    {
      std::cerr << "Both ends hooked to same multigon in edge::check.\n";
      std::exit(1);
    }
  return true;
}

} // namespace ttauto

#endif // EDGE_HPP
