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

#include <iostream>
#include <algorithm>
#include <set>
#include "traintracks/edge.hpp"
#include "traintracks/multigon.hpp"
#include "traintracks/traintrack.hpp"

namespace traintracks {

// Make a train track from its coding.
traintrack::traintrack(const traintrack::intVec& code)
{
  isnormalised = false;

  // Train track always starts with an uncusped monogon.
#if __cplusplus > 201103L && !defined(TRAINTRACKS_NO_SHARED_PTR)
  // C++14 has make_unique.
  mgv.push_back(mgonp(std::make_unique<multigon>(1)));
#else
  mgv.push_back(mgonp(new multigon(1)));
#endif
  Multigon(0).attach_edge();
  if (label_multigons)
    Multigon(0).set_label(code[2]); // Copy the label of first monogon.

  // Iterator for coding: skip initial uncusped monogon marker.
  intVec::const_iterator ci(code.begin() + detail::coding_block::length);

  // Recurse down and build train track.
  recursive_build(Multigon(0).Edge(0,0),ci,code.end());

  normalise();
}

// Make a train track from a string of its coding, in either of the two
// printed forms (see parse_coding in coding.hpp).
traintrack::traintrack(const char* codes) : traintrack(parse_coding(codes))
{
}

// Recursively consume coding blocks and attach descendant multigons.
void traintrack::recursive_build(traintrack::edgep& ee,
				 traintrack::intVec::const_iterator& ci,
				 const traintrack::intVec::const_iterator& cend)
{
  detail::coding_block outb;

  // A truncated coding would otherwise be walked off the end.
  if (cend - ci < detail::coding_block::length)
    {
      std::cerr << "Error in traintrack::recursive_build(): ";
      std::cerr << "coding ends in the middle of the walk.\n";
      exit(1);
    }

  // Extract the block in the coding corresponding to the entry edge.
  // A block consists of detail::coding_block::length digits.
  // The entry edge info is not strictly needed as it can be deduced,
  // but makes things easier.
  detail::coding_block inb(ci);

  if (debug)
    {
      std::cerr << "In: " << inb.nprongs << inb.prong << inb.nedges;
      std::cerr	<< inb.edge << std::endl;
    }

  // Add a new multigon, copying label.
  // Could replace new by std::make_unique<multigon> in C++14.
#if __cplusplus > 201103L && !defined(TRAINTRACKS_NO_SHARED_PTR)
  // C++14 has make_unique.
  mgv.push_back(mgonp(std::make_unique<multigon>(inb.nprongs, inb.label)));
#else
  mgv.push_back(mgonp(new multigon(inb.nprongs, inb.label)));
#endif

  // Attach the multigon to the edge we came down.
  mgv.back()->attach_edge(ee,inb.prong,inb.edge);

  // If it's an uncusped monogon, we're done.
  if (inb.nprongs == 1 && inb.nedges == 1) return;

  // Remember the index of the current (in) multigon.
  int in_m = mgv.size()-1;

  do
    {
      // Extract the block in the coding corresponding to the next exit edge.
      // A block consists of detail::coding_block::length digits.
      if (cend - ci < detail::coding_block::length)
	{
	  std::cerr << "Error in traintrack::recursive_build(): ";
	  std::cerr << "coding ends in the middle of the walk.\n";
	  exit(1);
	}
      outb = detail::coding_block(ci);

      if (debug)
	{
	  std::cerr << "Out: " << outb.nprongs << outb.prong << outb.nedges;
	  std::cerr << outb.edge << std::endl;
	}

      mgv[in_m]->attach_edge(outb.prong,outb.edge);
      // Recurse!
      recursive_build(mgv[in_m]->Edge(outb.prong,outb.edge),ci,cend);

      // Which edge is next?
      if (++outb.edge == outb.nedges)
	{
	  outb.edge = 0;
	  outb.prong = (outb.prong+1) % outb.nprongs;
	}
    }
  while (!(outb.prong == inb.prong && outb.edge == inb.edge));
}

// Copy a traintrack ttexist onto ttnew.
void traintrack::copy(traintrack& ttnew, const traintrack& ttexist)
{
  ttnew.isnormalised = ttexist.isnormalised;

  // Copy the multigons and create new edges.

  // Have to carefully hook the new edges in the proper places.

  for (int m = 0; m < ttexist.multigons(); ++m)
    {
      // Allocate the new multigons.
#if __cplusplus > 201103L && !defined(TRAINTRACKS_NO_SHARED_PTR)
      // C++14 has make_unique.
      ttnew.mgv.push_back(mgonp(std::make_unique<multigon>
				(ttexist.Multigon(m).prongs(),
				 ttexist.Multigon(m).label())));
#else
      ttnew.mgv.push_back(mgonp(new multigon(ttexist.Multigon(m).prongs(),
					     ttexist.Multigon(m).label())));
#endif
      for (int p = 0; p < ttexist.Multigon(m).prongs(); ++p)
	{
	  for (int e = 0; e < ttexist.Multigon(m).edges(p); ++e)
	    {
	      // Find where the edge is hooked to.
	      const multigon* mp = &ttexist.Multigon(m);

	      // Unfortunately from the old edge we only get a pointer
	      // to the multigon.  We need to find what index number
	      // that corresponds to, so we can use it for the new
	      // multigons.
	      int t_pr, t_pre;
	      multigon* t_mp =
		mp->Edge(p,e)->target_multigon(mp,t_pr,t_pre);

	      // Search through the multigons for the pointer.
	      int t_m = 0;
	      for (; t_m < (int)ttexist.multigons(); ++t_m)
		{
#if __cplusplus > 199711L && !defined(TRAINTRACKS_NO_SHARED_PTR)
		  if (ttexist.mgv[t_m].get() == t_mp) break;
#else
		  if (ttexist.mgv[t_m] == t_mp) break;
#endif
		}
	      if (t_m == (int)ttexist.multigons())
		{
		  std::cerr << "Couldn't find multigon";
		  std::cerr << " in traintrack::traintrack::copy.\n";
		  std::exit(1);
		  break;
		}

	      // Attach to the target edge.
	      if (t_m > m)
		{
		  // We haven't allocated the target multigon yet, so
		  // just create a new edge.
		  ttnew.Multigon(m).attach_edge(p,e);
		}
	      else
		{
		  // We've already been there, so there is a dangling
		  // edge waiting for us.
		  edgep ep = ttnew.Multigon(t_m).Edge(t_pr,t_pre);
		  ttnew.Multigon(m).attach_edge(ep,p,e);
		}
	    }
	}
    }
}

void traintrack::require_normalised(const char* where) const
{
  if (isnormalised) return;

  std::cerr << "traintrack must be normalised in " << where << ".\n";
  std::exit(1);
}

// Check that everything is hooked
void traintrack::check() const
{
  // Should also check that no other multigons are created in the
  // graph (simply connected).  Not so easy!  But I guess the
  // recursive methods won't terminate in that case.  Calculate a
  // maximum recursion depth?  Write a recursive_check.

  // Check that there are no unhooked prongs.
  for (cmit it = mgv.begin(); it != mgv.end(); ++it)
    {
      (*it)->check();
    }
}

void traintrack::set_label(const int m, const int lb)
{
  isnormalised = false;

  // If this is called explicitly, then we are labeling multigons.
  if (label_multigons)
    mgv[m]->set_label(lb);
  else
    {
      std::cerr << "Error in traintrack::traintrack::set_label(): ";
      std::cerr << "flag label_multigons must be set.\n";
      exit(1);
    }
  normalise();
}

// Give unique label to each punctured multigon; label 0 if unpunctured.
void traintrack::pure_braid()
{
  isnormalised = false;

  // If this is called explicitly, then we are labeling multigons.
  if (!label_multigons)
    {
      std::cerr << "Error in traintrack::traintrack::pure_braid(): ";
      std::cerr << "flag label_multigons must be set.\n";
      exit(1);
    }

  // Loop over punctured multigons, labeling as we go.
  int lb = 0;
  for (int m = 0; m < (int)mgv.size(); ++m)
    {
      if (mgv[m]->punctured())
	mgv[m]->set_label(++lb);
      else
	mgv[m]->set_label(0);
    }

  normalise();
}

// Sort ascending using the strict order relation for multigons.
void traintrack::sort()
{
  bool swapped;
  do {
    swapped = false;
    for (int m = 0; m < (int)mgv.size()-1; ++m)
      {
	// Compare multigons.
	if (*mgv[m] > *mgv[m+1])
	  {
	    traintrack::swap(m,m+1); swapped = true;
	  }
      }
  } while (swapped);
}

//  f even = fold clockwise
//  f odd  = fold counterclockwise
//
// Cusps are numbered by the canonical numbering (ttnumbering::cusp), i.e.
// in the order of the depth-first walk from monogon 0 that also numbers
// edges and prongs.
bool traintrack::fold(const int f)
{
  // Need normalised train track.
  require_normalised("traintrack::fold");

  if (f < 0 || f >= foldings())
    {
      std::cerr << "Nonexistent cusp in traintrack::traintrack::fold.\n";
      std::exit(1);
    }

  const ttnumbering num = detail::coding_engine::numbering(*this,0);
  int m, p, slot;
  num.fold_cusp(f,m,p,slot);
  return fold(Multigon(m),p,slot,1 - 2*(f % 2));
}

// Resolve fold index f to the concrete cusp location on the current track.
void traintrack::fold_cusp_location(const int f, multigon*& mmc, int& pc, int& ec) const
{
  require_normalised("traintrack::fold_cusp_location");

  if (f < 0 || f >= foldings())
    {
      std::cerr << "Illegal folding index in traintrack::traintrack::fold_cusp_location.\n";
      std::exit(1);
    }

  const ttnumbering num = detail::coding_engine::numbering(*this,0);
  int m;
  num.fold_cusp(f,m,pc,ec);
  mmc = &*mgv[m];
}

// Fold cusp c of prong p of multigon m in direction dir.
//   dir = 1 clockwise, dir = -1 anticlockwise.
//
// A cusp is specified by c, the first of its two edges encountered
// clockwise.
bool traintrack::fold(multigon& mm, const int p, const int c, const int dir)
{
  int e0, e1;

  if (abs(dir) != 1)
    {
      std::cerr << "Bad folding direction in traintrack::traintrack::fold.\n";
      std::exit(1);
    }

  if (c > mm.edges(p)-2)
    {
      std::cerr << "Not a cusp in traintrack::traintrack::fold.\n";
      std::exit(1);
    }

  if (dir == 1)
    { e0 = c; e1 = c+1; }
  else
    { e0 = c+1; e1 = c; }

  // Find target of edge: multigon pointer, prong, and edge.
  int t_pr, t_pre;
  multigon *t_mm = mm.Edge(p,e1)->target_multigon(&mm,t_pr,t_pre);

  if (debug)
    {
      std::cerr << "traintrack::fold dir=" << dir;
      std::cerr << " t_prongs=" << t_mm->prongs() << "\n";
    }

  // Can only fold if target edge is the last (first) of its prong.
  if (t_pre == (dir == 1 ? t_mm->edges(t_pr)-1 : 0))
    {
      // Fold e0 onto e1.  This means that e0 ends up on the next
      // multigon, and adds its weight to e1.
      //
      // We fold until the next prong clockwise (anticlockwise) on the
      // target multigon.
      // The new weight on e1.
      double w0 = mm.Edge(p,e0)->weight();
      double w1 = mm.Edge(p,e1)->weight();
      mm.Edge(p,e1)->weight(w0 + w1);

      if (debug) std::cerr << "traintrack::fold t_pr=" << t_pr;
      // Find the prong we're folding until, which means adding dir to
      // the target.  Then take the mod.
      //
      // Do not use the % function here!!
      int t2_pr = traintracks::mod(t_pr+dir,t_mm->prongs());
      if (debug) std::cerr << " t2_pr=" << t2_pr;

      // e0 will become the new first (last) prong at t2_pr.
      int t2_pre = (dir == 1 ? 0 : t_mm->edges(t2_pr));
      if (debug) std::cerr << " t2_pre=" << t2_pre << std::endl;
      // Now detach edge e0 from mm and re-attach to *t_mm.
      edgep ep = mm.Edge(p,e0);
      ep->detach_from_multigon(&mm);
      t_mm->insert_edge(ep,t2_pr,t2_pre);
    }
  else
    {
      if (debug)
	{
	  std::cerr << "Can't fold " << (dir == 1 ? "" : "anti");
	  std::cerr << "clockwise in traintrack::traintrack::fold.\n";
	}
      return false;
    }

  // We've mucked up the normalisation, so re-normalise.
  normalise();
  return true;
}

// Weights in canonical edge order (ttnumbering rooted at monogon mono).
// Assumes a track is normalised.
traintrack::dblVec traintrack::weights(const int mono) const
{
  require_normalised("traintrack::weights(get)");

  const ttnumbering num = detail::coding_engine::numbering(*this,mono);
  dblVec wv(num.nedges());
  for (int e = 0; e < num.nedges(); ++e) wv[e] = num.edge_ptr[e]->weight();
  return wv;
}

// Set weights from an iterator, in canonical edge order (monogon 0).
// Assumes a track is normalised.
traintrack::dblVec::const_iterator
traintrack::weights(traintrack::dblVec::const_iterator wi)
{
  require_normalised("traintrack::weights(set)");

  const ttnumbering num = detail::coding_engine::numbering(*this,0);
  for (int e = 0; e < num.nedges(); ++e)
    {
      // Reach the edge through its tail slot, which is the only non-const
      // handle on it.
      const int q = num.edge_tail[e];
      const int m = num.prong[q].multigon, p = num.prong[q].prong;
      const std::vector<int>& letters = num.prong_letters[q];
      int slot = -1;
      for (std::size_t k = 0; k < letters.size(); ++k)
	if (letters[k] == num.main_letter(e)) slot = k;
      if (slot < 0)
	{
	  std::cerr << "Edge not found at its tail prong in traintrack::traintrack::weights(set).\n";
	  std::exit(1);
	}
      Multigon(m).Edge(p,slot)->weight(*wi++);
    }
  return wi;
}

// Print some information about the traintrack.
std::ostream& traintrack::print(std::ostream& strm) const
{
  for (int m = 0; m < (int)mgv.size(); ++m)
    {
      strm << "multigon " << m << " is a ";
      mgv[m]->print(strm);
    }
  return strm;
}

// Print singularity data of the train track.
std::ostream& traintrack::print_singularity_data(std::ostream& strm) const
{
  for (int m = 0; m < (int)mgv.size(); ++m)
    {
      strm << mgv[m]->prongs();
      if (mgv[m]->punctured()) strm << ".";
      strm << " ";
    }
  strm << "(" << cusps() << ")";
  return strm;
}

// Print train track in Mathematica-friendly graph-edge format.
std::ostream& printMathematicaForm(std::ostream& strm,
				   const traintrack& tt)
{
  std::set<std::pair<int,int> > ttg;

  for (int m = 0; m < tt.multigons(); ++m)
    {
      for (int p = 0; p < tt.Multigon(m).prongs(); ++p)
	{
	  for (int e = 0; e < tt.Multigon(m).edges(p); ++e)
	    {
	      int t_p, t_e;
	      multigon *t_mm =
		tt.Multigon(m).Edge(p,e)->target_multigon(&tt.Multigon(m),t_p,t_e);
	      int t_m = tt.multigon_index(t_mm);
	      int v1 = tt.multigon_prong_index(m,p);
	      int v2 = tt.multigon_prong_index(t_m,t_p);
	      // The problem is that this does not record the
	      // clockwise order of the edges.
	      ttg.insert(std::make_pair(std::min(v1,v2),std::max(v1,v2)));
	    }
	}
    }
  // Draw multigons as closed polygons in the graph.
  for (int m = 0; m < tt.multigons(); ++m)
    {
      int np = tt.Multigon(m).prongs();
      if (np > 1)
	{
	  for (int p = 0; p < np; ++p)
	    {
	      int v1 = tt.multigon_prong_index(m,p);
	      int v2 = tt.multigon_prong_index(m,(p+1)%np);
	      ttg.insert(std::make_pair(std::min(v1,v2),std::max(v1,v2)));
	    }
	}
    }
  strm << "{";
  for (auto i = ttg.begin(); i != ttg.end(); ++i)
    {
      strm << i->first+1 << "->" << i->second+1;
      if (std::distance(i,ttg.end()) > 1) strm << ",";
    }
  strm << "}";
  return strm;
}

} // namespace traintracks
