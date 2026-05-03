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
  if (label_multiprongs)
    Multigon(0).set_label(code[2]); // Copy the label of first monogon.

  // Iterator for coding: skip initial uncusped monogon marker.
  intVec::const_iterator ci(code.begin() + detail::coding_block::length);

  // Recurse down and build train track.
  recursive_build(Multigon(0).Edge(0,0),ci);

  normalise();
}

// Make a train track from a string of its coding.
// Only works when the numbers are < 10.
// Train track must be one-indexed: e.g.:
//   1111 1311 2311 1111 3312 1111 3322 1111
// as output by print_coding().
traintrack::traintrack(const char* codes)
{
  isnormalised = false;

  std::cerr << "traintrack(char *) constructor: Broken?: ";
  std::cerr << "Reads in an \"old style\" (unlabeled) coding.";
  exit(1);

  intVec code;

  int i = 0;
  while (codes[i] != '\0')
    {
      // Make a number from 0 to 8 from char.
      int num = (int)codes[i]-49;
      // Skip characters that are not numbers.
      if (num >= 0 && num <= 8) code.push_back(num);
      ++i;
    }
  // Add 1 to the number of prongs/edges.
  for (int i = 1; i < (int)code.size(); i+=2) ++code[i];

  /* From here duplicate code from constructor from coding. */
  /* This is not very nice but problems otherwise (tmp vector?). */
  /* traintracl(code); */

  // Train track always starts with an uncusped monogon.
#if __cplusplus > 201103L && !defined(TRAINTRACKS_NO_SHARED_PTR)
  // C++14 has make_unique.
  mgv.push_back(mgonp(std::make_unique<multigon>(1)));
#else
  mgv.push_back(mgonp(new multigon(1)));
#endif
  Multigon(0).attach_edge();

  // Iterator for coding: skip initial uncusped monogon marker.
  intVec::const_iterator ci(code.begin() + detail::coding_block::length);

  // Recurse down and build train track.
  recursive_build(Multigon(0).Edge(0,0),ci);

  normalise();
}

// Recursively consume coding blocks and attach descendant multigons.
void traintrack::recursive_build(traintrack::edgep& ee,
				 traintrack::intVec::const_iterator& ci)
{
  detail::coding_block outb;

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
      outb = detail::coding_block(ci);

      if (debug)
	{
	  std::cerr << "Out: " << outb.nprongs << outb.prong << outb.nedges;
	  std::cerr << outb.edge << std::endl;
	}

      mgv[in_m]->attach_edge(outb.prong,outb.edge);
      // Recurse!
      recursive_build(mgv[in_m]->Edge(outb.prong,outb.edge),ci);

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

  // If this is called explicitly, then we are labeling multiprongs.
  if (label_multiprongs)
    mgv[m]->set_label(lb);
  else
    {
      std::cerr << "Error in traintrack::traintrack::set_label(): ";
      std::cerr << "flag label_multiprongs must be set.\n";
      exit(1);
    }
  normalise();
}

// Give unique label to each punctured multigon; label 0 if unpunctured.
void traintrack::pure_braid()
{
  isnormalised = false;

  // If this is called explicitly, then we are labeling multiprongs.
  if (!label_multiprongs)
    {
      std::cerr << "Error in traintrack::traintrack::pure_braid(): ";
      std::cerr << "flag label_multiprongs must be set.\n";
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
bool traintrack::fold(const int f)
{
  // Need normalised train track.
  require_normalised("traintrack::fold");

  int fcusp = f/2;
  int fdir = 1 - 2*(f % 2);

  if (f >= foldings())
    {
      std::cerr << "Nonexistent cusp in traintrack::traintrack::fold.\n";
      std::exit(1);
    }
  // Loop over cusps from first uncusped monogon, assuming a
  // normalised train track. (Ambiguous for multigons with cyclic
  // symmetry, such as 2121?)

  // Since normalised, number cusps from 0th monogon.
  int mono = 0;

  // Recurse down and count cusps.
  // Start by finding the multigon the edge is attached to, and which
  // prong.
  int pmono, pemono;
  multigon *egmono =
    Multigon(mono).Edge(0,0)->target_multigon(&Multigon(mono),pmono,pemono);

  // Find the monogon, prong, edge of the cusp with number fcusp.
#if __cplusplus > 199711L
  multigon* mmc = nullptr;
#else
  multigon* mmc = 0;
#endif
  int pc, ec, cs = fcusp;
  recursive_find_cusp(*egmono,pmono,pemono,cs,mmc,pc,ec);

#if __cplusplus > 199711L
  if (mmc == nullptr)
#else
  if (mmc == 0)
#endif
    {
      std::cerr << "Could not find cusp in traintrack::traintrack::fold.\n";
      std::exit(1);
    }

  // Fold at the cusp.
  return fold(*mmc,pc,ec,fdir);
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

  int cs = f/2;
  int mono = 0;

  // Find the 0th monogon, where coding starts.
  for (; mono < (int)mgv.size(); ++mono)
    {
      if (Multigon(mono).prongs() == 1 && Multigon(mono).edges() == 1) break;
    }
  if (mono == (int)mgv.size())
    {
      std::cerr << "Could not find monogon in traintrack::traintrack::fold_cusp_location.\n";
      std::exit(1);
    }

  // Start by finding the multigon the monogon edge is attached to.
  int pmono, pemono;
  multigon* egmono =
    Multigon(mono).Edge(0,0)->target_multigon(&Multigon(mono),pmono,pemono);

  // Reuse recursive cusp ordering to locate cusp cs.
  mmc = 0;
  pc = -1;
  ec = -1;
  recursive_find_cusp(*egmono,pmono,pemono,cs,mmc,pc,ec);

  if (mmc == 0 || pc < 0 || ec < 0)
    {
      std::cerr << "Could not resolve cusp in traintrack::traintrack::fold_cusp_location.\n";
      std::exit(1);
    }
}

// Return zero-based infinitesimal edge index selected by fold f.
int traintrack::fold_infinitesimal_index(const int f) const
{
  multigon* mmc = 0;
  int pc = -1, ec = -1;
  fold_cusp_location(f,mmc,pc,ec);

  int mi = multigon_index(mmc);
  return multigon_prong_index(mi,pc);
}

// Return signed-generator label for fold f after nmain main generators.
int traintrack::fold_infinitesimal_generator(const int f, const int nmain) const
{
  multigon* mmc = 0;
  int pc = -1, ec = -1;
  fold_cusp_location(f,mmc,pc,ec);

  const int mi = multigon_index(mmc);
  const int infix = multigon_prong_index(mi,pc);
  return nmain + infix + 1;
}

// Depth-first cusp locator used by fold() and fold_cusp_location().
bool traintrack::recursive_find_cusp(multigon& mm,
				     const int pin, const int ein,
				     int& fcusp, multigon*& mmc,
				     int& pc, int& ec) const
{
  int p = pin, e = ein;

  do
    {
      // Is there a cusp here?  There is, as long as it's not the
      // last edge on a prong.
      if (e < mm.edges(p)-1)
	{
	  if (fcusp-- == 0)
	    {
	      // We've found the cusp!  We're done here.  Send
	      // back the location of the cusp.
	      mmc = &mm;
	      pc = p;
	      ec = e;
	      return true;
	    }
	}

      // Find the next multigon down.
      int pout, eout;
      multigon *ed = mm.Edge(p,e)->target_multigon(&mm,pout,eout);

      // Don't recurse down the entry edge.
      if (!(p == pin && e == ein))
	{
	  // Recurse only if it's not an uncusped monogon down there.
	  if (ed->edges() > 1)
	    {
	      if (recursive_find_cusp(*ed,pout,eout,fcusp,mmc,pc,ec))
		return true;
	    }
	}
      mm.cycle_edges(p,e);	// Increment the edge and prong.
    }
  while (!(p == pin && e == ein));

  return false;
}

// Fold cusp c of prong p of multigon m in direction dir.
//   dir = 1 clockwise, dir = 1 anticlockwise.
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

// Assumes a track is normalised.
traintrack::dblVec traintrack::weights(const int mono) const
{
  require_normalised("traintrack::weights(get)");

  // Start the weights vector with the monogon.
  dblVec wv;
  wv.push_back(Multigon(mono).Edge(0,0)->weight());

  // Recurse down and collect weights.
  // Start by finding the multigon the edge is attached to, and which
  // prong.
  int pmono, pemono;
  multigon *egmono =
    Multigon(mono).Edge(0,0)->target_multigon(&Multigon(mono),pmono,pemono);

  recursive_get_weights(*egmono,pmono,pemono,wv);

  return wv;
}

void traintrack::recursive_get_weights(const multigon& mm,
				       const int pin, const int ein,
				       dblVec& wv) const
{
  int p = pin, e = ein;
  mm.cycle_edges(p,e);	// Increment the edge and prong.

  do
    {
      // Record the weight.
      wv.push_back(mm.Edge(p,e)->weight());

      // Find the next multigon down.
      int pout, eout;
      multigon *ed = mm.Edge(p,e)->target_multigon(&mm,pout,eout);

      // Recurse only if it's not an uncusped monogon down there.
      if (ed->edges() > 1)
	{
	  recursive_get_weights(*ed,pout,eout,wv);
	}
      mm.cycle_edges(p,e);	// Increment the edge and prong.
    }
  while (!(p == pin && e == ein));
}

// Assumes a track is normalised.
traintrack::dblVec::const_iterator
traintrack::weights(traintrack::dblVec::const_iterator wi)
{
  require_normalised("traintrack::weights(set)");

  // Since normalised, distribute weights from 0th monogon.
  int mono = 0;

  // Start the weights vector with the monogon.
  Multigon(mono).Edge(0,0)->weight(*wi++);

  // Recurse down and collect weights.
  // Start by finding the multigon the edge is attached to, and which
  // prong.
  int pmono, pemono;
  multigon *egmono =
    Multigon(mono).Edge(0,0)->target_multigon(&Multigon(mono),pmono,pemono);

  recursive_set_weights(*egmono,pmono,pemono,wi);

  return wi;
}

void traintrack::recursive_set_weights(const multigon& mm,
				       const int pin, const int ein,
				       dblVec::const_iterator& wi)
{
  int p = pin, e = ein;
  mm.cycle_edges(p,e);	// Increment the edge and prong.

  do
    {
      // Copy the weight.
      mm.Edge(p,e)->weight(*wi++);

      // Find the next multigon down.
      int pout, eout;
      multigon *ed = mm.Edge(p,e)->target_multigon(&mm,pout,eout);

      // Recurse only if it's not an uncusped monogon down there.
      if (ed->edges() > 1)
	{
	  recursive_set_weights(*ed,pout,eout,wi);
	}
      mm.cycle_edges(p,e);	// Increment the edge and prong.
    }
  while (!(p == pin && e == ein));
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
