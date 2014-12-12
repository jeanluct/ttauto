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

#include <iostream>
#include <jlt/vector.hpp>
#include "edge.hpp"
#include "multigon.hpp"
#include "traintrack.hpp"

// Build predefined train tracks.

// These are constructors, methods, and helper functions for the class
// traintrack.

namespace traintracks {

//
// Constructors
//
//   These just call the appropriate helper function.
//

// Track with N monogons around punctures and a (N-2)-gon on the
// boundary (all cusps).
traintrack::traintrack(const int N)
{
  ttbuild_all_monogons(N);
}

// Track with N monogons around punctures, an unpunctured K-gon,
// and an (N-K)-gon on the boundary.
traintrack::traintrack(const int N, const int K)
{
  if (K == 1)
    {
      std::cerr << "Specify monogons via N, the number of punctures,";
      std::cerr << " in traintrack::traintrack(N,K).\n";
      std::exit(1);
    }

  if (N < K+1)
    {
      std::cerr << "Not a valid train track for K > N-1";
      std::cerr << " in traintrack::traintrack(N,K).\n";
      std::exit(1);
    }

  intVec Kv(1);
  Kv[0] = K;

  ttbuild_monogoncusps(N,Kv);
}

// Track with N monogons around punctures, an unpunctured K1-gon, an
// unpunctured K2-gon, and an (N-(K1+K2)+2)-gon on the boundary.
traintrack::traintrack(const int N, const int K1, const int K2)
{
  if (K1 == 1 || K2 == 1)
    {
      std::cerr << "Specify monogons via N, the number of punctures,";
      std::cerr << " in traintrack::traintrack(N,K1,K2).\n";
      std::exit(1);
    }

  if (N < K1+K2-1)
    {
      std::cerr << "Not a valid train track for N < K1+K2-1";
      std::cerr << " in traintrack::traintrack(N,K1,K2).\n";
      std::exit(1);
    }

  intVec Kv(2);
  Kv[0] = K1;
  Kv[1] = K2;

  traintrack(N,Kv);
}

// Track with N monogons around punctures, unpunctured K1-, K2-,
// and K3-gons, and an (N-(K1+K2+K3)+4)-gon on the boundary.
traintrack::traintrack(const int N, const int K1, const int K2, const int K3)
{
  if (K1 == 1 || K2 == 1 || K3 == 1)
    {
      std::cerr << "Specify monogons via N, the number of punctures,";
      std::cerr << " in traintrack::traintrack(N,K1,K2,K3).\n";
      std::exit(1);
    }

  if (N < K1+K2+K3-2)
    {
      std::cerr << "Not a valid train track for N < K1+K2+K3-2";
      std::cerr << " in traintrack::traintrack(N,K1,K2,K3).\n";
      std::exit(1);
    }

  intVec Kv(3);
  Kv[0] = K1;
  Kv[1] = K2;
  Kv[2] = K3;

  traintrack(N,Kv);
}

// Track with N monogons around punctures and unpunctured K-gons
// specified by a vector Kv.
traintrack::traintrack(const int N, const intVec& Kv)
{
  ttbuild_monogons_and_multigons(N,Kv);
}

#if 0
// Track with all multigons specified by Kv, including monogons.
/* Cannot do this since it clashes with creating a track from coding. */
/* Introduce coding class. */
traintrack::traintrack(const intVec& Kv)
{
  int N = 0;
  intVec Kv2;

  for (int m = 0; m < (int)Kv.size(); ++m)
    {
      if (Kv[m] == 1)
	++N;
      else
	Kv2.push_back(Kv[m]);
    }

  ttbuild_monogons_and_multigons(N,Kv2);
}
#endif

//
// Helper functions
//

// Track with N monogons around punctures and a (N-2)-gon on the
// boundary (all cusps).
void traintrack::ttbuild_all_monogons(const int N)
{
  mgv.push_back(mgonp(new multigon(1)));	// A monogon (puncture)

  // All the edges are attached to monogon 0.
  for (int i = 0; i < N-1; ++i)
    Multigon(0).attach_edge(0,i);

  // Cap-off by attaching monogons to the free edge endings.
  capoff();

  check();
  normalise();
}

// Track with N monogons around punctures and unpunctured K-gons
// specified by a vector Kv.
void traintrack::ttbuild_monogons_and_multigons(const int N,
						const intVec& Kv)
{
  if (Kv.empty())
    {
      ttbuild_all_monogons(N);
      return;
    }

  for (int i = 0; i < (int)Kv.size(); ++i)
    {
      if (Kv[i] == 1)
	{
	  std::cerr << "Specify monogons via N, the number of punctures,";
	  std::cerr << " in traintrack::traintrack(N,Kv)\n";
	  std::exit(1);
	}
    }

  // Start with the first two multigons.
  int N0 = boundary_prongs(N,Kv) + (Kv[0]-1) + (Kv[1]-1);
  intVec Kv0(2);
  Kv0[0] = Kv[0]; Kv0[1] = Kv[1];
  ttbuild_monogoncusps(N0,Kv0);

  // Now successively replace a monogon by a multigon in the list.
  for (intVec::const_iterator i = Kv.begin()+2; i != Kv.end(); ++i)
    {
      // If we normalise, uncusped monogons are always first.
      int monog = 0;
      monogon_to_multigon(monog,*i);
      capoff();
      normalise();
    }

  check();
}

// Track with N monogons around punctures and unpunctured K-gons
// specified by a vector K.
//
// All the multigons are anchored at a monogon, so not all tracks
// can be made this way.  See traintrack::traintrack(N,K).
void traintrack::ttbuild_monogoncusps(const int N, const intVec& Kv)
{
  // See diagram in Marseille ring NB 2 for notation.

  // Tracks with one cusp cannot have all cusps on a monogon unless
  // they only have two (non-monogon) multigons.
  if (boundary_prongs(N,Kv) == 1 && Kv.size() > 2)
    {
      std::cerr << "Use traintrack(N,Kv) for 1 cusp and >2 multigons";
      std::cerr << " in traintrack::ttbuild_monogoncusps.\n";
      std::exit(1);
    }

  const int L = Kv.size();

  for (int i = 0; i < L; ++i)
    {
      if (Kv[i] == 1)
	{
	  std::cerr << "Specify monogons via N, the number of punctures,";
	  std::cerr << " in traintrack::ttbuild_monogoncusps.\n";
	  std::exit(1);
	}
    }

  int Nmin = 1;
  for (int i = 0; i < L; ++i) Nmin += Kv[i]-1;
  if (N < Nmin)
    {
      std::cerr << "Need at least " << Nmin << " monogons";
      std::cerr << " in traintrack::ttbuild_monogoncusps.\n";
      std::exit(1);
    }

  intVec P(L), Q(L);

  Q[L-1] = N-1;
  P[L-1] = Q[L-1] - Kv[L-1] + 2;
  for (int i = L-2; i >= 0; --i)
    {
      Q[i] = P[i+1] - 1;
      P[i] = Q[i] - Kv[i] + 2;
    }

  traintrack tt;
  mgv.push_back(mgonp(new multigon(1)));		// monogon
  for (int i = 0; i < L; ++i)
    {
      mgv.push_back(mgonp(new multigon(Kv[i])));	// Kv[i]-gon
    }

  // Attach edges to monogon 0.
  for (int i = 0; i <= (P[0]-1)+(L-1); ++i)
    Multigon(0).attach_edge(0,i);

  for (int i = 0; i < L; ++i)
    {
      // Edges on prong 0, edges [P0-1 + 0,P0-1 + (L-1)] of monogon 0
      // are attached to the 0th prong of multigons [1,L].
      multigon::edgep ep = Multigon(0).Edge(0,P[0]-1+i);
      Multigon(i+1).attach_edge(ep);
    }

  // Cap-off by attaching monogons to the free edge endings.
  capoff();

  check();
  normalise();
}

// Add a monogon to all the free edge endings, creating new edges if
// some prong edges are still free.
void traintrack::capoff()
{
  int mgin = multigons();

  for (int m = 0; m < mgin; ++m)
    {
      for (int p = 0; p < Multigon(m).prongs(); ++p)
	{
	  for (int e = 0; e < Multigon(m).edges(p); ++e)
	    {
	      if (Multigon(m).prong_edge_is_unattached(p,e))
		{
		  // There's no edge there at all, so add one.
		  // It will get capped off in a bit...
		  Multigon(m).attach_edge(p,e);
		}
	      bool ua0 = Multigon(m).Edge(p,e)->is_unattached(0);
	      bool ua1 = Multigon(m).Edge(p,e)->is_unattached(1);
	      if (ua0 && ua1)
		{
		  std::cerr << "Both ends unattached";
		  std::cerr << " in traintrack::capoff.\n";
		  std::exit(1);
		}
	      if (ua0 || ua1)
		{
		  mgv.push_back(mgonp(new multigon(1)));
		  multigon::edgep ep = Multigon(m).Edge(p,e);
		  mgv.back()->attach_edge(ep);
		}
	    }
	}
    }
}

void traintrack::monogon_to_multigon(const int L, const int K)
{
  typedef traintrack::mgonp	mgonp;
  typedef multigon::edgep	edgep;

  if (Multigon(L).prongs() != 1)
    {
      std::cerr << "Not a monogon";
      std::cerr << " in traintracks::ttbuild_monogon_to_multigon\n";
      std::exit(1);
    }

  // Save a pointer to the edge multigon L was pointing to.
  edgep ep = Multigon(L).Edge(0,0);
  // Detach the edge.
#ifdef TTAUTO_NO_BOOST
  Multigon(L).Edge(0,0)->detach_from_multigon(mgv[L]);
#else
  Multigon(L).Edge(0,0)->detach_from_multigon(mgv[L].get());
#endif
  // Now make a new multigon, overwriting the old one.
  mgv[L] = mgonp(new multigon(K));

  // Attach it to the edge we freed up.
  Multigon(L).attach_edge(ep,0,0);
}


//
// Non-member helper functions for building tracks
//

jlt::vector<traintrack> ttbuild_list2(const int N)
{
  typedef traintrack::intVec intVec;
  typedef std::list<intVec>::const_iterator cit;

  jlt::vector<traintrack> ttv;

  // Tracks with punctured monogons and bigons only.
  for (int N2 = 0; N2 <= N-3; ++N2)
    {
      jlt::vector<traintrack> ttv2 = ttbuild_list(N-N2,N2);
      ttv.insert(ttv.end(),ttv2.begin(),ttv2.end());
    }

  return ttv;
}

jlt::vector<traintrack> ttbuild_list(const int N, const int N2)
{
  typedef traintrack::intVec intVec;
  typedef std::list<intVec>::const_iterator cit;

  if (N < 3)
    {
      std::cerr << "traintracks::ttbuild_list: Need at least 3 puntures.\n";
      std::exit(1);
    }

  std::list<intVec> sdata_list;

  // The maximum prong order is N-1.
  intVec sdata(N-1);
  sdata[0] = N;    // 1-prongs
  sdata[1] = N2;   // 2-prongs
  sdata[2] = -1;

  // The maximum possible number of singularities of a given
  // prongness.  Obtained from Euler-Poincare formula by setting
  // Nsep=1 (number of boundary separatrices).
  intVec sdatamax(N-1);
  sdatamax[0] = N;    // 1-prongs
  sdatamax[1] = N2;   // 2-prongs
  for (int p = 2; p < N-1; ++p) sdatamax[p] = (N-3)/((p+1)-2);

  bool incr = true;
  if (N == 3) { incr = false; sdata[2] = 0; sdata_list.push_back(sdata); }
  while (incr)
    {
      for (int p = 2; p < N-1; ++p)
	{
	  ++sdata[p];
	  if (sdata[p] > sdatamax[p])
	    {
	      if (p == N-2) incr = false;
	      sdata[p] = 0;
	    }
	  else
	    {
	      // Verify that Euler-Poincare can be satisfied with at
	      // least a monogon on the boundary.
	      if (boundary_prongs(sdata) >= 1) sdata_list.push_back(sdata);
	      break;
	    }
	}
    }

  // Convert to individual multigon format.
  std::list<intVec> Klist;
  for (cit i = sdata_list.begin(); i != sdata_list.end(); ++i)
    {
      intVec K;
      for (int p = (int)i->size()-1; p > 0; --p)
	{
	  int nsing = (*i)[p];
	  if (nsing != 0) K.insert(K.end(),nsing,p+1);
	}
      Klist.push_back(K);
    }

  // Sort to agree with convention (not so great: increasing order of
  // rightmost digit!).
  Klist.sort(compare_strata);
#if 0
  copy(Klist.begin(),Klist.end(),
       std::ostream_iterator<intVec>(std::cout,"\n"));
#endif

  // Build list of train tracks.
  jlt::vector<traintrack> ttv;
  for (cit i = Klist.begin(); i != Klist.end(); ++i)
      if (i->empty())
	  ttv.push_back(traintrack(N));
      else if (i->size() == 1)
	  ttv.push_back(traintrack(N,(*i)[0]));
      else
	  ttv.push_back(traintrack(N,*i));

  if (traintrack::debug)
    {
      std::cout << "Built list of " << ttv.size() << " train tracks";
      std::cout << " for " << N << " monogons";

      if (N2 > 0)
	{
	  std::cout << " and " << N2 << " bigon";
	  std::cout << (N2 > 1 ? "s." : ".");
	}
      std::cout << std::endl;
    }

  return ttv;
}

// Comparison function for sorting strata.
bool compare_strata(const traintrack::intVec& first,
		    const traintrack::intVec& second)
{
  if (second.empty()) return false;
  if (first.empty()) return true;

  if (first.size() < second.size()) return true;
  if (first.size() > second.size()) return false;

  // Sort the two vectors.
  traintrack::intVec sfirst = first;
  traintrack::intVec ssecond = second;
  std::sort(sfirst.begin(),sfirst.end());
  std::sort(ssecond.begin(),ssecond.end());

  // Use lexicographical compare for vectors of equal length on the
  // sorted vectors.
  if (sfirst < ssecond) return true;

  return false;
}

// How many prongs on the boundary for singularity data sdata?
int boundary_prongs(traintrack::intVec sdata)
{
  const int EC_disk = 1;
  int id = 0;

  for (int p = 0; p < (int)sdata.size(); ++p) id += sdata[p]*((p+1)-2);

  return -(2*EC_disk+id);
}

// How many prongs on the boundary for N punctures and multigons given
// by the vector Kv?
int boundary_prongs(const int N, const traintrack::intVec Kv)
{
  const int EC_sphere = 2;

  int EP = -N;

  for (int i = 0; i < (int)Kv.size(); ++i) EP += Kv[i]-2;

  return (-2*EC_sphere - EP) + 2;
}

} // namespace traintracks
