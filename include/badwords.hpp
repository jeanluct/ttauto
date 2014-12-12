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

#ifndef BADWORDS_HPP
#define BADWORDS_HPP

#include <iostream>
#undef TTAUTO_HASH_BADWORDS
#ifdef TTAUTO_HASH_BADWORDS
#include <ext/hash_set>
#else
#include <list>
#endif
#include <jlt/vector.hpp>
#include <jlt/matrix.hpp>
#include "ttfoldgraph.hpp"
#include "traintracks_util.hpp"
#include "folding_path.hpp"

namespace traintracks {

// Find paths with matrices of equal patterns of zeros-and-nonzeros.
//
// See Ham&Song preprint (2006).

#ifdef TTAUTO_HASH_BADWORDS
template<class TrTr>
struct direct_equal_to
{
  bool operator()(const folding_path<TrTr>& p1, const folding_path<TrTr>& p2)
    const
  {
    return traintracks::direct_equal(p1,p2);
  }
};

template<class TrTr>
jlt::matrix<__gnu_cxx::hash_set<folding_path<TrTr>,
				typename folding_path<TrTr>::hash,
				direct_equal_to<TrTr> > >
#else
template<class TrTr>
jlt::matrix<std::list<folding_path<TrTr> > >
#endif
badwords(const ttfoldgraph<TrTr>& ttg, const int maxplen)
{
  using std::cout;
  using std::endl;
  using std::setw;

  typedef jlt::mathmatrix<int> Mat;
  typedef folding_path<TrTr> fpath;
#ifdef TTAUTO_HASH_BADWORDS
  typedef __gnu_cxx::hash_set<fpath,typename fpath::hash,
			    direct_equal_to<TrTr> > badword_list;
#else
  typedef std::list<folding_path<TrTr> > badword_list;
#endif
  typedef typename badword_list::iterator pit;
  typedef typename badword_list::const_iterator cpit;

  if (maxplen == 0) return jlt::matrix<badword_list>();

  cout << "Finding bad words...\n";

  const int minplen = 1;
  int total_bad_words = 0;

  jlt::matrix<badword_list> pbm(ttg.vertices(),maxplen);

  for (int v = 0; v < ttg.vertices(); ++v)
    {
      badword_list pb;

      cout << "  Vertex " << setw(4) << v << ",  length";

      for (int plen = minplen; plen <= maxplen; ++plen)
	{
	  // A path from vertex v of length plen, with all 0 foldings.
	  fpath p(ttg,v,plen);
	  // Make a different list for the new words we find.
	  badword_list pbnew;

	  cout << " " << plen; cout.flush();

	  do
	    {
	      // Check for closed paths.
	      if (p.closed())
		{
		  // Use this as the reference pattern.
		  Mat TM(p.transition_matrix());
		  // Now iterate path once and check for the same
		  // pattern.
		  /* Only looking for square: make this more general? */
		  Mat TMit(TM*TM);
		  if (pattern_equal(TM,TMit))
		    {
		      // We've found a repeated pattern.
		      //
		      // Does it match any of the previous ones? Check
		      // for this.
		      bool match_prev = false;
		      if (plen > 1 && !pb.empty())
			{
			  // Multiply up to check_len words together.
			  int check_len = plen;
			  // Make a vector of iterators.
			  jlt::vector<cpit> iv(check_len);
			  // All the iterators begin at the start of
			  // the list of words, which all start and end at
			  // the current vertex.
			  for (int l = 0; l < check_len; ++l)
			    iv[l] = pb.begin();
			  fpath pw(ttg,v);
			  bool incr;
			  do {
			    incr = false;

			    // Form words using the list iterators.
			    pw.clear();
			    for (int l = 0; l < check_len; ++l)
			      {
				// If the path is already too long,
				// don't bother.
				if ((pw.length() + iv[l]->length())
				    > p.length()) break;

				pw *= *(iv[l]);

				if (direct_equal(p,pw))
				  {
				    // We matched the path... break out!
				    match_prev = true;
				    break;
				  }
			      }
			    // We matched the path... break out again!
			    if (match_prev) break;

			    // We haven't found a match... increment.
			    for (int l = 0; l < check_len; ++l)
			      {
				// iterate through list of bad words.
				++iv[l];
				if (iv[l] != pb.end())
				  {
				    // We haven't reached the end.
				    incr = true;
				    break;
				  }
				else
				  {
				    // We've reached the end... reset
				    // iterator at that position, and
				    // let the loop move on to the
				    // next position.
				    iv[l] = pb.begin();
				  }
			      }
			  } while (incr);
			}
		      // I think the Ham&Song condition on the matrix
		      // is always satisfied.
		      if (!match_prev)
			{
			  cout << "."; cout.flush();
#ifdef TTAUTO_HASH_BADWORDS
			  pbnew.insert(p);
#else
			  pbnew.push_back(p);
#endif
			}
		    }
		}
	    }
	  while (++p); // Try all possible foldings at a fixed length.

	  // Copy the bad words we found for this length to the master list.
#ifdef TTAUTO_HASH_BADWORDS
	  pb.insert(pbnew.begin(),pbnew.end());
#else
	  std::copy(pbnew.begin(),pbnew.end(),back_inserter(pb));
#endif
	}

      total_bad_words += pb.size();
      cout << " total " << pb.size() << endl;

      // Now classify the paths according to length.
      for (pit i = pb.begin(); i != pb.end(); ++i)
	{
	  // Square all the paths.  This is because the bad words are
	  // actually the ones that are repeated twice or more, which
	  // means they can be collapsed to one instance.
#ifdef TTAUTO_HASH_BADWORDS
	  pbm(v,i->length()-1).insert((*i)*(*i));
#else
	  pbm(v,i->length()-1).push_back((*i)*(*i));
#endif
	}
    }

  cout << "Identified a grand total of " << total_bad_words;
  cout << " bad words." << endl;

  return pbm;
}

} // namespace traintracks

#endif // BADWORDS_HPP
