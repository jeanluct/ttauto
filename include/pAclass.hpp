// <LICENSE
//   ttauto: a C++ library for building train track automata
//
//   Copyright (C) 2010--2014 Jean-Luc Thiffeault and Erwan Lanneau
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

#ifndef PACLASS_HPP
#define PACLASS_HPP

#include <iostream>
#include <ext/hash_map>
#include <jlt/mathmatrix.hpp>
#include <jlt/polynomial.hpp>
#include "folding_path.hpp"
#include "traintracks_util.hpp"

namespace traintracks {

// Info about a pA class, that is pA's with the same dilatation.
template <class TrTr>
class pAclass
{
  static const int debug = 0;

public:
  typedef jlt::mathmatrix<int> Mat;
  typedef jlt::polynomial<int> Polynomial;

private:
  Polynomial charpoly;		// Characteristic polynomial of the class.
  double lambda;		// Dilatation of the class.
  int max_paths_to_save;	// Cap the max # of paths to save in
				// the class.  This is necessary since
				// there can be thousands, and often
				// we just want a representative.

  // Get hash function from folding_path<TrTr>.
  typedef typename folding_path<TrTr>::hash path_hash;
  // Container for the list of paths and transition matrices.
  //   Note that we use hash_map instead of map: it would be nice to
  //   use map since we could easily keep the list sorted by length,
  //   but when there are millions of pAs it's better to use
  //   folding_path's hash function to speed up checking.
  typedef typename __gnu_cxx::hash_map<folding_path<TrTr>,Mat,path_hash>
  pathlist;

  pathlist pathl;	// List of vertex paths/transition matrix
                        // pairs in this pseudo-Anosov class.

private:

  void copy(pAclass& pacnew, const pAclass& pacexist)
  {
    pacnew.pathl = pacexist.pathl;
    pacnew.charpoly = pacexist.charpoly;
    pacnew.lambda = pacexist.lambda;
    pacnew.max_paths_to_save = pacexist.max_paths_to_save;
  }

public:
  /* Would it be better to figure out the dilatation from the polynomial? */
  pAclass(const Polynomial& charpoly_ = Polynomial(), const double lambda_ = 1,
	  const int max_paths_to_save_ = 1)
    : charpoly(charpoly_), lambda(lambda_),
      max_paths_to_save(max_paths_to_save_) {}

  pAclass(const pAclass& pac) { copy(*this,pac); }

  pAclass& operator=(const pAclass& pac)  { copy(*this,pac); return *this; }

  const Polynomial& polynomial() const { return charpoly; }

  double dilatation() const { return lambda; }

  int number_of_paths() const { return pathl.size(); }

  // The minimum path length over all the paths in this class.
  int shortest() const
  {
    int minlen = -1;
    for (typename pathlist::const_iterator pit = pathl.begin();
	 pit != pathl.end(); ++pit)
      {
	int len = pit->first.length();
	if (minlen == -1 || len < minlen) minlen = len;
      }
    return minlen;
  }

  // The maximum path length over all the paths in this class.
  int longest() const
  {
    int maxlen = 0;
    for (typename pathlist::const_iterator pit = pathl.begin();
	 pit != pathl.end(); ++pit)
      {
	int len = pit->first.length();
	if (maxlen == 0 || len > maxlen) maxlen = len;
      }
    return maxlen;
  }

  // Add a pA path to the list.
  void add_path(const folding_path<TrTr>& p, const Mat& TM)
  {
    // Is the path closed?  It better be.
    if (!p.closed())
      {
	std::cerr << "Vertex path not closed in pAclass::add_path.\n";
	std::exit(1);
      }

    // Do we have a representative with this length yet?
    bool got_that_length = false;
    for (typename pathlist::const_iterator pit = pathl.begin();
	 pit != pathl.end(); ++pit)
      {
	if (p.length() == pit->first.length())
	  {
	    got_that_length = true;
	    break;
	  }
      }

    // Cap the number of paths kept at max_paths_to_save.  However,
    // make an exception if it's a new "unique" length: we don't want
    // to miss those.
    if ((int)pathl.size() < max_paths_to_save || max_paths_to_save == 0
	|| !got_that_length)
      {
	// Insert new path.  Hashed container ensures we don't duplicate.
	pathl[p] = TM;
      }
  }

  std::ostream& print_paths(std::ostream& strm = std::cout,
			    const int max_paths_to_print = 0) const
  {
    int printed = 0;
    for (typename pathlist::const_iterator pit = pathl.begin();
	 pit != pathl.end(); ++pit, ++printed)
      {
	if (printed < max_paths_to_print || max_paths_to_print == 0)
	  {
	    strm << pit->first << " [" << pit->first.vertices() << "]";
	    if (printed != max_paths_to_print-1 &&
		printed != number_of_paths()-1) strm << ", ";
	  }
	else
	  { strm << " ..."; break; }
      }
    return strm;
  }

  std::ostream& print_pA_MathematicaForm(std::ostream& strm = std::cout,
					 const int max_paths_to_print = 0)
    const
  {
    using std::setw;

    int prec = strm.precision();
    strm.precision(12);
    strm.setf(std::ios::showpoint);
    strm << "{polynomial -> " << charpoly;
    strm << ", lambda -> " << lambda << ", pathlist -> {";
    strm.unsetf(std::ios::showpoint);
    strm.precision(prec);
    int printed = 0;
    typename pathlist::const_iterator pit = pathl.begin();
    for (; pit != pathl.end(); ++pit, ++printed)
      {
	if (printed < max_paths_to_print || max_paths_to_print == 0)
	  {
	    strm << "{";
	    pit->first.vertices().printMathematicaForm(strm) << ",";
	    pit->second.printMathematicaForm(strm) << "}";
	    if (printed != max_paths_to_print-1 &&
		printed != number_of_paths()-1) strm << ",";
	  }
      }
    strm << "}}";
    return strm;
  }
};

} // namespace traintracks

#endif // PACLASS_HPP
