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

#ifndef TTAUTO_TTAUTO_HPP
#define TTAUTO_TTAUTO_HPP

#include <iostream>
#include <fstream>
#include <list>
#include <map>
#include <stack>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <jlt/matrix.hpp>
#include <jlt/mathmatrix.hpp>
#include <jlt/vector.hpp>
#ifdef TTAUTO_USE_FORTRAN
#include <jlt/eigensystem.hpp>
#endif
#include <jlt/polynomial.hpp>
#include "ttauto/folding_path.hpp"
#include "ttauto/ttfoldgraph.hpp"
#include "ttauto/pAclass.hpp"
#include "ttauto/badwords.hpp"

namespace ttauto {

// Find the largest real root of a polynomial via Newton iteration.
//
// The polynomial comes from a Perron-Frobenius matrix, so the largest
// real root is guaranteed to be real and unique.
inline double findroot(const jlt::polynomial<int>& p,
		       const double x0, const double tol)
{
  double px(p(x0)), x(x0);

  int i = 0;
  const int itmax = 100;

  while (std::abs(px) > tol && i++ < itmax)
    {
      x = x - px / p.derivative_at(x);
      px = p(x);
    }

  if (i == itmax)
    throw
      jlt::failed_to_converge<double>
      ("Failed to converge to specified accuracy.\n",std::abs(px));

  return x;
}


template<class TrTr>
class ttauto
{
  struct ltPoly;
public:
  //
  // Typedefs
  //

  typedef jlt::mathmatrix<int>			Mat;
  typedef jlt::vector<int>			Vec;
  typedef typename pAclass<TrTr>::Polynomial	Poly;

  // The container for the pseudo-Anosovs and their paths in the automaton.
  // They are kept sorted according to the above function ltPoly.
  typedef std::map<Poly,pAclass<TrTr>,ltPoly>	pAlist;

  typedef folding_path<TrTr>			fpath;
  typedef long long int				llint;
#ifdef TTAUTO_HASH_BADWORDS
  typedef std::unordered_set<fpath,typename fpath::hash,
			     direct_equal_to<TrTr> > badword_list;
#else
  typedef std::list<fpath>			badword_list;
#endif

private:
  static constexpr int debug = 0;
  static const bool exploit_symmetries = TrTr::exploit_symmetries;

  ttfoldgraph<TrTr> ttg;	// Train track graph.
  const int nfoldsmax;		// Maximum # of foldings at each vertex.
  const int n;			// Number of edges in train track.

  // Store at each vertex as we descend the graph.
  fpath p;			// Folding path to current vertex.
  std::stack<int,Vec> tryv;	// Next folding to try.
  std::stack<Mat> TMl;		// Current transition matrix.

  // Information about pA's.
  Poly charpoly;		// Characteristic polynomial of current pA path.
  double lambda;		// Dilatation of the current pA path.
  int tt0;			// Current initial vertex.
  Mat TM;			// Transition matrix of the current pA path.
  pAlist pAl;			// List of pAs.
  pAlist rejl;			// Closed irreducible paths rejected by the gate test.
  static constexpr double tol = 1e-5;// Tolerance for a match of the dilatation.

  // Criteria to reject paths.
  bool do_check_norms;		// Limit path lengths from matrix norms.
  bool do_check_gates;		// Reject candidates failing the gate test.
  double lambdamin;		// Minimum dilatation to keep (0: keep all).
  double lambdamax;		// Maximum dilatation to keep (0: keep all).
  double maxnorm;		// Maximum norm of matrix before giving up.
  int max_path_length;		// Maximum path length before giving up.
  int max_badword_length;	// How long badwords can be.
  jlt::matrix<badword_list> pbad;// Bad words list.
  Vec done;			// Vertices we have previously checked.
  Vec todo_list;		// Vertices to check.

  // Output control.
  int max_paths_save;		// How many path of each dilatation to keep.
  int max_paths_print;		// How many total paths to print.
  int print_every;		// How often to print the current path.
  const char *pAfile;		// Filename for pA data file.

  // Counters.
  llint totalpathstried;	// How many total paths tried?
  double totalpathlength;	// Total pathlength traversed?
  llint closed;			// Total closed paths encountered?
  llint primitive;		// Total closed paths with primitive matrix?
  llint pseudoAnosov;		// Total pA paths encountered?
  llint maxpathlengthexceeded;	// Total times exceeded max_path_length?
  llint normexceeded;		// Total times exceeded matrix norm?
  llint colsumexceeded;		// Total times exceeded matrix column sum?
  llint rowsumexceeded;		// Total times exceeded matrix row sum?
#ifdef TTAUTO_CHECK_SYMMETRIC_NORM
  llint symmetricnormexceeded;	// Total times exceeded matrix symmetric norm?
#endif
  llint badwordsomitted;	// Total times we omitted bad words?
  llint gatecandidates;		// Candidates reaching the gate test (closed,
				// primitive, inside the dilatation window),
				// cumulative over the whole search (the
				// other counters are per initial vertex).
  llint gaterejected;		// Of those, rejected by the gate test.

public:
  // Contructor: given an initial configuration in a train track
  // folding graph, follow all paths and record pseudo-Anosovs up to a
  // certain dilatation dilmax.
  ttauto(const ttfoldgraph<TrTr>& ttg_)
    :
      ttg(ttg_),
      nfoldsmax(ttg_.foldings()),
      n(ttg_.edges()),
      p(ttg),
      TM(n,n),
      do_check_norms(false),
      do_check_gates(true),
      lambdamin(0),
      lambdamax(0),
      maxnorm(0),
      max_path_length(15),
      max_badword_length(0),
      max_paths_save(1),
      max_paths_print(3),
      print_every(2000000),
      pAfile(0),
      gatecandidates(0),
      gaterejected(0)
  {
    eliminate_pairs();
    if (debug)
      {
	std::cout << "Checking " << todo_list.size();
	std::cout << " tracks out of " << ttg.vertices() << std::endl;
      }
  }

  //
  // Methods for setting parameters
  //

  /*

     There are basically two modes of operation: search up to a fixed
     length without checking norms, or search until the norms are
     exceeded.  The default is to search by length with a max
     pathlength of 15.  (Not a good idea to default to the other mode,
     check_norms, because when the object is first created it is
     unusable until a max dilatation is set.)

     Defaults:

       do_check_norms = false
       max_path_length = 15
       maxnorm = 0 (not used)
       max_badword_length = 0 (no bad-word pruning)

     After calling check_norms(true):

       do_check_norms = true
       max_path_length = <set by find_maxnorm>
       maxnorm = <set by find_maxnorm>
       max_badword_length = unchanged

     Calling check_norms(false) resets do_check_norms, but leaves
     everything else untouched (though only max_path_length will get used).

     Bad-word pruning is independent of both modes: badword_length()
     alone turns it on, in either.  It is off by default because it is
     only valid when hunting a minimum -- see badword_length().

  */

  // Switch to the norm-bounded mode.  This recomputes max_path_length
  // from the dilatation window (see find_maxnorm), discarding any value
  // set by max_pathlength(); so does a later max_dilatation().  To keep
  // an explicit length bound, call max_pathlength() after both.  Turning
  // this back off does not restore a length bound it replaced.
  ttauto<TrTr>& check_norms(const bool do_check_norms_ = true)
  {
    do_check_norms = do_check_norms_;
    if (do_check_norms) find_maxnorm();
    return *this;
  }

  // Bestvina-Handel gate test (issue #2): an irreducible matrix with
  // dilatation > 1 is not enough for a closed path to be pseudo-Anosov;
  // the gates at every vertex must also be connected by infinitesimal
  // edges (see traintracks/gates.hpp and devel/iss002/issue2_gates.tex).
  // On by default.  Rejected candidates are counted and kept in
  // rejected_pA_list().
  ttauto<TrTr>& check_gates(const bool do_check_gates_ = true)
  {
    do_check_gates = do_check_gates_;
    return *this;
  }

  ttauto<TrTr>& min_dilatation(const double lambdamin_)
  {
    lambdamin = lambdamin_;
    return *this;
  }

  ttauto<TrTr>& max_dilatation(const double lambdamax_)
  {
    lambdamax = lambdamax_;
    // If lambdamax is set to 0, then it makes no sense to check
    // norms.  Otherwise the next call to find_maxnorm would give a
    // ridiculous bound on path lengths.
    if (lambdamax == 0) do_check_norms = false;
    if (do_check_norms) find_maxnorm();
    return *this;
  }

  // Set to 0 for no length bound, which is only usable together with
  // check_norms(): the automaton has cycles, so the norm tests are what
  // stop the search.  search() refuses the unbounded combination rather
  // than hang.  Note check_norms() and, once it is on, max_dilatation()
  // both overwrite this; call it after them.
  ttauto<TrTr>& max_pathlength(const int max_path_length_)
  {
    max_path_length = max_path_length_;
    return *this;
  }

  // Prune folding paths that traverse a closed subpath twice in a row,
  // when that subpath's transition matrix has the same pattern of zeros
  // as its square.  Set to 0 (the default) to disable.
  //
  // This is only valid when looking for a MINIMUM dilatation.  Repeating
  // such a loop can only raise the dilatation, so the minimiser never
  // contains one; but the repeated path is usually a perfectly good
  // pseudo-Anosov class of its own, just not a minimal one, and pruning
  // silently drops it.  Do not use this when enumerating every class
  // below a bound.  See doc/ttauto.tex and Ham & Song, Exp. Math. 16
  // (2007).
  ttauto<TrTr>& badword_length(const int max_badword_length_)
  {
    max_badword_length = max_badword_length_;
    if (max_badword_length > 0)
      {
	pbad = badwords(ttg,max_badword_length);
      }
    else
      {
	pbad = jlt::matrix<badword_list>(0,0);
      }
    return *this;
  }

  ttauto<TrTr>& skip_vertices(const Vec& done_)
  {
    done = done_;
    return *this;
  }

  // Set to 0 to save all paths.
  ttauto<TrTr>& max_paths_to_save(const int max_paths_save_)
  {
    max_paths_save = max_paths_save_;
    return *this;
  }

  // Set to 0 to print all paths found and saved.
  ttauto<TrTr>& max_paths_to_print(const int max_paths_print_)
  {
    max_paths_print = max_paths_print_;
    return *this;
  }

  ttauto<TrTr>& print_path_every(const int print_every_)
  {
    print_every = print_every_;
    return *this;
  }

  ttauto<TrTr>& output_file(const char* pAfile_)
  {
    pAfile = pAfile_;
    return *this;
  }


  //
  // Output of pAs
  //

  // Access the map of recorded pseudo-Anosov classes.
  const pAlist& pA_list() const { return pAl; }

  // Candidates that passed the matrix test but failed the gate test.
  const pAlist& rejected_pA_list() const { return rejl; }

  // How many branches the bad-word prune cut off in the last search.
  llint badwords_omitted() const { return badwordsomitted; }

  // Gate-test statistics, cumulative over the whole search: candidates
  // that reached the test (closed path, primitive matrix, dilatation in
  // the window), those rejected, and the rejection rate.
  llint gate_candidates() const { return gatecandidates; }
  llint gate_rejected() const { return gaterejected; }
  double gate_rejection_rate() const
  {
    return (gatecandidates > 0 ? (double)gaterejected/gatecandidates : 0.0);
  }

  std::ostream& print_pA_list(std::ostream& strm = std::cout) const;

  std::ostream& print_rejected_pA_list(std::ostream& strm = std::cout) const;

  std::ostream&
  print_pA_list_MathematicaForm(std::ostream& strm = std::cout) const;

  //
  // Do it: Search the automaton starting from vertex tt00
  //

  // Explore folding paths from starting vertex tt00 and record pA classes.
  void search(const int tt00 = 0);

private:
  void find_maxnorm()
  {
    // A bound on the norm (total entry sum) of a transition matrix whose
    // spectral radius is at most lambdamax.  See Ham & Song, Exp. Math.
    // 16 (2007), page 170.  Note n is the number of edges, so this grows
    // very fast as the dilatation window widens.
    maxnorm = std::floor(std::pow(lambdamax,(double)n)) + n - 1;
    // Doesn't make much sense to manually limit the path length if
    // check_norms is on, so go as far as conceivable.  User will have
    // to override pathlength explicitly after calling check_norms()
    // if for some reason that's needed.
    //
    // This reuses a bound on the matrix norm as a bound on the path
    // length.  They are different quantities, and the cap is only meant
    // as a stop: the pruning in this mode is the row and column sum test
    // in check_all_norms(), which grows monotonically along a path.
    // Clamp, since lambdamax^n overflows an int for a wide window.
    const double intmax = std::numeric_limits<int>::max();
    max_path_length = (int)(maxnorm < intmax ? maxnorm : intmax);
  }

  // Build todo_list while removing one side of reflection-symmetric pairs.
  void eliminate_pairs();

  void reset_counters()
  {
    // Reset counters.
    totalpathstried = 0;
    totalpathlength = 0;
    closed = 0;
    primitive = 0;
    pseudoAnosov = 0;
    maxpathlengthexceeded = 0;
    normexceeded = 0;
    colsumexceeded = 0;
    rowsumexceeded = 0;
#ifdef TTAUTO_CHECK_SYMMETRIC_NORM
    symmetricnormexceeded = 0;
#endif
    badwordsomitted = 0;
  }

  // Run depth-first search from current tt0 and collect qualifying classes.
  bool find_pAs();

  // Advance one DFS step, including pruning and closed-path checks.
  bool descend_graph();

  // Apply matrix-based lower-bound pruning for current path.
  bool check_all_norms();

  void check_OstrovskiSchneider();

  bool backtrack(const int stps = 1)
  {
    for (int i = 0; i < stps; ++i)
      {
	p.pop_back();
	if (debug)
	  std::cerr << "Backtracking...  pathlength=" << p.length() << "\n";
	tryv.pop();
	if (lambdamax != 0 && do_check_norms) TMl.pop();
	if (p.zerolength()) return false;
      }
    return true;
  }

  // Extend path by one fold and update traversal/matrix stacks.
  void new_vertex();

  // Save info about a pA class to a list.
  // For a fixed polynomial, save paths and transition matrices.
  void record_pA();

  // Add the current path and matrix to a class list keyed by charpoly.
  void add_current_path(pAlist& l);

}; // class ttauto


//
// Public methods definitions
//

template<class TrTr>
std::ostream& ttauto<TrTr>::print_pA_list(std::ostream& strm) const
{
  using std::setw;

  int prec = strm.precision();
  strm.precision(6);
  strm.setf(std::ios::showpoint);
  int npAs = 0;

  strm << "     dilatation #  paths\n";

  for (auto it = pAl.begin(); it != pAl.end(); ++it, ++npAs)
    {
      strm << setw(3) << npAs << "  ";
      strm << it->second.dilatation() << "  ";
      strm << setw(3) << it->second.number_of_paths() << "  ";
      it->second.print_paths(strm,max_paths_print);
      strm << std::endl;
    }
  strm.unsetf(std::ios::showpoint);
  strm.precision(prec);
  return strm;
}

template<class TrTr>
std::ostream& ttauto<TrTr>::print_rejected_pA_list(std::ostream& strm) const
{
  using std::setw;

  int prec = strm.precision();
  strm.precision(6);
  strm.setf(std::ios::showpoint);
  int npAs = 0;

  strm << "     dilatation #  paths\n";

  for (auto it = rejl.begin(); it != rejl.end(); ++it, ++npAs)
    {
      strm << setw(3) << npAs << "  ";
      strm << it->second.dilatation() << "  ";
      strm << setw(3) << it->second.number_of_paths() << "  ";
      it->second.print_paths(strm,max_paths_print);
      strm << std::endl;
    }
  strm.unsetf(std::ios::showpoint);
  strm.precision(prec);
  return strm;
}


template<class TrTr>
std::ostream&
ttauto<TrTr>::print_pA_list_MathematicaForm(std::ostream& strm) const
{
  if (pAl.empty())
    {
      strm << "{}\n";
      return strm;
    }

  strm << "{";
  auto it = pAl.begin();
  do
    {
      it->second.print_pA_MathematicaForm(strm);
      if (++it != pAl.end()) strm << ",\n";
    }
  while (it != pAl.end());
  strm << "}\n";

  return strm;
}


template<class TrTr>
void ttauto<TrTr>::search(const int tt00)
{
  using std::cout;
  using std::endl;

  // An unbounded path length only terminates because the norm tests stop
  // it: the automaton has cycles, so with neither bound the search runs
  // forever.  Catch that here rather than hang.
  if (max_path_length == 0 && !do_check_norms)
    {
      std::cerr << "Error in ttauto::search(): ";
      std::cerr << "max_pathlength(0) means no length bound, which only ";
      std::cerr << "terminates with check_norms() on.\n";
      std::exit(1);
    }

  gatecandidates = 0;
  gaterejected = 0;
  rejl.clear();

  // Loop over selected vertices as initial vertex to search for pAs,
  // starting from tt00.
  for (int i = 0; i < ttg.vertices(); ++i)
    {
      tt0 = (tt00 + i) % ttg.vertices();
      cout << "\n\nVERTEX " << tt0;

      auto it =	std::find(done.begin(),done.end(),tt0);
      if (it != done.end())
	{
	  cout << "  Already done ...skipping";
	  continue;
	}

      if (exploit_symmetries)
	{
	  it = std::find(done.begin(),done.end(),ttg.symmetric_double(tt0));
	  if (it != done.end())
	    {
	      cout << "  Already done its symmetric double " << std::setw(4);
	      cout << ttg.symmetric_double(tt0) << " ...skipping";
	      continue;
	    }

	  it = std::find(todo_list.begin(),todo_list.end(),tt0);
	  if (it == todo_list.end())
	    {
	      cout << "  Reflection-symmetric to " << std::setw(4);
	      cout << ttg.symmetric_double(tt0) << " ...skipping";
	      continue;
	    }
	}

      cout << endl << endl;
      find_pAs();

      // Display some intermediate results.
      if (!pA_list().empty())
	{
	  cout << "pA candidates found";
	  if (lambdamax != 0)
	    {
	      cout << " with " << (lambdamin == 0 ? 1 : lambdamin)
		   << " < dilatation < " << lambdamax;
	      if (max_path_length != 0) cout << " and ";
	    }
	  else
	    {
	      if (max_path_length != 0) cout << " with ";
	    }
	  if (max_path_length != 0)
	    cout << "path length <= " << max_path_length;
	  cout << ":\n";
	  print_pA_list();

	  // Save all the cumulative pA path info to a file,
	  // overwriting each time.
	  if (pAfile)
	    {
	      std::ofstream pAstrm(pAfile);
	      print_pA_list_MathematicaForm(pAstrm);
	    }
	}
      else
	{
	  cout << "NO pseudo-Anosov candidates found";
	  if (lambdamax != 0)
	    {
	      cout << " with " << (lambdamin == 0 ? 1 : lambdamin)
		   << " < dilatation < " << lambdamax;
	      if (max_path_length != 0) cout << " and ";
	    }
	  else
	    {
	      if (max_path_length != 0) cout << " with ";
	    }
	  if (max_path_length != 0)
	    cout << "path length <= " << max_path_length;
	  cout << ".\n";
	}
    }

  cout << endl;

  // Bad-word pruning is only valid when hunting a minimum (see
  // badword_length).  Say so when it was on.  With a dilatation window we
  // can guess at the intent: if the lowest dilatation found is nowhere
  // near the window, the user was probably collecting every class below
  // it, and pruning may have cost them some.  Without a window there is
  // nothing to compare against, so just state the caveat.
  if (max_badword_length > 0)
    {
      if (lambdamax != 0)
	{
	  if (!pA_list().empty() &&
	      (lambdamax - pA_list().begin()->second.dilatation()) > tol)
	    {
	      cout << "\nWarning: skipping badwords might miss some pAs";
	      cout << " if not a strict minimum.\n";
	    }
	}
      else
	{
	  cout << "\nNote: badwords were skipped, so pAs that are not";
	  cout << " minimal may be missing.\n";
	}
    }
}


//
// Private methods definitions
//

//
// Sorting comparison function for the polynomials (not generally reciprocal).
//
//   Ideally, we would use the dilatation; however, we don't know that
//   easily from the polynomial.  More importantly, some polynomials
//   have the same dilatation.  So instead sort according to
//   coefficients.
//
template<class TrTr>
struct ttauto<TrTr>::ltPoly
{
  bool operator()(const ttauto<TrTr>::Poly& p1,
		  const ttauto<TrTr>::Poly& p2) const
  {
    // Smaller degree means less than.
    if (p1.degree() < p2.degree()) return true;
    if (p1.degree() > p2.degree()) return false;

    const int m = p1.degree();
    for (int i = m-1; i >= 0; --i)
      {
	// Make the leading coefficient have negative sign.
	typename ttauto<TrTr>::Poly::coeff_type
	  p1m = -p1[m]*p1[i], p2m = -p2[m]*p2[i];
	if (p1m != p2m) return (p1m < p2m);
      }
    return false;
  }
};


template<class TrTr>
void ttauto<TrTr>::eliminate_pairs()
{
  // Eliminate 1/2 of each "symmetric pair".
  if (exploit_symmetries)
    {
      for (int tt0 = 0; tt0 < ttg.vertices(); ++tt0)
	{
	  const int symm_double = ttg.symmetric_double(tt0);
	  auto it =
	    std::find(todo_list.begin(),todo_list.end(),symm_double);
	  if (it == todo_list.end()) todo_list.push_back(tt0);
	}
    }
  else
    {
      for (int i = 0; i < ttg.vertices(); ++i)
	todo_list.push_back(i);
    }
}


template<class TrTr>
bool ttauto<TrTr>::find_pAs()
{
  using std::cout;
  using std::endl;
  using std::setw;

  // Don't keep track of transition matrices if no max dilatation is
  // specified, since we then don't use the criteria to reject
  // paths.
  if (lambdamax != 0 && do_check_norms)
    {
      // Start with identity matrix.
      TMl.push(jlt::identity_matrix<int>(n));
    }

  tryv.push(0);
  p = fpath(ttg,tt0);

  reset_counters();

  // Descend graph, letting the function increment the initial folding.
  for (int i = 0; i < ttg.foldings(tt0); ++i) while (descend_graph()) {};

  int prec = cout.precision();
  cout.precision(10);

  cout << "---------------------------------------------------\n";
  cout << "Statistics\n\n";
  cout << "Total paths tried  = " << totalpathstried << endl;
  cout << "Total path length  = " << totalpathlength << endl;
  cout << "Mean path length   = " << totalpathlength/totalpathstried << endl;
  cout << "Closed paths       = " << closed << endl;
  cout << "Primitive paths    = " << primitive << endl;
  if (do_check_gates)
    {
      cout << "Gate test          = " << gaterejected << " rejected of "
	   << gatecandidates << " candidates";
      if (gatecandidates > 0)
	cout << " (" << std::setprecision(3) << 100*gate_rejection_rate()
	     << "%)" << std::setprecision(prec);
      cout << " (cumulative)" << endl;
    }
#if 0 /* This is redundant.  We don't yet know if it's really pA. */
  cout << "pseudoAnosov paths = " << pseudoAnosov << endl;
#endif
  if (lambdamax != 0 && do_check_norms)
    {
      cout << "\n       Exceeded max norm " << setw(9);
      cout << normexceeded << " times\n";
      cout << "     Exceeded column sum " << setw(9);
      cout << colsumexceeded << " times\n";
      cout << "        Exceeded row sum " << setw(9);
      cout << rowsumexceeded << " times\n";
#ifdef TTAUTO_CHECK_SYMMETRIC_NORM
      cout << " Exceeded symmetric norm " << setw(9);
      cout << symmetricnormexceeded << " times\n";
#endif
    }
  if (max_badword_length > 0)
    {
      cout << "       Omitted bad words " << setw(9);
      cout << badwordsomitted << " times\n";
    }
  cout << "Exceeded max path length " << setw(9);
  cout << maxpathlengthexceeded << " times\n";
  cout << "---------------------------------------------------\n";

  cout.precision(prec);

  return true;
}


template<class TrTr>
bool ttauto<TrTr>::descend_graph()
{
  if (tryv.top() < p.number_of_foldings())
    {
      new_vertex();
      if (debug)
	{
	  std::cerr << "Folding...  pathlength=" << p.length();
	  std::cerr << "\tpath=" << p << "\n";
	}
    }
  else
    {
      // We've already tried all the foldings.
      return backtrack();
    }

  // Check if maximum path length exceeded.
  /* Something's bothering me here... p is already greater that
     max_path_length, so why don't I see lengths larger than
     max_path_length printed? */
  // max_path_length 0 means no bound, as max_pathlength() documents and
  // as the summary printout already assumed.
  if (max_path_length != 0 && (int)p.length() > max_path_length)
    {
      if (debug)
	std::cerr << "Exceeded max path length " << max_path_length << "\n";
      ++maxpathlengthexceeded;
      return backtrack();
    }

  // Check for bad words.  Controlled by badword_length() alone: see the
  // caveat there about what this assumes.
  if (max_badword_length > 0)
    {
      typedef typename badword_list::const_iterator cpit;

      for (int bwl = 1; bwl <= max_badword_length; ++bwl)
	{
	  // The bad word is longer than the actual path, so move on.
	  if ((int)p.length() < 2*bwl) break;

	  // The first vertex to match in the current word.
	  int vv = p.vertices()[p.length()-2*bwl];
#ifdef TTAUTO_HASH_BADWORDS
	  // The last 2*pwl foldings of the path.
	  fpath endword(p.subpath(-2*bwl));

	  // Go through list of bad words at fixed length, and check
	  // if the ending of our current path matches a bad word.
	  if (pbad(vv,bwl-1).count(endword))
	    {
	      // It matches!  Backtrack by half the length of
	      // the bad word.
	      ++badwordsomitted;
	      return backtrack(bwl);
	    }
#else
	  // Go through list of bad words at fixed length, and check
	  // if the ending of our current path matches a bad word.
	  for (cpit i = pbad(vv,bwl-1).begin();
	       i != pbad(vv,bwl-1).end(); ++i)
	    {
	      if (p.ending_equals(*i))
		{
		  // It matches!  Backtrack by half the length of
		  // the bad word.
		  ++badwordsomitted;
		  return backtrack(bwl);
		}
	    }
#endif
	}
    }

  if (lambdamax != 0 && do_check_norms && !check_all_norms())
    return backtrack();

  // Is path closed?
  if (p.closed())
    {
      if (debug) std::cerr << "Found closed path...";
      ++closed;

      // Compute transition matrix.
      //
      // If lambdamax==0, it's actually much cheaper to do this by
      // starting over, rather than keeping a running transition
      // matrix, because there are so few closed paths.
      if (lambdamax != 0 && do_check_norms)
	TM = TMl.top();
      else
	TM = p.transition_matrix();

      // Is path a pA?  The transition matrix of a pseudo-Anosov is
      // primitive (irreducible and aperiodic): its dilatation is a
      // simple dominant eigenvalue.  An irreducible but imprimitive
      // matrix (eigenvalues +-lambda, characteristic polynomial in
      // x^2) cannot come from a pseudo-Anosov; the paper's definition
      // of a primitive closed path already excludes it, and the search
      // used to admit it (issue #2).
      if (TM.is_primitive())
	{
	  if (debug) std::cerr << " primitive";
	  ++primitive;

	  // Find the characteristic polynomial.
	  charpoly = TM.charpoly();

	  // Get an upper bound on the spectral radius.
	  int colsummax = 0, rowsummax = 0;
	  for (int i = 0; i < n; ++i)
	    {
	      int colsum = 0, rowsum = 0;
	      for (int j = 0; j < n; ++j)
		{
		  colsum += TM(j,i);
		  rowsum += TM(i,j);
		}
	      if (colsum > colsummax || colsummax == 0) colsummax = colsum;
	      if (rowsum > rowsummax || rowsummax == 0) rowsummax = rowsum;
	    }
	  double lambdaupper = std::min(colsummax,rowsummax);

	  // Find the Perron root of the polynomial, using the upper
	  // bound as an initial guess.  Convexity guarantees
	  // convergence.
	  lambda = findroot(charpoly,lambdaupper,.01*tol);

	  if (lambda > lambdaupper)
	    {
	      std::cerr << "Oops!  lambda = " << lambda << " shouldn't ";
	      std::cerr << "exceed upper bound " << lambdaupper << std::endl;
	      exit(1);
	    }

#ifdef TTAUTO_USE_FORTRAN
	  // Check the root by computing the spectral radius.

	  // Copy transition matrix to a double matrix, so we can find its
	  // spectral radius..
	  jlt::mathmatrix<double> TMd(n,n);
	  for (int i = 0; i < n; ++i)
	    {
	      for (int j = 0; j < n; ++j)
		{
		  TMd(i,j) = TM(i,j);
		}
	    }
	  // Spectral radius destroys matrix, so don't use after this.
	  double lambda2 = jlt::spectral_radius(TMd);

	  if (std::abs(lambda-lambda2) > .01*tol)
	    {
	      std::cout << lambda << "\t" << lambda2;
	      std::cout << "\t" << lambda-lambda2 << std::endl;
	    }
#endif

	  // This is somewhat redundant: a primitive matrix has unit
	  // spectral radius iff it is the 1x1 matrix (1), and none of
	  // our matrices can be permutation matrices since there is at
	  // least one fold.  Nevertheless, it is very cheap to do this
	  // and we need to find the spectral radius anyways.
	  if (lambda-1 > tol) // Avoid false pA's.  Tricky.
	    {
	      if (debug) std::cerr << ", pseudo-Anosov";
	      ++pseudoAnosov;

	      // Save the pA to our list.
	      record_pA();
	    }
	}
    }
  return true;
}


template<class TrTr>
bool ttauto<TrTr>::check_all_norms()
{
  // Check column sums, row sums, norm of TM.
#ifdef TTAUTO_CHECK_SYMMETRIC_NORM
  // Also check norm of G_ij = sqrt(TM_ij TM_ji).
  // This only works for Rauzy classes, where all the matrices are of
  // the form Identity+1.
  double Gnorm = 0;
#endif
  int Mnorm = 0, colsummin = 0, rowsummin = 0;
  for (int i = 0; i < n; ++i)
    {
      int colsum = 0, rowsum = 0;
      for (int j = 0; j < n; ++j)
	{
	  colsum += TMl.top()(j,i);
	  rowsum += TMl.top()(i,j);
#ifdef TTAUTO_CHECK_SYMMETRIC_NORM
	  Gnorm += std::sqrt((double)TMl.top()(i,j)*TMl.top()(j,i));
#endif
	}
      Mnorm += colsum;
      if (Mnorm > maxnorm)
	{
	  if (debug)
	    {
	      std::cerr << "Exceeded matrix norm at pathlength ";
	      std::cerr << p.length() << "\n";
	    }
	  ++normexceeded;
	  return false;
	}
      if (colsum < colsummin || colsummin == 0) colsummin = colsum;
      if (rowsum < rowsummin || rowsummin == 0) rowsummin = rowsum;
    }

  // The minimum column sum is a lower bound on the spectral radius.
  if (colsummin > lambdamax)
    {
      if (debug)
	{
	  std::cerr << "Exceeded column sum at pathlength ";
	  std::cerr << p.length() << "\n";
	}
      ++colsumexceeded;
      return false;
    }
  // The minimum row sum is also a lower bound on the spectral radius.
  if (rowsummin > lambdamax)
    {
      if (debug)
	{
	  std::cerr << "Exceeded row sum at pathlength ";
	  std::cerr << p.length() << "\n";
	}
      ++rowsumexceeded;
      return false;
    }
#ifdef TTAUTO_CHECK_SYMMETRIC_NORM
  // Bound on the norm of the symmetrised form G_ij = sqrt(TM_ij TM_ji).
  Gnorm /= n;
  if (Gnorm > lambdamax)
    {
      if (debug)
	{
	  std::cerr << "Exceeded symmetrised norm bound at pathlength ";
	  std::cerr << p.length() << "\n";
	}
      ++symmetricnormexceeded;
      return false;
    }
#endif
  return true;
}


template<class TrTr>
inline void ttauto<TrTr>::new_vertex()
{
  using std::cout;
  using std::setw;

  ++totalpathstried;
  totalpathlength += p.length();
  if (totalpathstried % print_every == 0)
    {
      // Five digits for floating point.
      int prec = cout.precision();
      cout.precision(5);
      // Print trailing zeros.
      cout.setf(std::ios::showpoint);
      // Adjust to the right.
      cout.setf(std::ios::left,std::ios::adjustfield);
      cout << setw(11) << totalpathstried;
      cout << " (length=" << setw(3) << p.length();
      cout << " mean=" << setw(6) << totalpathlength/totalpathstried;
      cout << "  max=" << max_path_length << ")  " << p << "\n";
      cout.precision(prec);
      cout.unsetf(std::ios::showpoint);
      cout.flush();
    }

  // Copy the current fold to try, then increment it.
  int f = tryv.top()++;
  // Start with 0th fold on the new vertex.
  tryv.push(0);
  if (lambdamax != 0 && do_check_norms)
    {
      // New transition matrix.
      TMl.push(ttg.transition_matrix(p.final_vertex(),f) * TMl.top());
    }
  // Save the fold to path and find the new vertex.
  p.push_back(f);
}


// Save info about a pA class to a list.
// For a fixed lambda, save paths and transition matrices.
template<class TrTr>
inline void ttauto<TrTr>::record_pA()
{
  // Don't keep pA's above a certain dilatation, unless lambdamax=0,
  // in which case keep them all.
  if (lambda > lambdamax+tol && lambdamax != 0) return;
  // Don't keep pA's below a certain dilatation, unless lambdamin=0,
  // in which case keep them all.
  if (lambda < lambdamin-tol && lambdamin != 0) return;

  // Bestvina-Handel gate test: an irreducible matrix is not sufficient
  // (issue #2).  Rejected candidates go to a separate list.
  ++gatecandidates;
  if (do_check_gates && !p.gates().connected)
    {
      ++gaterejected;
      add_current_path(rejl);
      return;
    }

  add_current_path(pAl);
}

template<class TrTr>
inline void ttauto<TrTr>::add_current_path(pAlist& l)
{
  // Find the first polynomial not less than charpoly.
  auto pAit = l.lower_bound(charpoly);
  if (pAit == l.end() || charpoly != pAit->first)
    {
      // A new dilatation: add it to the list.
      pAclass<TrTr> newpA(charpoly,lambda,max_paths_save);
      newpA.add_path(p,TM);
      // Use pAit as a hint: for fastest results, it should point to
      // the element before the one to be inserted.  So decrement
      // pAit, unless it's already at the end.
      if (pAit != l.end()) --pAit;
      l.insert(pAit,std::pair<Poly,pAclass<TrTr> >(charpoly,newpA));
    }
  else
    {
      // Add the path to existing pAclass.
      pAit->second.add_path(p,TM);
    }
}


} // namespace ttauto

#endif // TTAUTO_TTAUTO_HPP
