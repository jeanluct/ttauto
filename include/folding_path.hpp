// <LICENSE
// LICENSE>

#ifndef FOLDING_PATH_HPP
#define FOLDING_PATH_HPP

#include <iostream>
#include <jlt/vector.hpp>
#include "path.hpp"
#include "traintracks_util.hpp"
#include "ttfoldgraph.hpp"

namespace traintracks {

template<class TrTr> class folding_path;

template<class TrTr> bool
direct_equal(const folding_path<TrTr>& p1, const folding_path<TrTr>& p2);

template<class TrTr> folding_path<TrTr>
operator*(const folding_path<TrTr>& p1, const folding_path<TrTr>& p2);

template<class TrTr>
std::ostream& operator<<(std::ostream& strm, const folding_path<TrTr>& pp);

template<class TrTr>
class folding_path
{
  static const int debug = 0;

  typedef typename ttfoldgraph<TrTr>::Mat	Mat;

  const ttfoldgraph<TrTr>* ttg;	// Pointer to train track graph.
  path fp;			// List of foldings.
  path vp;			// List of vertices.

public:
  folding_path(const ttfoldgraph<TrTr>& ttg_, const int vi_ = 0,
	       const int n_ = 0);

  // Copy constructor.
  folding_path(const folding_path& pp) : ttg(pp.ttg), fp(pp.fp), vp(pp.vp) {}

  void clear();

  folding_path& operator=(const folding_path& pp);

  int initial_vertex() const { return vp.front(); }

  void initial_vertex(const int vi_);

  // Add a fold to the path.
  void push_back(const int f);

  // Remove a fold from the path.
  void pop_back() { fp.pop_back(); vp.pop_back(); }

  path::size_type length() const { return fp.size(); }

  bool zerolength() const { return fp.empty(); }

  // Hash function for paths.
  struct hash;

  // Equality of paths.  Must check if they're closed.
  bool operator==(const folding_path& pp) const;

  bool operator!=(const folding_path& pp) const { return !operator==(pp); }

  // "Rotate" a closed path by shifting the starting point.
  void cycle_path(const int i);

  int final_vertex() const { return vp.back(); }

  bool closed() const;

  int number_of_foldings() const { return ttg->foldings(final_vertex()); }

  const path& foldings() const { return fp; }

  const path& vertices() const { return vp; }

  const Mat transition_matrix() const;

  // Return a subpath of length |nf|.
  // nf > 0 measure from initial vertex;
  // nf < 0 measure from final vertex;
  // nf = 0 return path of length 0 from initial vertex.
  folding_path subpath(int nf) const;

  // Check if a path pp matches the end of our path.
  bool ending_equals(const folding_path<TrTr>& pp) const;

  // Cycle through all allowable foldings at fixed length: 0000, 0001,
  // etc.
  bool operator++();

  // Catenate path pp to *this.
  const folding_path<TrTr>& operator*=(const folding_path<TrTr>& pp);

private:
  // Find the vertices along the path starting from the i0-th folding.
  void find_vertices(const int i0 = 0);

  // Direct equality of paths, without checking for cyclic equality
  // for closed paths.
  friend bool
  direct_equal<>(const folding_path<TrTr>& p1, const folding_path<TrTr>& p2);

  // Multiplication catenates two paths.
  // Must have (final vertex of p1) == (initial vertex of p2).
  friend folding_path<TrTr>
  operator*<>(const folding_path<TrTr>& p1, const folding_path<TrTr>& p2);

  friend std::ostream&
  operator<< <>(std::ostream& strm, const folding_path<TrTr>& pp);
};


//
// Class member definitions
//

// Hash function for paths.
//
// Just sum the vertices (but not the last one) and the foldings,
// offsetting the folding sum so it appears in higher digits.
//
// This is an invariant for cyclically equal paths, so it's good for
// distinguishing them.
//
// I also tried the 0th vertex (bad -- equal closed paths don't
// necessarly have the same), the path length (too coarse), and the
// minimum vertex (too coarse).  Nowhere near as good.
template<class TrTr>
struct folding_path<TrTr>::hash
{ 
  std::size_t operator()(const folding_path<TrTr>& pp) const
  {
    // An offset from the vertex sum.
    const int foffset = 1000000;

    if (!pp.zerolength())
      {
	// Sum all vertices but last and all foldings.
	int sumv = pp.vertices()[0] + foffset*pp.foldings()[0];
	for (int i = 1; i < (int)pp.length(); ++i)
	  {
	    sumv += pp.vertices()[i] + foffset*pp.foldings()[i];
	  }
	return sumv;
      }
    else
      return 0;
  }
};


//
// Inline method definitions
//

// Find the vertices along the path starting from the i0-th folding.
template<class TrTr>
inline void folding_path<TrTr>::find_vertices(const int i0)
{
  for (int i = i0; i < (int)length(); ++i)
    {
      if (fp[i] < 0 || fp[i] > (ttg->foldings(vp[i])-1))
	{
	  std::cerr << "Illegal folding " << fp[i] << " at vertex " << vp[i];
	  std::cerr << " in folding_path::find_vertices.\n";
	  std::exit(1);
	}
      vp[i+1] = ttg->target_vertex(vp[i],fp[i]);
    }
}

template<class TrTr>
inline folding_path<TrTr>::folding_path(const ttfoldgraph<TrTr>& ttg_,
					const int vi_, const int n_)
  : ttg(&ttg_),
    fp(0,ttg_.foldings()-1,n_),
    vp(0,ttg_.vertices()-1,n_+1)
{
  vp[0] = vi_;
  if (!zerolength()) find_vertices();
}

template<class TrTr> inline
void folding_path<TrTr>::clear()
{
  // Save initial vertex.
  int v0 = initial_vertex();
  fp.clear();
  vp.clear();
  vp.push_back(v0);
}

template<class TrTr> inline
folding_path<TrTr>& folding_path<TrTr>::operator=(const folding_path& pp)
{
  if (ttg != pp.ttg)
    {
      std::cerr << "The two paths are associated with different graphs";
      std::cerr << " in folding_path::operator=.\n";
      std::exit(1);
    }
  fp = pp.fp;
  vp = pp.vp;

  return *this;
}

template<class TrTr> inline
void folding_path<TrTr>::initial_vertex(const int vi_)
{
  if (fp.empty())
    vp.front() = vi_;
  else
    {
      std::cerr << "Can only change initial vertex for an empty path";
      std::cerr << " in folding_path::initial_vertex.\n";
      std::exit(1);
    }
}

// Add a fold to the path.
template<class TrTr> inline
void folding_path<TrTr>::push_back(const int f)
{
  if (debug && (f < 0 || f > ttg->foldings(vp.back())-1))
    {
      std::cerr << "Illegal folding " << f << " at vertex " << vp.back();
      std::cerr << " in folding_path::push_back.\n";
      std::exit(1);
    }
  fp.push_back(f);
  vp.push_back(ttg->target_vertex(vp.back(),f));
}

template<class TrTr> inline
bool folding_path<TrTr>::closed() const
{
  if (!fp.empty())
    return (initial_vertex() == final_vertex());
  else
    return false;
}


//
// Non-inline method definitions
//

// Equality of paths.  Must check if they're closed.
template<class TrTr>
bool folding_path<TrTr>::operator==(const folding_path& pp) const
{
  // They should both be open or both closed.
  if (closed() != pp.closed()) return false;

  // Are the paths open?  Then just compare foldings.
  if (!closed()) return direct_equal(*this,pp);

  // They must at least have the same length.
  if (this->length() != pp.length()) return false;

  // If the paths are closed then equality could follow by
  // cyclically permuting the path.  Check for this.

  // Copy vertex paths.
  /* Change: pop/push, or loop to length-1? */
  // Remove the last (final = initial) vertex of the closed paths.
  path vp1(vp); vp1.pop_back();
  path vp2(pp.vp); vp2.pop_back();

  // Shift the initial point of the second path to test for 'cyclic
  // equality'.
  const int k = length();
  jlt::vector<int> off0(cyclic_shift(vp1,vp2));

  // off0 nonempty means they are cyclically equivalent.
  if (!off0.empty())
    {
      // Check that the foldings are the same after offset.  Some
      // closed paths visit the same vertices, but are not equal.
      // For example, in the B_3 track all vertex paths are 0 0 0 0
      // ..., irrespective of the folding sequence.
      //
      // There could be more than one offset that leads to a match,
      // which is why off0 is a vector, so the foldings must be
      // tested from all those potential matches.
      /* Put cyclic_shift code here so that both vp and fp are
	 compared together. */
      for (int i = 0; i < (int)off0.size(); ++i)
	{
	  bool fmatch = true;
	  for (int v = 0; v < k; ++v)
	    {
	      if (!(fp[v] == pp.fp[(v+off0[i]) % k]))
		{ fmatch = false; break; }
	    }
	  if (fmatch) return true;
	}
    }

  // Not cyclically equal or, if so, foldings are different.
  return false;
}

// "Rotate" a closed path by shifting the starting point.
template<class TrTr>
void folding_path<TrTr>::cycle_path(const int i)
{
  // Can't cycle an open path.
  if (!closed() || i == 0) return;

  const int k = length();

  // Remove the last vertex.  We'll put it back later.
  vp.pop_back();

  // Use -i to get clockwise cycling.
  std::rotate(fp.begin(),fp.begin()+traintracks::mod(-i,k),fp.end());
  std::rotate(vp.begin(),vp.begin()+traintracks::mod(-i,k),vp.end());

  // Reinstate final vertex.
  vp.push_back(ttg->target_vertex(vp.back(),fp.back()));
}

template<class TrTr>
const typename folding_path<TrTr>::Mat
folding_path<TrTr>::transition_matrix() const
{
  const int n = ttg->edges();
  Mat TM(jlt::identity_matrix<int>(n));
  int v = initial_vertex();

  for (path::const_iterator i = fp.begin(); i != fp.end(); ++i)
    {
      if (*i >= ttg->foldings(v))
	{
	  std::cerr << "Illegal folding " << *i << " at vertex " << v;
	  std::cerr << " in folding_path::transition_matrix.\n";
	  std::exit(1);
	}
      TM = ttg->transition_matrix(v,*i) * TM;
      v = ttg->target_vertex(v,*i);
    }
  return TM;
}

// Return a subpath of length |nf|.
// nf > 0 measure from initial vertex;
// nf < 0 measure from final vertex;
// nf = 0 return path of length 0 from initial vertex.
template<class TrTr>
folding_path<TrTr> folding_path<TrTr>::subpath(int nf) const
{
  if ((int)length() < abs(nf)) return *this;

  if (nf < 0)
    {
      folding_path<TrTr> ep(*ttg,vp[length()+nf]);
      std::copy(fp.end()+nf,fp.end(),back_inserter(ep.fp));
      std::copy(vp.end()+nf,vp.end(),back_inserter(ep.vp));
      return ep;
    }
  else
    {
      folding_path<TrTr> ep(*ttg,vp[0]);
      std::copy(fp.begin(),fp.begin()+nf,back_inserter(ep.fp));
      std::copy(vp.begin()+1,vp.begin()+nf+1,back_inserter(ep.vp));
      return ep;
    }
}

// Check if a path pp matches the end of our path.
template<class TrTr>
bool folding_path<TrTr>::ending_equals(const folding_path<TrTr>& pp) const
{
  int nf = pp.length();

  if (nf > (int)length()) return false;

  for (int i = 0; i < nf; ++i)
    {
      if (fp[(length()-1)-i] != pp.fp[(nf-1)-i] ||
	  vp[length()-(i+1)] != pp.vp[nf-(i+1)]) return false;
    }
  return true;
}

// Cycle through all allowable foldings at fixed length: 0000, 0001,
// etc.
template<class TrTr>
bool folding_path<TrTr>::operator++()
{
  bool incr = false;

  for (int b = length()-1; b >= 0; --b)
    {
      int maxfold = ttg->foldings(vp[b])-1;

      // The sequence that branches cycle through is 0 ... maxfold.
      if (fp[b] != maxfold)
	{
	  incr = true;
	  // Increment the branch.
	  ++fp[b];
	  // Update the vertices from this new fold.
	  find_vertices(b);
	  break;
	}
      else
	{
	  // Otherwise reset branch at that position, and let
	  // the loop move on to the next (previous) position.
	  fp[b] = 0;
	}
    }

  // If nothing was changed, we're done.
  return incr;
}

template<class TrTr>
const folding_path<TrTr>&
folding_path<TrTr>::operator*=(const folding_path<TrTr>& pp)
{
  if (ttg != pp.ttg)
    {
      std::cerr << "The two paths are associated with different graphs";
      std::cerr << " in folding_path::operator*=.\n";
      std::exit(1);
    }
  if (final_vertex() != pp.initial_vertex())
    {
      std::cerr << "Second path must begin at final vertex of first path";
      std::cerr << " in folding_path::operator*=.\n";
      std::exit(1);
    }

  // Copy pp.fp to fp.
  std::copy(pp.fp.begin(),pp.fp.end(),back_inserter(fp));
  // Copy pp.vp to vp, skipping the first entry.
  std::copy(pp.vp.begin()+1,pp.vp.end(),back_inserter(vp));

  return *this;
}

//
// Friend functions
//

// Direct equality of paths: don't test for cyclic equality.
template<class TrTr> inline
bool direct_equal(const folding_path<TrTr>& p1, const folding_path<TrTr>& p2)
{
  // They must at least have the same length.
  if (p1.length() != p2.length()) return false;

  // Their initial vertices must agree.
  if (p1.initial_vertex() != p2.initial_vertex()) return false;

  // Then their folding sequences must agree.
  return (p1.fp == p2.fp);
}

template<class TrTr> inline folding_path<TrTr>
operator*(const folding_path<TrTr>& p1, const folding_path<TrTr>& p2)
{
  return folding_path<TrTr>(p1) *= p2;
}

// Print a folding path.
template<class TrTr>
std::ostream& operator<<(std::ostream& strm, const folding_path<TrTr>& pp)
{
  return (strm << pp.fp);
}

} // namespace traintracks

#endif // FOLDING_PATH_HPP
