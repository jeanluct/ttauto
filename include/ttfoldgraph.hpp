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

#ifndef TTFOLDGRAPH_HPP
#define TTFOLDGRAPH_HPP

#include <iostream>
#include <list>
#include <algorithm>
#include <jlt/mathmatrix.hpp>
#include <jlt/vector.hpp>
#include <jlt/csparse.hpp>
#include "mathmatrix_permplus1.hpp"

// Graph of train tracks linked by foldings.

//
// Train track automaton class.
//
// The template argument TrTr is a train track class.
//
// class TrTr must provide the following methods:
//
// int edges()				# of edges in train track
// int foldings()			max # of possible foldings
// bool operator==(const TrTr&)		isotopy of train tracks
// TrTr(const TrTr&)			copy constructor
// fold_transition_matrix(const int f)	fold and return transition matrix
// coding
// print_coding
//

namespace ttauto {

template<class TrTr> class ttfoldgraph;

template<class TrTr>
std::list<ttfoldgraph<TrTr> > subgraphs(const ttfoldgraph<TrTr>& tt);

template<class TrTr>
class ttfoldgraph
{
public:
  static const bool exploit_symmetries = TrTr::exploit_symmetries;
  typedef jlt::mathmatrix<int> Mat;
  typedef mathmatrix_permplus1 Matpp1;	// Sparse matrix type.

private:
  static const int debug = 0;

  // The train track associated with the vertex.
  jlt::vector<TrTr> trtrv;
  // Target vertex of a branch.
  jlt::vector<jlt::vector<int> > tv;
  // Transition matrix on a branch.
  jlt::vector<jlt::vector<Matpp1> > TMv;
  // Number of foldings (outgoing branches) at each vertex.
  jlt::vector<int> nfoldsv;
  /* The train track should have this info.  Ditch the sequential
     numbering of folds in traintrack?  In other words, use only
     fold(cusp,dir) from here, not fold(f)? */
  // Is the train track reflection-symmetric to another or itself?
  jlt::vector<int> symtov;
  // Is the train track cyclically symmetric?
  jlt::vector<int> symcycv;
  // Number of edges in each train track.
  const int n;
  // Maximum number of foldings at each vertex.
  const int nfoldsmax;

  Mat TM;		// Transition matrix.
  const Mat id;		// Identity matrix.

public:

  // Make a folding graph from an initial train track.
  ttfoldgraph(const TrTr& trtr)
    : n(trtr.edges()), nfoldsmax(trtr.foldings()), TM(n,n),
      id(jlt::identity_matrix<int>(n))
  {
    // Build the graph.
    add_vertex(trtr);

    if (exploit_symmetries) find_symmetries();
  }

  // Make a folding graph from nonempty vectors of train tracks,
  // target vertices, transition matrices, and number of folds.
  ttfoldgraph(const jlt::vector<TrTr>& trtrv_,
	      const jlt::vector<jlt::vector<int> >& tv_,
	      const jlt::vector<jlt::vector<Matpp1> >& TMv_,
	      const jlt::vector<int>& nfoldsv_)
    : trtrv(trtrv_), tv(tv_), TMv(TMv_), nfoldsv(nfoldsv_),
      n(trtrv.front().edges()), nfoldsmax(trtrv.front().foldings()), TM(n,n),
      id(jlt::identity_matrix<int>(n))
  {
    if (exploit_symmetries) find_symmetries();
  }

private:
  // Add a vertex.  This recursively builds the whole graph.
  int add_vertex(const TrTr& trtr)
  {
    // See if the vertex is already in the graph.
    int idx = std::distance(trtrv.begin(),
			    std::find(trtrv.begin(),trtrv.end(),trtr));

    if (idx == (int)vertices())
      {
	if (debug)
	  std::cerr << "Adding new vertex " << vertices() << std::endl;

	// It's not already in the graph, so add it.
	trtrv.push_back(trtr);
	// Add an outgoing branch.
	tv.push_back(jlt::vector<int>());
	// Add an outgoing matrix.
	TMv.push_back(jlt::vector<Matpp1>());
	// Initalise the vector giving the number of foldings.
	nfoldsv.push_back(0);
	// The index of the vertex we just added.
	idx = vertices()-1;
      }
    else
      {
	if (debug)
	  std::cerr << "Not adding vertex " << idx << std::endl;

	// It's already in the graph.  Just return its id.
	return idx;
      }

    // Try all nfolds foldings from this vertex.
    for (int f = 0; f < nfoldsmax; ++f)
      {
	// Fold and find the transition matrix.
	TrTr trtr0(trtr);
	TM = trtr0.fold_transition_matrix(f);
	if (TM != id)
	  {
	    ++nfoldsv[idx];
	    // Convert matrix to sparse type and add to list.
	    TMv[idx].push_back(Matpp1(TM));
	    // Add the target vertex, if it's not already in there, and
	    // point to it.
	    int tidx = add_vertex(trtr0);
	    tv[idx].push_back(tidx);
	  }
      }
    // Return the index, which allows the recursion.
    return idx;
  }

  // Delete a vertex.
  /* This function is no longer in use? */
  void delete_vertex(const int vv)
  {
    if (debug)
      {
	std::cerr << "Deleting vertex " << vv;
	std::cerr << " in ttfoldgraph::delete_vertex.\n";
      }

    // Erase vertex vv from graph.
    trtrv.erase(trtrv.begin()+vv);
    tv.erase(tv.begin()+vv);
    TMv.erase(TMv.begin()+vv);
    nfoldsv.erase(nfoldsv.begin()+vv);
    if (exploit_symmetries)
      {
	// If there was a vertex symmetric to this one, mark it as
	// unsymmetric.
	if (symtov[vv] != -1) symtov[symtov[vv]] = -1;
	symtov.erase(symtov.begin()+vv);
	symcycv.erase(symcycv.begin()+vv);
      }

    // Find and erase all other branches pointing to vv.
    //
    // Also adjust branches pointing to vertices > vv, since they have
    // to be decremented by one.
    //
    // Use iterators, since things might get erased as we go and it is
    // better to check against end().
    for (int v = 0; v < vertices(); ++v)
      {
	if (exploit_symmetries)
	  {
	    if (symtov[v] == vv)
	      {
		// Shouldn't happen: we did that above.
		std::cerr << "Error: we already erased this.\n";
		std::exit(1);
	      }
	    // Adjust symmetry partners past this one.
	    if (symtov[v] > vv)
	      {
		--symtov[v];
	      }
	  }

	// Loop over branches.
	jlt::vector<int>::iterator it = tv[v].begin();
	jlt::vector<Matpp1>::iterator iTM = TMv[v].begin();
	while (it != tv[v].end())
	  {
	    if (*it == vv)
	      {
		// Delete this branch.
		if (debug)
		  {
		    std::cerr << "Deleting branch pointing to deletee";
		    std::cerr << " in ttfoldgraph::delete_vertex.\n";
		  }
		// it and iTM will be pointing to the next item, so
		// don't incremement them.
		it = tv[v].erase(it);
		iTM = TMv[v].erase(iTM);
		--nfoldsv[v];
		continue;
	      }

	    // Adjust branches pointing to vertices > vv, since we
	    // have removed the vertex vv.
	    if (*it > vv) (*it)--;

	    ++it; ++iTM;
	  }
	// Possible that another vertex had only branches pointing to
	// vv?  Then delete the vertex.  But wait until after we're
	// done so as not to screw up numbering.
      }
  }

  /* This function is no longer in use? */
  void delete_vertices(std::list<int>& vl)
  {
    if (vl.empty()) return;

    // Delete the vertices in reverse order, so as to not invalidate
    // the earlier members of the list!
    std::for_each(vl.rbegin(),vl.rend(),delete_vertex);

    // Clear the list since it contains oudated numbers after deletion.
    vl.clear();
  }

  void swap_vertices(const int v1, const int v2)
  {
    if (v1 == v2) return;

    if (v1 < 0 || v1 >= vertices() || v2 < 0 || v2 >= vertices())
      {
	std::cerr << "Error: v1 or v2 out of range in swap_vertices().\n";
	exit(1);
      }

    // Swap the data for each vertex.
    std::swap(trtrv[v1],trtrv[v2]);
    std::swap(tv[v1],tv[v2]);
    std::swap(TMv[v1],TMv[v2]);
    std::swap(nfoldsv[v1],nfoldsv[v2]);
    if (exploit_symmetries)
      {
	std::swap(symtov[v1],symtov[v2]);
	std::swap(symcycv[v1],symcycv[v2]);
      }

    // Update branches pointing to v1 and v2, as well as symmetric double.
    for (int v = 0; v < vertices(); ++v)
      {
	if (exploit_symmetries)
	  {
	    if (symtov[v] == v1) symtov[v] = v2;
	    else if (symtov[v] == v2) symtov[v] = v1;
	  }

	// Loop over branches.
	jlt::vector<int>::iterator it = tv[v].begin();
	for (; it != tv[v].end(); ++it)
	  {
	    if (*it == v1) *it = v2;
	    else if (*it == v2) *it = v1;
	  }
      }
  }

  // The adjancency matrix: which edges (cols) are reachable from each
  // edge (rows).
  csparse::cs* adjacency_matrix() const
  {
    jlt::cs_auto_ptr am(csparse::cs_spalloc(0,0,1,1,1));

    for (int v = 0; v < vertices(); ++v)
      {
      	for (int b = 0; b < (int)tv[v].size(); ++b)
	  {
	    // Set row v, column tv[v][b] to 1.
	    csparse::cs_entry(am,v,tv[v][b],1);
	  }
      }

    // Return column-compress sparse matrix.
    return csparse::cs_compress(am);
  }

  void find_symmetries()
  {
    if (!exploit_symmetries)
      {
	std::cerr << "find_symmetries() error: ";
	std::cerr << "exploit_symmetries is false.\n";
	std::exit(1);
      }

    // Reflection symmetry.
    for (int v = 0; v < vertices(); ++v)
      {
	// The reverse coding.
	jlt::vector<int> rcoding = trtrv[v].coding(-1);
	// Check for symmetry with the previous tracks.
	bool issym = false;
	for (int vv = 0; vv <= v; ++vv)
	  {
	    if (trtrv[vv].coding() == rcoding)
	      {
		if (debug)
		  {
		    std::cerr << v << " symmetric to "  << vv << "!\n";
		  }
		// Make the two symmetric tracks point to each other.
		// Note that it's possible to have vv = v
		// (track is self-symmetric).
		symtov.push_back(vv);
		symtov[vv] = v;
		issym = true;
		break;
	      }
	  }
	// A -1 indicates no symmetry (yet!)
	if (!issym) symtov.push_back(-1);
      }

    // Cyclic symmetry.
    for (int v = 0; v < vertices(); ++v)
      {
	symcycv.push_back(trtrv[v].cyclic_symmetry().order()-1);
      }

    sort_by_symmetries();
  }

  void sort_by_symmetries()
  {
    if (!exploit_symmetries)
      {
	std::cerr << "sort_by_symmetries() error: ";
	std::cerr << "exploit_symmetries is false.\n";
	std::exit(1);
      }

    int frnt = 0;
    int rs = 0;
    for (int i = 0; i < (int)vertices(); ++i)
      {
	if (symtov[i] == i)
	  {
	    // Reflection-symmetric to itself: move to front.
	    swap_vertices(frnt++,i);
	  }
	else if (symtov[i] != -1)
	  {
	    // Otherwise, count the ones that are reflection-symmetric.
	    ++rs;
	  }
      }
    if ((double)rs/2 != (double)(rs/2))
      {
	std::cerr << "Error in sort_by_symmetries(): ";
	std::cerr << "odd number of reflection-symmetric doubles.\n";
	exit(1);
      }
    rs /= 2;

    // Two ways of sorting:
#if 1
    //
    // 1) Self-sym -> all left-sym -> all right-sym (default).
    //
    // The advantage of this order is that to test vertices we only
    // need to do self-sym and left-sym, which form a contiguous
    // block.
    //
    for (int i = frnt; i < frnt+rs; ++i)
      {
	if (symtov[i] != -1)
	  {
	    // Reflection-symmetric to symtov[i]: move to rs+i.
	    swap_vertices(rs+i,symtov[i]);
	  }
      }
#else
    //
    // 2) Self-sym -> left/right-sym pairs.
    //
    // No real advantage: included for completeness.
    //
    for (int i = frnt; i < (int)vertices(); )
      {
	if (symtov[i] != -1)
	  {
	    // Reflection-symmetric to symtov[i]: move to i+1.
	    swap_vertices(i+1,symtov[i]);
	    // Skip over i+1.
	    i += 2;
	  }
	else
	  {
	    ++i;
	  }
      }
#endif
  }

public:

  int vertices() const { return trtrv.size(); }

  int foldings() const { return nfoldsmax; }

  int foldings(const int v) const { return nfoldsv[v]; }

  int minfoldings() const
  {
    return *std::min_element(nfoldsv.begin(),nfoldsv.end());
  }

  int maxfoldings() const
  {
    return *std::max_element(nfoldsv.begin(),nfoldsv.end());
  }

  int edges() const { return n; }

  const Matpp1& transition_matrix(const int idx, const int br) const
  {
    return TMv[idx][br];
  }

  int target_vertex(const int idx, const int br) const
  {
    return tv[idx][br];
  }

  const TrTr& traintrack(const int idx) const { return trtrv[idx]; }

  const jlt::vector<TrTr>& traintracks() const { return trtrv; }

  const int symmetric_double(const int idx) const
  {
    if (!exploit_symmetries)
      {
	std::cerr << "symmetric_double() error: ";
	std::cerr << "exploit_symmetries is false.\n";
	std::exit(1);
      }

    return symtov[idx];
  }

  std::ostream& print(std::ostream& strm = std::cout) const
  {
    using std::endl;

    strm << "=======================================================\n";
    for (int i = 0; i < (int)vertices(); ++i)
      {
	strm << "Vertex " << i;
	if (exploit_symmetries)
	  {
	    bool anysym = false;
	    if (symtov[i] == i)
	      {
		anysym = true;
		strm << "  (reflection-sym to itself";
	      }
	    else if (symtov[i] != -1)
	      {
		anysym = true;
		strm << "  (reflection-sym to " << symtov[i];
	      }
	    if (symcycv[i])
	      {
		if (!anysym) strm << "  ("; else strm << ",";
		anysym = true;
		strm << " cyclic-sym order " << symcycv[i]+1;
	      }
	    if (anysym) strm << ")";
	  }
	strm << "\nCoding: ";
	trtrv[i].print_coding(strm) << endl;
	if (debug)
	  {
	    strm << "Rev Co: ";
	    trtrv[i].print_coding(strm,-1) << endl;
	  }
	for (int j = 0; j < nfoldsv[i]; ++j)
	  {
	    strm << "Branch " << j << " has target vertex " << tv[i][j];
	    strm << " with matrix ";
	    strm << TMv[i][j] << endl;
	  }
	strm << "=======================================================\n";
      }
    return strm;
  }

  std::ostream& printMathematicaForm(std::ostream& strm = std::cout) const
  {
    using std::endl;

    strm << "(*\nTrain track graph for n=" << trtrv.front().punctures();
    strm << " punctures and " << trtrv.front().edges() << " edges.";
    strm << "\n\nInitial train track:\n\n";
    trtrv.front().print(strm) << "\n";
    strm << "Subgraph has " << vertices();
    strm << (vertices() == 1 ? " vertex" : " vertices") << ".\n*)\n";

    // Print train track codings.
    strm << "{\n {\n";
    for (int i = 0; i < (int)vertices(); ++i)
      {
	strm << "  {\"";
	trtrv[i].print_coding(strm) << "\",";
	trtrv[i].printMathematicaForm(strm) << "}";
	if (i != (int)vertices()-1) strm << ",";
	strm << endl;
      }
    strm << " },\n {\n";

    // Print connections and transition matrices.
    for (int i = 0; i < (int)vertices(); ++i)
      {
	for (int j = 0; j < nfoldsv[i]; ++j)
	  {
	    strm << "  {" << i+1 << "->" << tv[i][j]+1 << ",";
	    TMv[i][j].printMathematicaForm(strm) << "}";
	    if (j != nfoldsv[i]-1) strm << ",\n";
	  }
	if (i != (int)vertices()-1) strm << ",";
	strm << endl;
      }
    strm << " }\n}\n";

    return strm;
  }

  friend std::list<ttfoldgraph<TrTr> >
  subgraphs<>(const ttfoldgraph<TrTr>& tt);
};


// Break up a train track graph into subgraphs.
// The main graph is the first element.
template<class TrTr>
std::list<ttfoldgraph<TrTr> > subgraphs(const ttfoldgraph<TrTr>& ttg)
{
  static const int debug = 0;
  typedef jlt::vector<int>			Vec;
  typedef jlt::vector<Vec>::const_iterator	cit;
  typedef Vec::const_iterator			vcit;
  using std::cerr;
  using std::endl;
  using std::copy;
  using std::sort;
  using std::back_inserter;

  if (debug)
    {
      cerr << "Initial graph has " << ttg.vertices() << " vertices\n";
    }

  // The adjancency matrix: which edges (cols) are reachable from each
  // edge (rows).
  // Store as a csparse triplet sparse matrix.
  jlt::cs_auto_ptr am(ttg.adjacency_matrix());

  // Dulmage-Mendelsohn decomposition of sparse matrix into
  // block-triangular form.
  jlt::csd_auto_ptr dm(csparse::cs_dmperm(am,1));
  int nb = dm->nb;

  if (debug)
    {
      cerr << "p = ";
      for (int i = 0; i < ttg.vertices(); ++i) cerr << dm->p[i] << " ";
      cerr << "\nq = ";
      for (int i = 0; i < ttg.vertices(); ++i) cerr << dm->q[i] << " ";
      cerr << "\nr = ";
      for (int i = 0; i <= nb; ++i) cerr << dm->r[i] << " ";
      cerr << endl;
    }

  // Extract blocks of vertices representing subgraphs.
  jlt::vector<Vec> blks, badrowblks, badcolblks;
  for (int i = 0; i < nb; ++i)
    {
      Vec rowblock;
      copy(dm->p+dm->r[i],dm->p+dm->r[i+1],back_inserter(rowblock));
      sort(rowblock.begin(),rowblock.end());

      Vec colblock;
      copy(dm->q+dm->r[i],dm->q+dm->r[i+1],back_inserter(colblock));
      sort(colblock.begin(),colblock.end());

      if (rowblock == colblock)
	{
	  // If the rows and columns are the same, we can easily
	  // extract this block.
	  blks.push_back(rowblock);
	}
      else
	{
	  // Otherwise we will need to treat these separately.
	  badrowblks.push_back(rowblock);
	  badcolblks.push_back(colblock);
	}
    }

  if (debug)
    {
      cerr << "Bad blocks\n";
      cit i = badrowblks.begin();
      cit j = badcolblks.begin();
      for (; i != badrowblks.end(); ++i, ++j)
	{
	  cerr << "p = ";
	  copy(i->begin(),i->end(),std::ostream_iterator<int>(cerr," "));
	  cerr << endl;

	  cerr << "q = ";
	  copy(j->begin(),j->end(),std::ostream_iterator<int>(cerr," "));
	  cerr << endl;
	}
    }

  // Deal with the bad blocks: group them together to make row/column
  // blocks that share all the same vertices.
  while (!badrowblks.empty())
    {
      // Loop through columns blocks.
      //
      // Skip the col block corresponding to the 0th row block, since
      // we want to group together indices from different blocks.
      for (int j = 1; j < (int)badcolblks.size(); ++j)
	{
	  // Does bad column block j contain an index from bad row
	  // block 0?
	  bool inthere = false;
	  for (int kk = 0; kk < (int)badcolblks[j].size(); ++kk)
	    {
	      if (std::count(badrowblks[0].begin(),badrowblks[0].end(),
			     badcolblks[j][kk]) != 0)
		{
		  inthere = true;
		  break;
		}
	    }

	  if (inthere)
	    {
	      // The 0th row vector contains the 0th element of the
	      // jth column.  We have to splice the 0 and j entries
	      // into one.
	      Vec splrow, splcol;

	      copy(badrowblks[0].begin(),badrowblks[0].end(),
		   back_inserter(splrow));
	      copy(badrowblks[j].begin(),badrowblks[j].end(),
		   back_inserter(splrow));
	      sort(splrow.begin(),splrow.end());

	      copy(badcolblks[0].begin(),badcolblks[0].end(),
		   back_inserter(splcol));
	      copy(badcolblks[j].begin(),badcolblks[j].end(),
		   back_inserter(splcol));
	      sort(splcol.begin(),splcol.end());

	      if (splrow == splcol)
		{
		  // We've formed a complete invariant block.
		  // Add to the list of blocks.
		  blks.push_back(splrow);
		  // Erase the row/col bad blocks (j>0, so erase j first).
		  badrowblks.erase(badrowblks.begin()+j);
		  badrowblks.erase(badrowblks.begin()+0);
		  badcolblks.erase(badcolblks.begin()+j);
		  badcolblks.erase(badcolblks.begin()+0);
		  break;
		}
	      else
		{
		  // The blocks don't match.  We will need to add
		  // more blocks.  Combine those two and continue.
		  badrowblks[0] = splrow;
		  badcolblks[0] = splcol;
		  badrowblks.erase(badrowblks.begin()+j);
		  badcolblks.erase(badcolblks.begin()+j);
		}
	    }
	}
    }

  // The big list of subgraphs.
  std::list<ttfoldgraph<TrTr> > ttgl;

  // Loop over the blocks.
  for (cit i = blks.begin(); i != blks.end(); ++i)
    {
      // We will extract the subgraph from the big one.
      jlt::vector<TrTr> trtrv_sub;
      jlt::vector<Vec> tv_sub;
      jlt::vector<jlt::vector<mathmatrix_permplus1> > TMv_sub;
      Vec nfoldsv_sub;

      // Loop over train tracks in this block.
      for (int j = 0; j < (int)i->size(); ++j)
	{
	  int tr = (*i)[j];

	  trtrv_sub.push_back(ttg.trtrv[tr]);
	  tv_sub.push_back(ttg.tv[tr]);
	  TMv_sub.push_back(ttg.TMv[tr]);
	  nfoldsv_sub.push_back(ttg.nfoldsv[tr]);

	  // Remove outgoing arrows with target outside subgraph.
	  // Loop backwards so we can erase as we go.
	  for (int b = ttg.nfoldsv[tr]-1; b >= 0; --b)
	    {
	      vcit tvit = std::find(i->begin(),i->end(),ttg.tv[tr][b]);
	      if (tvit == i->end())
		{
		  // The arrow is outside the graph: remove it and the
		  // associated transition matrix.
		  tv_sub[j].erase(tv_sub[j].begin() + b);
		  TMv_sub[j].erase(TMv_sub[j].begin() + b);
		  --nfoldsv_sub[j];
		}
	      else
		{
		  // The arrow is inside the subgraph:
		  // Relabel the target vertex of the branch to fit
		  // within subgraph.
		  int idx = std::distance(i->begin(),tvit);
		  tv_sub[j][b] = idx;
		}
	    }
	}
      // Add this subgraph to the list.
      ttgl.push_back(ttfoldgraph<TrTr>(trtrv_sub,tv_sub,TMv_sub,nfoldsv_sub));
    }

  return ttgl;
}

// Eliminate subgraphs consisting entirely of "multihumps" from graph list.
// These have more than 1 cusp but have only 2 foldings.
template<class TrTr>
void prune_multihumps(std::list<ttfoldgraph<TrTr> >& ttg)
{
  typedef typename std::list<ttfoldgraph<TrTr> >::iterator	ttgit;

  for (ttgit i = ++ttg.begin(); i != ttg.end(); ++i)
    {
      if (i->minfoldings() == 2  && i->maxfoldings() == 2 && i->foldings() > 2)
	{
	  ttg.erase(i--);
	}
    }
}

template<class TrTr> std::ostream&
print_subgraphs(const std::list<ttfoldgraph<TrTr> >& ttg,
		std::ostream& strm = std::cout)
{
  using std::endl;

  typedef typename std::list<ttfoldgraph<TrTr> >::const_iterator cttgit;

  int totalvertices = 0;
  strm << "===============================================\n";
  strm << "#\t#vertx\tf_min\tf_max\tcoding of initial track\n";
  strm << "-----------------------------------------------\n";
  for (cttgit i = ttg.begin(); i != ttg.end(); ++i)
    {
      totalvertices += i->vertices();
      strm << std::distance(cttgit(ttg.begin()),i)+1 << "\t";
      strm << i->vertices() << "\t";
      strm << i->minfoldings() << "\t" << i->maxfoldings() << "\t";
      i->traintrack(0).print_coding(strm);
      strm << endl;
    }
  strm << "-----------------------------------------------\n";
  strm << "total\t" << totalvertices << "\n";
  strm << "===============================================\n";

  return strm;
}

} // namespace ttauto

#endif // TTFOLDGRAPH_HPP
