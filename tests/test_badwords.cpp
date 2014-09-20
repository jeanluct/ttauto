// <LICENSE
// LICENSE>

#include <iostream>
#include <list>
#include <jlt/matrix.hpp>
#ifdef DOUBLEHUMP_TRACK
#include "doublehump.hpp"
#else
#include "traintrack.hpp"
#include "ttfoldgraph.hpp"
#include "ttauto.hpp"
#endif
#include "traintracks_util.hpp"
#include "folding_path.hpp"
#include "ttauto.hpp"
#include "badwords.hpp"

using namespace traintracks;

int main()
{
  using std::cout;
  using std::endl;

#ifdef DOUBLEHUMP_TRACK
  typedef doublehump TrTr;
  doublehump tt(1,3,1, 1,1);
#if 1
  // Append a 3-prong to the 1st free prong of b branch, then another 3-prong.
  tt.b.append(1,3);//->append(1,3);//->append(1,3)->append(1,3)->append(1,3);
  // Append a 3-prong to the 2nd free prong of b branch, then another 3-prong.
  //tt.b.append(2,3);//->append(1,3);
#endif

  tt.check();
#else
  typedef traintrack TrTr;
#if 1
  const int n = 6;
  traintrack tt(n,3);
  int maxplen = 5;
#else
  const int n = 5; // Ham&Song
  traintrack tt(n,3,3);
  int maxplen = 6;
#endif
#endif

  typedef std::list<folding_path<TrTr> >::const_iterator path_it;

  cout << "Train track has " << tt.punctures() << " punctures and ";
  cout << tt.edges() << " edges\n";

  ttfoldgraph<TrTr> ttg(tt);

  cout << "TT folding graph has " << ttg.vertices() << " vertices\n";

  jlt::matrix<std::list<folding_path<TrTr> > > pbad(badwords(ttg,maxplen));

  cout << "Bad words:\n";

  for (int plen = 0; plen < maxplen; ++ plen)
    {
      for (int v = 0; v < ttg.vertices(); ++v)
	{
	  for (path_it i = pbad(v,plen).begin(); i != pbad(v,plen).end(); ++i)
	    {
	      cout << "Vertex " << v << " " << "\tfolding path " << *i;
	      cout << "\t(vertices " << i->vertices() << ")\n";
	    }
	}
    }

  return 0;
}
