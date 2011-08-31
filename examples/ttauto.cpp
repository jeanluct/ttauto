#include <iostream>
#include <sstream>
#include <string>
#include <jlt/vector.hpp>
#include <jlt/prompt.hpp>
#include "traintrack.hpp"
#include "ttfoldgraph.hpp"
#include "ttauto.hpp"

using namespace traintracks;

int main()
{
  using std::cout;
  using std::endl;

  typedef ttfoldgraph<traintrack>			ttgraph;
  typedef std::list<ttgraph>::iterator			ttgit;
  typedef std::list<ttgraph>::const_iterator		cttgit;
  typedef jlt::vector<traintrack>			ttVec;
  typedef ttVec::const_iterator				ttVeccit;

  std::string svnId("$Id$");
  std::string svnVersion("$LastChangedDate$");
  std::string svnDate("$LastChangedRevision$");

  cout << "-------------------------------------------------------------\n";
  cout << "---------------- ttauto r4 (2011-08-31) ALPHA ---------------\n";
  cout << "-------------------------------------------------------------\n";
#define TTAUTO_RELEASE
#ifdef TTAUTO_RELEASE
  cout << endl;
  cout << "Copyright 2010,2011 Jean-Luc Thiffeault (jeanluc@mailaps.org)\n";
  cout << "                    Erwan Lanneau       (lanneau@cpt.univ-mrs.fr)\n";
  cout << endl;
  cout << "This program is a preliminary ALPHA version and should\n";
  cout << "be used with caution.  The authors make no guarantees\n";
  cout << "regarding the validity of the results.\n";
  cout << endl;
#endif

  int n = jlt::read_number("\nNumber of punctures",3,20,5);

  ttVec ttv = ttbuild_list(n);
  cout << "\nList of strata:\n";
  for (ttVeccit i = ttv.begin(); i != ttv.end(); ++i)
    {
      cout << std::setw(2) << std::distance(ttVeccit(ttv.begin()),i)+1;
      cout << ")\t";
      i->print_singularity_data(cout) << endl;
    }
  // Which train track to search.
  int trk = jlt::read_number("\nChoose stratum",1,(int)ttv.size(),1);
  cout << "\nSelected stratum ";
  --trk;
  ttv[trk].print_singularity_data(cout) << endl;

  cout << "\nTrain track has " << ttv[trk].punctures() << " punctures and ";
  cout << ttv[trk].edges() << " edges\n";

  cout << "\nConstructing automaton...\n";

  ttgraph ttg0(ttv[trk]);
  cout << "\nFolding automaton has " << ttg0.vertices() << " vertices\n";

  std::list<ttgraph> ttg;

  bool div = jlt::yesno("\nDivide into invariant subgraphs",true);
  if (div)
    {
      cout << "\nDividing automaton into invariant subgraphs...\n";
      // Make a list of train track graphs. The main graph is the 1st element.
      std::list<ttgraph> ttgdiv(subgraphs(ttg0));
      // Copy to master list.
      std::copy(ttgdiv.begin(),ttgdiv.end(),back_inserter(ttg));
    }
  else
    {
      // Just add a single graph to the master list.
      ttg.push_back(ttg0);
    }

  cout << "\nSubgraphs:\n";
  print_subgraphs(ttg);

  if (jlt::yesno("\nWrite subgraphs to file",false))
    {
      // Print subgraphs to file.
      // Make filename.
      cout << "\nWriting subgraphs to Mathematica file ";
      std::ostringstream ostr;
      ostr << "_n=" << n << "_" << trk+1;
      if (div) ostr << "_inv";
      std::string ttautofile = "ttauto" + ostr.str() + ".m";
      cout << ttautofile << "...\n";
      std::ofstream ttastrm(ttautofile.c_str());
      ttastrm << "{\n";
      for (cttgit i = ttg.begin(); i != ttg.end(); ++i)
	{
	  i->printMathematicaForm(ttastrm);
	  if (std::distance(i,cttgit(ttg.end())) > 1) ttastrm << ",\n";
	}
      ttastrm << "}\n";
      ttastrm.close();
    }

  int fg = jlt::read_number("\nSubgraph to search",1,(int)ttg.size(),1);
  --fg;
  int len = jlt::read_number("\nMax path length in graph",0,100,5);
  if (len == 0) exit(0);
  cttgit ifg = ttg.begin();
  std::advance(ifg,fg);		// Go to the correct position.

  // Make filename.
  std::ostringstream ostr;
  ostr << "_n=" << n << "_" << trk+1;
  if (div) ostr << "_" << fg+1 << "_inv";
  std::string pAfile = "pA" + ostr.str() + ".m";

  ttauto<traintrack> tta(*ifg);
  tta.output_file(pAfile.c_str());
  tta.max_pathlength(len);
  tta.search();

  cout << "\nThe pseudo-Anosov candidates are in file " << pAfile << endl;
}
