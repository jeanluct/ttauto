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

// Pins what the pruning tests of ttauto::search do, so that restructuring
// them cannot change it unnoticed (issue #23).
//
// In norm-bounded mode (max_dilatation + check_norms) the search abandons
// a path when the norm of its matrix exceeds the Ham-Song bound, or its
// smallest column sum or row sum exceeds the window.  For each search
// below this records how many paths it tried, how many each test
// abandoned, and the characteristic polynomials of the classes it found,
// read from the statistics the search prints, and checks that
// ttauto::pruned() agrees with them.  With check_OstSch(false) the
// expected values are those of the code before the restructuring of issue
// #23.  The same searches are then run with the Ostrowski-Schneider test
// on, the default: it is safe, so each must find exactly the classes of
// its twin without it, and it must try no more paths.

#include <iostream>
#include <list>
#include <set>
#include <sstream>
#include <string>
#include "check.hpp"
#include "traintracks/build.hpp"
#include "traintracks/traintrack.hpp"
#include "ttauto/ttauto.hpp"
#include "ttauto/ttfoldgraph.hpp"

using traintracks::traintrack;
typedef ttauto::ttfoldgraph<traintrack> ttgraph;

struct counts
{
  long long tried = 0, HamSong = 0, colsum = 0, rowsum = 0, OstSch = 0;
  // The same four, from ttauto::pruned() rather than the printout.
  long long pHamSong = 0, pcolsum = 0, prowsum = 0, pOstSch = 0;
  int headings = 0;		// Printed "Norms exceeded:" headings.
  int OstSchlines = 0;		// Printed Ostrowski-Schneider lines.
  std::set<std::string> polys;
};

// The number after the label on a statistics line, or -1.
static long long after(const std::string& line, const std::string& label)
{
  const std::size_t k = line.find(label);
  if (k == std::string::npos) return -1;
  std::istringstream in(line.substr(k + label.size()));
  std::string eq;
  long long v = -1;
  if (label == "Total paths tried") in >> eq;   // "= N"
  in >> v;
  return v;
}

// Search every subgraph of stratum s of n punctures, window Lambda.  The
// statistics block is printed once per initial vertex; the counts add up.
static counts run(const int n, const int s, const double lam, const bool OstSch)
{
  jlt::vector<traintrack> ttv = traintracks::build_traintrack_list(n);
  std::list<ttgraph> sg(ttauto::subgraphs(ttgraph(ttv[s])));
  counts c;
  for (const ttgraph& g : sg)
    {
      ttauto::ttauto<traintrack> tta(g);
      tta.max_dilatation(lam).check_norms().check_OstSch(OstSch);
      std::ostringstream out;
      std::streambuf* saved = std::cout.rdbuf(out.rdbuf());
      tta.search();
      std::cout.rdbuf(saved);
      typedef ttauto::ttauto<traintrack> tt;
      c.pHamSong += tta.pruned(tt::prune_HamSong);
      c.pcolsum += tta.pruned(tt::prune_colsum);
      c.prowsum += tta.pruned(tt::prune_rowsum);
      c.pOstSch += tta.pruned(tt::prune_OstSch);

      std::istringstream in(out.str());
      std::string line;
      while (std::getline(in,line))
        {
          long long v;
          if ((v = after(line,"Total paths tried")) >= 0) c.tried += v;
          if (line == "Norms exceeded:") ++c.headings;
          if ((v = after(line,"Ham-Song norm")) >= 0) c.HamSong += v;
          if ((v = after(line,"column sum")) >= 0) c.colsum += v;
          if ((v = after(line,"row sum")) >= 0) c.rowsum += v;
          if ((v = after(line,"Ostrowski-Schneider")) >= 0)
            { c.OstSch += v; ++c.OstSchlines; }
        }
      for (auto it = tta.pA_list().begin(); it != tta.pA_list().end(); ++it)
        {
          std::ostringstream p;
          p << it->first;
          c.polys.insert(p.str());
        }
    }
  return c;
}

static void print(const char* what, const int n, const int s,
                  const double lam, const counts& c)
{
  std::cout << "n=" << n << " stratum " << s << " Lambda " << lam << " "
            << what << ": tried " << c.tried << ", abandoned by Ham-Song "
            << c.HamSong << ", column sum " << c.colsum << ", row sum "
            << c.rowsum << ", Ostrowski-Schneider " << c.OstSch << ", "
            << c.polys.size() << " classes\n";
}

static std::string join(const std::set<std::string>& polys)
{
  std::string joined;
  for (const std::string& p : polys)
    joined += (joined.empty() ? "" : "; ") + p;
  return joined;
}

int main()
{
  // Classes are listed by characteristic polynomial, in the order the
  // set sorts their printed forms, joined by "; ".  The first counts are
  // without the Ostrowski-Schneider test, the last five with it.
  struct expected { int n, s; double lam;
                    long long tried, HamSong, colsum, rowsum;
                    const char* polys;
                    long long wtried, wHamSong, wcolsum, wrowsum, wOstSch; };
  const expected cases[] = {
    { 3, 0,  5.0,     206,     50,    30,     24,
      "x^2 - 3 x + 1; x^2 - 4 x + 1; x^2 - 5 x + 1",
         142,    28,     6,    20,    18 },
    { 3, 0, 10.0,    1090,    176,   200,    170,
      "x^2 - 10 x + 1; x^2 - 3 x + 1; x^2 - 4 x + 1; x^2 - 5 x + 1; "
      "x^2 - 6 x + 1; x^2 - 7 x + 1; x^2 - 8 x + 1; x^2 - 9 x + 1",
         766,   102,    88,   144,    50 },
    { 4, 0,  2.5,  216840, 142240,   608,   2292, "",
       54672, 16154,   256,  1894, 18152 },
    { 4, 1,  2.5,    8838,   2839,   343,   1239,
      "x^4 - 2 x^3 - 2 x + 1",
        3504,   245,    22,   837,   650 },
    { 4, 1,  3.0,   95602,  29395,  4082,  14326,
      "x^4 - 2 x^3 - 2 x + 1; x^4 - 2 x^3 - 2 x^2 - 2 x + 1",
       18212,  1264,     0,  4348,  3496 },
    { 5, 2,  2.0,   17920,   6945,   245,   1772, "",
        5264,   865,     0,   606,  1163 },
    { 5, 3,  2.0, 2557432, 934747, 47112, 296862, "",
      301428, 25432,     0, 44501, 80786 },
  };
  for (const expected& e : cases)
    {
      // Without the Ostrowski-Schneider test: the search as it was.
      const counts c = run(e.n,e.s,e.lam,false);
      print("without",e.n,e.s,e.lam,c);
      for (const std::string& p : c.polys) std::cout << "    " << p << "\n";
      CHECK(c.tried == e.tried);
      CHECK(c.HamSong == e.HamSong);
      CHECK(c.colsum == e.colsum);
      CHECK(c.rowsum == e.rowsum);
      // pruned() counts the whole search, all initial vertices together,
      // so it must equal the printed per-vertex counts added up.
      CHECK(c.pHamSong == c.HamSong);
      CHECK(c.pcolsum == c.colsum);
      CHECK(c.prowsum == c.rowsum);
      // Switched off, the test neither fires nor prints a line.
      CHECK(c.pOstSch == 0);
      CHECK(c.OstSchlines == 0);
      CHECK(c.headings > 0);
      CHECK_MSG(join(c.polys) == e.polys, join(c.polys));

      // With it, as by default.
      const counts d = run(e.n,e.s,e.lam,true);
      print("with",e.n,e.s,e.lam,d);
      CHECK(d.tried == e.wtried);
      CHECK(d.HamSong == e.wHamSong);
      CHECK(d.colsum == e.wcolsum);
      CHECK(d.rowsum == e.wrowsum);
      CHECK(d.OstSch == e.wOstSch);
      CHECK(d.pHamSong == d.HamSong);
      CHECK(d.pcolsum == d.colsum);
      CHECK(d.prowsum == d.rowsum);
      CHECK(d.pOstSch == d.OstSch);
      CHECK(d.OstSchlines > 0);
      CHECK(d.headings == c.headings);
      CHECK(d.tried <= c.tried);
      // Safe: the same classes as without it.
      CHECK_MSG(join(d.polys) == e.polys, join(d.polys));
    }
  std::cout << "\ntest_prune_counts: all checks passed" << std::endl;
  return 0;
}
