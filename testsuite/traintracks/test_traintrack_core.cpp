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

// Core train-track invariants on every stratum for n = 3, 4, 5 and every
// one-fold neighbour: check() after folding, idempotent normalisation
// (coding and numbering unchanged), coding round trip, symmetry accessors,
// puncture labels, and the text output routines.
//
// The string constructor traintrack(const char*) is deliberately not
// exercised: it exits with "Broken?" (lib/traintracks/traintrack.cpp) and
// reads an old, unlabelled coding.  Its intended behaviour is covered here
// by parsing print_coding() output back into a coding vector.

#include <iostream>
#include <sstream>
#include <string>
#include <jlt/vector.hpp>
#include "traintracks/build.hpp"
#include "traintracks/coding.hpp"
#include "traintracks/traintrack.hpp"
#include "check.hpp"

using traintracks::multigon;
using traintracks::traintrack;
using traintracks::ttnumbering;

// Parse print_coding() output (one-indexed 5-digit blocks) back into the
// zero-based coding vector, as traintrack(const char*) was meant to.
static traintrack::intVec parse_printed_coding(const std::string& s)
{
  CHECK(traintrack::label_multigons);
  traintrack::intVec code;
  std::string digits;
  for (char c : s) if (c >= '0' && c <= '9') digits += c;
  CHECK(digits.size() % 5 == 0);
  for (size_t i = 0; i < digits.size(); i += 5)
    {
      code.push_back(digits[i]-'0'-1);
      code.push_back(digits[i+1]-'0');
      code.push_back(digits[i+2]-'0'-1);
      code.push_back(digits[i+3]-'0'-1);
      code.push_back(digits[i+4]-'0');
    }
  return code;
}

// Const view, so Multigon() resolves to the public accessor.
static const multigon& MG(const traintrack& tt, const int m) { return tt.Multigon(m); }

static int check_track(const traintrack& tt)
{
  tt.check();

  // Copy equals original; coding round trip.
  traintrack tt_copy(tt);
  CHECK(tt == tt_copy);
  CHECK(tt.coding() == tt_copy.coding());
  traintrack tt_from_code(tt.coding());
  CHECK(tt == tt_from_code);
  CHECK(tt_from_code.coding() == tt.coding());

  // The reversed coding is the canonical coding of the mirror image: it
  // builds the mirror track, which equals tt iff tt is reflection symmetric.
  traintrack tt_from_rev(tt.coding(-1));
  CHECK(tt_from_rev.coding() == tt.coding(-1));
  CHECK(tt_from_rev.coding(-1) == tt.coding());
  CHECK((tt_from_rev == tt) == tt.is_reflection_symmetric());

  // Normalisation is idempotent on coding and numbering.
  const traintrack::intVec code0 = tt.coding();
  const ttnumbering N0 = tt.numbering();
  for (int i = 0; i < 3; ++i)
    {
      tt_copy.normalise();
      tt_copy.check();
      CHECK(tt_copy.coding() == code0);
      CHECK(tt_copy.numbering() == N0);
    }

  // Symmetry accessors agree with the coding.
  CHECK(tt.is_reflection_symmetric() == (tt.coding(1) == tt.coding(-1)));
  traintrack tt_sym(tt);
  const int cs = tt_sym.is_cyclically_symmetric();
  CHECK(cs >= 0);
  CHECK(tt_sym == tt);

  // Counts are consistent.
  int prongs = 0, punct = 0, mono = 0;
  for (int m = 0; m < tt.multigons(); ++m)
    {
      prongs += tt.Multigon(m).prongs();
      if (tt.Multigon(m).punctured()) ++punct;
      if (tt.Multigon(m).prongs() == 1) ++mono;
    }
  CHECK(prongs == tt.total_prongs());
  CHECK(punct == tt.punctures());
  CHECK(mono == tt.monogons());
  CHECK(tt.foldings() == 2*tt.cusps());
  CHECK((int)tt.weights().size() == tt.edges());
  CHECK(N0.nedges() == tt.edges());
  CHECK(N0.nprongs() == tt.total_prongs());
  CHECK(N0.ncusps() == tt.cusps());

  // Printed coding parses back to the coding vector.
  std::ostringstream oc;
  tt.print_coding(oc);
  CHECK(parse_printed_coding(oc.str()) == tt.coding());
  std::ostringstream ocr;
  tt.print_coding(ocr,-1);
  CHECK(parse_printed_coding(ocr.str()) == tt.coding(-1));

  // Other printers produce something sensible.
  std::ostringstream os;
  tt.print(os);
  CHECK(os.str().find("multigon 0 is a") == 0);
  std::ostringstream osd;
  tt.print_singularity_data(osd);
  CHECK(osd.str().find("(" + std::to_string(tt.cusps()) + ")") != std::string::npos);
  std::ostringstream om;
  traintracks::printMathematicaForm(om,tt);
  CHECK(!om.str().empty());

  // Weights set through the iterator are read back in the same order.
  traintrack tt_w(tt);
  jlt::vector<double> w(tt.edges());
  for (int i = 0; i < tt.edges(); ++i) w[i] = i + 1;
  tt_w.weights(w.begin());
  CHECK(tt_w.weights() == w);

  return cs;
}

int main()
{
  int ntracks = 0, nfolds = 0, nlegal = 0, ncyclic = 0, nreflect = 0;

  for (int n = 3; n <= 5; ++n)
    {
      jlt::vector<traintrack> ttv = traintracks::build_traintrack_list(n);
      CHECK(!ttv.empty());
      for (const traintrack& tt0 : ttv)
        {
          CHECK(tt0.punctures() == n);
          if (check_track(tt0) > 0) ++ncyclic;
          if (tt0.is_reflection_symmetric()) ++nreflect;
          ++ntracks;

          // Every fold index: legal folds give a valid normalised track;
          // illegal folds leave the track unchanged.
          for (int f = 0; f < tt0.foldings(); ++f)
            {
              traintrack tt(tt0);
              const bool ok = tt.fold(f);
              ++nfolds;
              if (!ok) { CHECK(tt == tt0); continue; }
              ++nlegal;
              CHECK(tt.edges() == tt0.edges());
              CHECK(tt.total_prongs() == tt0.total_prongs());
              CHECK(tt.punctures() == tt0.punctures());
              if (check_track(tt) > 0) ++ncyclic;
              if (tt.is_reflection_symmetric()) ++nreflect;
              ++ntracks;
            }
        }
    }
  CHECK(nlegal > 0);
  // (On these strata every fold index is legal; illegal ones are still
  // handled above when they occur.)

  // Puncture labels: set_label() changes the coding of a punctured monogon
  // and pure_braid() gives every puncture a distinct label.
  if (traintrack::label_multigons)
    {
      traintrack tt(5,3);
      const traintrack::intVec code0 = tt.coding();
      traintrack tl(tt);
      int m = 0;
      while (!MG(tl,m).punctured()) ++m;
      tl.set_label(m,7);
      tl.check();
      CHECK(tl.coding() != code0);
      CHECK(!(tl == tt));
      bool found = false;
      for (int k = 0; k < tl.multigons(); ++k)
        if (MG(tl,k).label() == 7) found = true;
      CHECK(found);
      traintrack tl_rt(tl.coding());
      CHECK(tl_rt == tl);

      traintrack tp(tt);
      tp.pure_braid();
      tp.check();
      std::vector<int> labels;
      for (int k = 0; k < tp.multigons(); ++k)
        {
          const int lb = MG(tp,k).label();
          if (MG(tp,k).punctured())
            {
              CHECK(lb >= 1 && lb <= tp.punctures());
              for (int l : labels) CHECK(l != lb);
              labels.push_back(lb);
            }
          else
            CHECK(lb == 0);
        }
      CHECK((int)labels.size() == tp.punctures());
      traintrack tp_rt(tp.coding());
      CHECK(tp_rt == tp);
    }

  // The original nine-puncture fixture.
  {
    jlt::vector<int> Kv(3);
    Kv[0] = 4; Kv[1] = 5; Kv[2] = 3;
    traintrack tt(9,Kv);
    if (traintrack::label_multigons) tt.set_label(0,8);
    check_track(tt);
    CHECK(tt.punctures() == 9);
  }

  std::cout << "test_traintrack_core: OK (" << ntracks << " tracks, "
            << nlegal << " legal of " << nfolds << " fold indices, "
            << ncyclic << " cyclically symmetric, " << nreflect
            << " reflection symmetric)\n";
  return 0;
}
