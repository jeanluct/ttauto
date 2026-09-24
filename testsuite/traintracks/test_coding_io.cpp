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

// Printing and reading codings: traintrack::print_coding,
// traintracks::parse_coding and traintrack(const char*).
//
// A block is four fields, (prong, nprongs, edge, nedges), and gains the
// label as a third field when some multigon carries one.  Fields run
// together as digits while each is a single digit, and are separated by
// '-' otherwise.  The width is carried per block, so reading never has to
// infer it.
//
// Two things here are not covered anywhere else.  The hyphenated form
// needs a coding field to reach ten, which no other test comes close to:
// block fields are local (nprongs is one multigon's prong count, nedges
// the edges at one prong), so they stay small even at seven punctures.
// And the parser's rejections are reachable from user input, unlike the
// library's other fail-fast paths, so they are worth testing; since they
// report by exiting, that has to happen in a child process.

#include <algorithm>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <jlt/vector.hpp>
#include "traintracks/build.hpp"
#include "traintracks/coding.hpp"
#include "traintracks/traintrack.hpp"
#include "check.hpp"

#if defined(__unix__) || defined(__APPLE__)
#define TTAUTO_CAN_FORK 1
#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

using traintracks::build_traintrack_list;
using traintracks::parse_coding;
using traintracks::traintrack;

// Fields in the first block of a printed coding: digits when the blocks
// run together, hyphen-separated parts otherwise.
static int count_fields(const std::string& s)
{
  const std::string tok = s.substr(0,s.find(' '));
  if (tok.find('-') == std::string::npos) return (int)tok.size();
  return 1 + (int)std::count(tok.begin(),tok.end(),'-');
}

static bool is_labelled(const traintrack& tt)
{
  for (int m = 0; m < tt.multigons(); ++m)
    if (tt.Multigon(m).label() != 0) return true;
  return false;
}

static bool hyphenated(const std::string& s)
{
  return s.find('-') != std::string::npos;
}

// Every printed form of one track reads back to the same coding, and
// rebuilds the same track.
static void check_forms(const traintrack& tt)
{
  std::ostringstream os;
  tt.print_coding(os);
  const std::string s = os.str();
  CHECK_MSG(parse_coding(s) == tt.coding(), s);
  CHECK_MSG(traintrack(s.c_str()) == tt, s);

  // Forcing the label always gives five fields, and loses nothing.
  std::ostringstream of;
  tt.print_coding(of,1,true);
  const std::string f = of.str();
  CHECK_MSG(parse_coding(f) == tt.coding(), f);
  CHECK_MSG(traintrack(f.c_str()) == tt, f);
  CHECK_MSG(count_fields(f) == 5, f);

  // The two widths of one track carry the same coding.
  CHECK_MSG(parse_coding(s) == parse_coding(f), s + " / " + f);

  // The reversed walk codes the mirror image, so it matches coding(-1)
  // rather than the track itself.
  std::ostringstream orv;
  tt.print_coding(orv,-1);
  CHECK_MSG(parse_coding(orv.str()) == tt.coding(-1), orv.str());

  // Without forcing, the width follows the labels.
  CHECK_MSG(count_fields(s) == (is_labelled(tt) ? 5 : 4), s);
}

#ifdef TTAUTO_CAN_FORK
// Build a track from coding in a child process, returning the child's
// exit status and whatever it wrote to stderr.  parse_coding() and
// recursive_build() report bad input by exiting, as the library does
// everywhere, so letting them is the only way to test what they reject.
static int run_isolated(const char* coding, std::string& err)
{
  err.clear();
  std::cout.flush();
  std::cerr.flush();

  int fd[2];
  CHECK(pipe(fd) == 0);

  const pid_t pid = fork();
  CHECK(pid >= 0);

  if (pid == 0)
    {
      // Child: stderr down the pipe, stdout discarded.
      close(fd[0]);
      dup2(fd[1],STDERR_FILENO);
      close(fd[1]);
      const int devnull = open("/dev/null",O_WRONLY);
      if (devnull >= 0) { dup2(devnull,STDOUT_FILENO); close(devnull); }
      traintrack tt(coding);
      tt.print_coding(std::cout);
      _exit(0);
    }

  close(fd[1]);
  char buf[512];
  ssize_t nr;
  while ((nr = read(fd[0],buf,sizeof(buf))) > 0) err.append(buf,(size_t)nr);
  close(fd[0]);

  int status = 0;
  CHECK(waitpid(pid,&status,0) == pid);
  return WIFEXITED(status) ? WEXITSTATUS(status) : -1;
}
#endif

int main()
{
  using std::cout;
  using std::endl;

  //
  // Round trips, every stratum and every printed form.
  //
  int ntracks = 0;
  for (int n = 3; n <= 11; ++n)
    {
      const jlt::vector<traintrack> ttv = build_traintrack_list(n);
      for (int t = 0; t < (int)ttv.size(); ++t)
	{
	  check_forms(ttv[t]);

	  // pure_braid() labels every puncture, which forces the wider
	  // block and, past nine punctures, the hyphenated form.
	  traintrack tp(ttv[t]);
	  tp.pure_braid();
	  CHECK(is_labelled(tp));
	  check_forms(tp);

	  ntracks += 2;
	}
    }
  cout << "round trips: " << ntracks << " tracks, n = 3 to 11" << endl;

  //
  // The hyphenated form, reached two independent ways.
  //
  {
    // A prong carrying ten edges: the eleven-puncture star track.
    const traintrack tt = build_traintrack_list(11)[0];
    std::ostringstream os;
    tt.print_coding(os);
    CHECK_MSG(hyphenated(os.str()), os.str());
    CHECK(count_fields(os.str()) == 4);
    CHECK(traintrack(os.str().c_str()) == tt);

    // The coding really does hold a field that no digit could express.
    const jlt::vector<int> code = tt.coding();
    bool big = false;
    for (int i = 0; i < (int)code.size(); i += 5)
      if (code[i]+1 > 9 || code[i+1] > 9 || code[i+3]+1 > 9 || code[i+4] > 9)
	big = true;
    CHECK(big);
    cout << "star track: " << os.str() << endl;
  }
  {
    // A label of ten, on a track small enough that nothing else is big.
    traintrack tt(5,3);
    // Const view, so Multigon() resolves to the public accessor.
    const traintrack& ctt = tt;
    int m = 0;
    while (!ctt.Multigon(m).punctured()) ++m;

    traintrack t8(tt);
    t8.set_label(m,8);
    std::ostringstream o8;
    t8.print_coding(o8);
    CHECK_MSG(!hyphenated(o8.str()), o8.str());
    CHECK(count_fields(o8.str()) == 5);
    CHECK(traintrack(o8.str().c_str()) == t8);

    traintrack t9(tt);
    t9.set_label(m,9);
    std::ostringstream o9;
    t9.print_coding(o9);
    CHECK_MSG(hyphenated(o9.str()), o9.str());
    CHECK(count_fields(o9.str()) == 5);
    CHECK(traintrack(o9.str().c_str()) == t9);
    cout << "label 9:    " << o9.str() << endl;
  }

  //
  // The published four-wide coding of ttauto_n=6_1.m.  This is the exact
  // string whose forty digits are divisible by five as well as four, so
  // any reader that counts globally takes it for eight five-digit blocks;
  // reading the width off each block cannot.
  //
  {
    const char* pub = "1111 1115 1125 1111 1135 1111 1145 1111 1155 1111";
    const traintrack tt(pub);
    CHECK(tt.punctures() == 6);
    CHECK(tt.edges() == 5);
    CHECK(!is_labelled(tt));
    std::ostringstream os;
    tt.print_coding(os);
    CHECK_MSG(os.str() == pub, os.str());
    CHECK(tt == build_traintrack_list(6)[0]);
  }

  //
  // What the reader rejects.  Every rejection exits 1, so the message is
  // checked too: the status alone would not show the right test fired.
  //
#ifdef TTAUTO_CAN_FORK
  {
    // A good coding must come through the same harness untouched, or the
    // cases below would pass for the wrong reason.
    std::string err;
    CHECK(run_isolated("1111 1115 1125 1111 1135 1111 1145 1111 1155 1111",
		       err) == 0);
    CHECK_MSG(err.empty(), err);
  }

  struct bad_case { const char* coding; const char* expect; };
  const bad_case bad[] = {
    // One integer per token: the form whose width could only be guessed.
    { "1 1 1 1 1 3 2 2",                  "expected 4 or 5" },
    // Blocks run together with no separator at all.
    { "1111111151125111113511111145111115511111 1111",
                                          "expected 4 or 5" },
    { "1111111151125111113511111145111115511111",
                                          "odd number of blocks" },
    // Widths that disagree between blocks.
    { "11111 11115 1125 11111",           "the coding is 5 wide" },
    // A prong or edge outside its own count.
    { "5111 1111",                        "prong 5 of 1" },
    { "1121 1111",                        "edge 2 of 1" },
    // A label is printed one-based, so zero cannot appear.
    { "11011 11011",                      "has label 0" },
    // Structurally short.
    { "1111 1115 1125",                   "odd number of blocks" },
    { "",                                 "no coding blocks" },
    // Even, per-block valid, but the walk runs out part way through.
    { "1111 1115 1125 1111 1135 1111",    "middle of the walk" },
    // Not a number at all.
    { "111a 1111",                        "is not a number" },
    { "1-1-x-1 1-1-1-1",                  "is not a number" },
  };

  for (int i = 0; i < (int)(sizeof(bad)/sizeof(bad[0])); ++i)
    {
      std::string err;
      const int status = run_isolated(bad[i].coding,err);
      CHECK_MSG(status == 1, std::string("accepted \"") + bad[i].coding + "\"");
      CHECK_MSG(err.find(bad[i].expect) != std::string::npos,
		std::string("\"") + bad[i].coding + "\" gave: " + err);
    }
  cout << "rejected:    " << (int)(sizeof(bad)/sizeof(bad[0]))
       << " malformed codings, each with the expected message" << endl;
#else
  cout << "rejection cases skipped: no fork() on this platform" << endl;
#endif

  cout << "\ntest_coding_io: all checks passed" << endl;
  return 0;
}
