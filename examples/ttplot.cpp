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


// Draw a train track from its coding, as TikZ.
//
// This is the "collapsed" representation of the ttauto paper, Fig. 12(b):
// every multigon is shrunk to a point and only the main edges are drawn,
// which the paper calls "a simplified, simply-connected version" of the
// track.  Style follows the paper's figures: plain black line art, no
// outer boundary, sparse italic labels.
//
// The layout is traintracks::make_collapsed_layout: punctures lie on an
// invisible horizontal axis in the order the boundary walk meets them,
// the track sits in the upper half plane, and the edges are routed as an
// arc diagram over the axis.  testsuite/traintracks/test_collapsed_layout
// checks every vertex of every automaton for n = 3..7, 3700 tracks, has
// no crossing.  For n = 8 some tracks still cross.
//
// The full representation of Fig. 12(a), with the monogon loops and the
// polygon sides drawn out, is not implemented here.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "traintracks/coding.hpp"
#include "traintracks/collapsed_layout.hpp"
#include "traintracks/embedding.hpp"
#include "traintracks/traintrack.hpp"

namespace {

using traintracks::collapsed_layout;
using traintracks::cubic;
using traintracks::multigon;
using traintracks::traintrack;
using traintracks::ttnumbering;
using traintracks::tt_embedding;
using traintracks::vec2;

struct options
{
  std::string coding;
  std::string coding_file;
  std::string output_file = "ttplot.tex";
  std::string labels = "none";      // none | edges | multigons | punctures
  std::string label_font = "small";
  bool snippet = false;
  double scale = 1.0;
  double label_offset = 0.22;
  double dot = 0.055;               // radius of a collapsed multigon
  double line_width = 0.9;          // pt, to match the paper's line art
};

//
// Command line
//

std::string fmt(const vec2& p)
{
  std::ostringstream o;
  o << std::fixed << std::setprecision(5) << "(" << p.x << "," << p.y << ")";
  return o.str();
}

void usage(std::ostream& out)
{
  out
    << "ttplot: draw a train track from its coding, as TikZ.\n"
    << "\n"
    << "Usage: ttplot [options]\n"
    << "  --coding STR          Coding of the train track\n"
    << "  --coding-file PATH    Read the coding from a file\n"
    << "  --output PATH         Output file (default: ttplot.tex)\n"
    << "  --snippet             Emit only the tikzpicture, no preamble\n"
    << "  --labels WHICH        none | edges | multigons | punctures\n"
    << "  --label-font NAME     LaTeX size, without the backslash\n"
    << "  --label-offset X      Label offset from what it labels\n"
    << "  --scale X             TikZ scale\n"
    << "  --help                Show this help\n"
    << "\n"
    << "With no --coding or --coding-file, the coding is read from stdin.\n"
    << "\n"
    << "Input coding, as printed by traintrack::print_coding:\n"
    << "  compact, e.g. \"1111 1311 2311 ...\" or \"11111 13111 ...\"\n"
    << "  hyphenated when a field reaches ten, e.g. \"1-1-12-1 1-15-2-2\"\n"
    << "\n"
    << "The multigons are drawn collapsed to points, as in Figure 12(b)\n"
    << "of the ttauto paper.  Punctures lie on a horizontal axis and the\n"
    << "track is drawn above it.  Tracks with up to seven punctures are\n"
    << "drawn without crossings; some with eight still cross.\n"
    << "\n"
    << "To make a PDF from standalone output:\n"
    << "  pdflatex ttplot.tex\n";
}

bool want(const std::string& labels, const char* which)
{
  return labels.find(which) != std::string::npos;
}

bool next_arg(const int argc, char** argv, int& i, std::string& v)
{
  if (i+1 >= argc) { std::cerr << "Missing value for " << argv[i] << "\n"; return false; }
  v = argv[++i];
  return true;
}

bool parse_args(const int argc, char** argv, options& opt)
{
  for (int i = 1; i < argc; ++i)
    {
      const std::string a(argv[i]);
      std::string v;
      if (a == "--help") { usage(std::cout); std::exit(0); }
      else if (a == "--snippet") opt.snippet = true;
      else if (a == "--coding")       { if (!next_arg(argc,argv,i,v)) return false; opt.coding = v; }
      else if (a == "--coding-file")  { if (!next_arg(argc,argv,i,v)) return false; opt.coding_file = v; }
      else if (a == "--output")       { if (!next_arg(argc,argv,i,v)) return false; opt.output_file = v; }
      else if (a == "--labels")       { if (!next_arg(argc,argv,i,v)) return false; opt.labels = v; }
      else if (a == "--label-font")   { if (!next_arg(argc,argv,i,v)) return false; opt.label_font = v; }
      else if (a == "--label-offset") { if (!next_arg(argc,argv,i,v)) return false; opt.label_offset = std::atof(v.c_str()); }
      else if (a == "--scale")        { if (!next_arg(argc,argv,i,v)) return false; opt.scale = std::atof(v.c_str()); }
      else { std::cerr << "Unknown option: " << a << "\n"; usage(std::cerr); return false; }
    }
  return true;
}

std::string read_all(std::istream& in)
{
  std::ostringstream o;
  o << in.rdbuf();
  return o.str();
}

} // namespace

int main(int argc, char** argv)
{
  options opt;
  if (!parse_args(argc,argv,opt)) return 1;

  std::string coding_text;
  if (!opt.coding_file.empty())
    {
      std::ifstream in(opt.coding_file.c_str());
      if (!in) { std::cerr << "Cannot read " << opt.coding_file << "\n"; return 1; }
      coding_text = read_all(in);
    }
  else if (!opt.coding.empty())
    coding_text = opt.coding;
  else
    coding_text = read_all(std::cin);

  // One parser, in the library: it takes each block's width from the
  // block itself rather than from the total count, so it reads both the
  // compact and the hyphenated forms and cannot mistake one width for
  // the other.
  traintrack tt(traintracks::parse_coding(coding_text));
  tt.check();

  // outer_embedding only places punctures that sit on monogons.  No
  // track the automaton builds has any other kind, but
  // build_traintrack_list(N,N2) with N2 > 0 does, so say so here rather
  // than let the library exit from further down.
  const traintrack& ctt = tt;
  for (int m = 0; m < tt.multigons(); ++m)
    if (ctt.Multigon(m).punctured() && ctt.Multigon(m).prongs() != 1)
      {
        std::cerr << "ttplot: multigon " << m << " is a punctured "
                  << ctt.Multigon(m).prongs() << "-gon.\n";
        std::cerr << "Only punctured monogons can be embedded, so this"
                  << " track cannot be drawn.\n";
        return 1;
      }

  const ttnumbering num = tt.numbering();
  const tt_embedding emb = traintracks::outer_embedding(num);
  const collapsed_layout L = traintracks::make_collapsed_layout(num,emb);

  std::ofstream fout;
  std::ostream* os = &std::cout;
  if (opt.output_file != "-")
    {
      fout.open(opt.output_file.c_str());
      if (!fout) { std::cerr << "Cannot write " << opt.output_file << "\n"; return 1; }
      os = &fout;
    }
  std::ostream& out = *os;

  if (!opt.snippet)
    out << "\\documentclass[tikz,border=4pt]{standalone}\n"
        << "\\usepackage{tikz}\n"
        << "\\begin{document}\n";

  out << "\\begin{tikzpicture}[x=1cm,y=1cm,scale=" << opt.scale
      << ",line cap=round,line join=round]\n";
  out << "  \\tikzset{ttedge/.style={line width=" << opt.line_width << "pt},"
      << "ttlab/.style={inner sep=1pt,font=\\" << opt.label_font << "}}\n";
  out << "  % Collapsed representation: every multigon is a point.\n"
      << "  % Punctures lie on the axis in boundary-walk order.\n";

  // Main edges first, so the dots sit on top of them.
  for (int e = 0; e < (int)L.edge_mg.size(); ++e)
    {
      out << "  \\draw[ttedge] " << fmt(L.arc[e].front().p0);
      for (const cubic& c : L.arc[e])
        out << " .. controls " << fmt(c.p1) << " and " << fmt(c.p2)
            << " .. " << fmt(c.p3);
      out << ";\n";
      if (want(opt.labels,"edges"))
        {
          vec2 mid = L.arc[e].front().p3;
          mid.y += opt.label_offset;
          out << "  \\node[ttlab] at " << fmt(mid) << " {$e_{" << e+1 << "}$};\n";
        }
    }

  for (int m = 0; m < (int)L.mg.size(); ++m)
    {
      out << "  \\fill[black] " << fmt(L.mg[m]) << " circle ("
          << opt.dot << ");\n";
      if (want(opt.labels,"multigons"))
        {
          vec2 p = L.mg[m];
          p.y -= opt.label_offset;
          out << "  \\node[ttlab] at " << fmt(p) << " {$m_{" << m << "}$};\n";
        }
      else if (want(opt.labels,"punctures") && L.puncture_pos[m])
        {
          vec2 p = L.mg[m];
          p.y -= opt.label_offset;
          out << "  \\node[ttlab] at " << fmt(p) << " {$"
              << L.puncture_pos[m] << "$};\n";
        }
    }

  out << "\\end{tikzpicture}\n";
  if (!opt.snippet) out << "\\end{document}\n";

  if (opt.output_file != "-")
    std::cerr << "Wrote TikZ output to " << opt.output_file << "\n";
  return 0;
}
