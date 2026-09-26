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
// This is the "collapsed" representation of the ttauto paper, Fig. 11(b):
// every multigon is shrunk to a point and only the main edges are drawn,
// which the paper calls "a simplified, simply-connected version" of the
// track.  Style follows the paper's figures: plain black line art, no
// outer boundary, sparse italic labels.
//
// The layout comes from traintracks::outer_embedding, so the picture is a
// proper embedding: punctures lie on an invisible horizontal axis in the
// order the boundary walk meets them, and the track sits in the upper
// half plane.  That makes it planar by construction rather than by
// searching for a layout with few crossings.  Checked over every vertex
// of every stratum for n = 3..7, 3700 tracks, with no crossing.
//
// The full representation of Fig. 11(a), with the monogon loops and the
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
#include "traintracks/embedding.hpp"
#include "traintracks/traintrack.hpp"

namespace {

using traintracks::multigon;
using traintracks::traintrack;
using traintracks::ttnumbering;
using traintracks::tt_embedding;

struct vec2 { double x; double y; };

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
// Layout
//

struct layout
{
  std::vector<vec2> mg;                       // by multigon index
  std::vector<std::pair<int,int> > edge_mg;   // by edge: its two multigons
  std::vector<std::pair<int,int> > edge_pr;   // by edge: its two prongs
  std::vector<vec2> prong_dir;                // by prong: outgoing tangent
  std::vector<int> puncture_pos;              // by multigon, 1..n, 0 if none
};

int multigon_of_prong(const ttnumbering& num, const int q)
{
  return num.prong[q].multigon;
}

// Place the multigons from the embedding.
//
// Punctures sit on the axis at the position the boundary walk gives them.
// An unpunctured multigon sits over the middle of the puncture positions
// in its subtree, at a height that grows with the width of that span.
// The walk order is the planar order, so those spans nest instead of
// interleaving, and a wider structure therefore sits above whatever it
// contains.
layout build_layout(const traintrack& tt, const tt_embedding& emb,
                    const ttnumbering& num)
{
  const int M = tt.multigons(), E = num.nedges();

  layout L;
  L.mg.assign(M,{0.0,0.0});
  L.puncture_pos.assign(M,0);
  L.edge_mg.assign(E,std::make_pair(-1,-1));
  L.edge_pr.assign(E,std::make_pair(-1,-1));
  for (int e = 0; e < E; ++e)
    {
      L.edge_mg[e] = std::make_pair(multigon_of_prong(num,num.edge_tail[e]),
                                    multigon_of_prong(num,num.edge_head[e]));
      L.edge_pr[e] = std::make_pair(num.edge_tail[e],num.edge_head[e]);
    }

  for (int i = 0; i < emb.npunctures(); ++i)
    L.puncture_pos[multigon_of_prong(num,emb.puncture_order[i])] = i+1;

  std::vector<std::vector<int> > adj(M);
  for (int e = 0; e < E; ++e)
    {
      adj[L.edge_mg[e].first].push_back(L.edge_mg[e].second);
      adj[L.edge_mg[e].second].push_back(L.edge_mg[e].first);
    }

  // The multigons and main edges form a tree, so rooting it and taking
  // the span of each subtree needs no search.
  const int root = multigon_of_prong(num,emb.puncture_order[0]);
  std::vector<int> lo(M,M+1), hi(M,0), parent(M,-2), order;
  {
    std::vector<int> stack(1,root);
    parent[root] = -1;
    while (!stack.empty())
      {
        const int m = stack.back(); stack.pop_back();
        order.push_back(m);
        for (int k = 0; k < (int)adj[m].size(); ++k)
          if (parent[adj[m][k]] == -2)
            { parent[adj[m][k]] = m; stack.push_back(adj[m][k]); }
      }
  }
  for (int i = (int)order.size()-1; i >= 0; --i)
    {
      const int m = order[i];
      if (L.puncture_pos[m])
        {
          lo[m] = std::min(lo[m],L.puncture_pos[m]);
          hi[m] = std::max(hi[m],L.puncture_pos[m]);
        }
      if (parent[m] >= 0)
        {
          lo[parent[m]] = std::min(lo[parent[m]],lo[m]);
          hi[parent[m]] = std::max(hi[parent[m]],hi[m]);
        }
    }

  for (int m = 0; m < M; ++m)
    {
      if (L.puncture_pos[m])
        L.mg[m] = { (double)L.puncture_pos[m], 0.0 };
      else
        L.mg[m] = { 0.5*(lo[m]+hi[m]), 0.55*(hi[m]-lo[m]) };
    }

  // A multigon's prongs are equally spaced by angle around it, and the
  // edges are then attached to the prongs, with every edge at a prong
  // leaving along that prong's direction.  Edges at the same prong are
  // therefore tangent.
  //
  // Only the overall rotation of the star is free.  Pointing each prong
  // at its own neighbours instead, which is what a naive fit does, lets
  // the prongs bunch together and the multigon stops looking like a
  // k-gon at all.
  L.prong_dir.assign(num.nprongs(),{0.0,1.0});

  // Where each prong would like to point, from the edges attached to it.
  std::vector<vec2> want(num.nprongs(),{0.0,0.0});
  for (int e = 0; e < E; ++e)
    {
      const int qt = L.edge_pr[e].first, qh = L.edge_pr[e].second;
      const vec2 a = L.mg[L.edge_mg[e].first], b = L.mg[L.edge_mg[e].second];
      want[qt].x += b.x-a.x; want[qt].y += b.y-a.y;
      want[qh].x += a.x-b.x; want[qh].y += a.y-b.y;
    }

  // The sense in which the prong index runs round a multigon in the
  // plane.  One convention for the whole picture, never per multigon.
  const double sense = 1.0;

  for (int m = 0; m < M; ++m)
    {
      // A puncture sits on the axis and the track lives above it, so its
      // one prong points straight up and the loop hangs below.
      if (L.puncture_pos[m]) continue;

      // Circular mean of what the prongs want, each shifted back by the
      // angle the equal spacing will give it.
      double sx = 0.0, sy = 0.0;
      int k = 0;
      for (int q = 0; q < num.nprongs(); ++q)
        {
          if (num.prong[q].multigon != m) continue;
          k = num.prong[q].nprongs;
          const double n2 = std::sqrt(want[q].x*want[q].x + want[q].y*want[q].y);
          if (n2 < 1e-9) continue;
          const double a = std::atan2(want[q].y,want[q].x)
                         - sense*2.0*M_PI*num.prong[q].prong/k;
          sx += std::cos(a); sy += std::sin(a);
        }
      if (k == 0) continue;
      const double theta = (sx*sx + sy*sy > 1e-18) ? std::atan2(sy,sx) : 0.5*M_PI;

      for (int q = 0; q < num.nprongs(); ++q)
        {
          if (num.prong[q].multigon != m) continue;
          const double a = theta + sense*2.0*M_PI*num.prong[q].prong/k;
          L.prong_dir[q] = { std::cos(a), std::sin(a) };
        }
    }

  return L;
}

// The control point at each end runs along that end's prong direction,
// so edges sharing a prong are tangent there.  Its length grows with the
// horizontal span, so an edge reaching over intervening structure still
// arcs above it rather than cutting through.
void edge_controls(const vec2& a, const vec2& b,
                   const vec2& da, const vec2& db, vec2& c1, vec2& c2)
{
  const double h = 0.55*std::fabs(a.x-b.x) + 0.30;
  c1 = { a.x + h*da.x, a.y + h*da.y };
  c2 = { b.x + h*db.x, b.y + h*db.y };
}

vec2 cubic_point(const vec2& p0, const vec2& p1, const vec2& p2,
                 const vec2& p3, const double t)
{
  const double u = 1.0-t;
  const double w0 = u*u*u, w1 = 3*u*u*t, w2 = 3*u*t*t, w3 = t*t*t;
  return { w0*p0.x + w1*p1.x + w2*p2.x + w3*p3.x,
           w0*p0.y + w1*p1.y + w2*p2.y + w3*p3.y };
}

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
    << "The multigons are drawn collapsed to points, as in Figure 11(b)\n"
    << "of the ttauto paper.  Punctures lie on a horizontal axis and the\n"
    << "track is drawn above it, so the picture has no crossings.\n"
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
  const layout L = build_layout(tt,emb,num);

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
      const vec2 a = L.mg[L.edge_mg[e].first], b = L.mg[L.edge_mg[e].second];
      vec2 c1, c2;
      edge_controls(a,b,L.prong_dir[L.edge_pr[e].first],
                    L.prong_dir[L.edge_pr[e].second],c1,c2);
      out << "  \\draw[ttedge] " << fmt(a) << " .. controls " << fmt(c1)
          << " and " << fmt(c2) << " .. " << fmt(b) << ";\n";
      if (want(opt.labels,"edges"))
        {
          vec2 mid = cubic_point(a,c1,c2,b,0.5);
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
