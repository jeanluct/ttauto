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

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include "traintracks/traintrack.hpp"

namespace {

using traintracks::edge;
using traintracks::multigon;
using traintracks::traintrack;

struct vec2
{
  double x;
  double y;
};

struct endpoint_ref
{
  int m;
  int p;
  int e;
};

struct edge_record
{
  endpoint_ref a;
  endpoint_ref b;
  int id;
};

struct cusp_record
{
  int m;
  int p;
  int e;
  int cusp_index;
  int fold_cw;
  int fold_ccw;
};

struct options
{
  std::string coding;
  std::string coding_file;
  std::string output_file = "ttplot.tex";
  std::string labels = "none";
  bool snippet = false;
  bool plot_punctures = true;
  std::string first_monogon_fill = "gray!12";
  double scale = 1.0;
  double curvature = 1.3;
  double slot_spacing = 0.18;
  double label_offset = 0.22;
  std::string label_font = "tiny";
  double fold_label_sep = 0.60;
  double fold_tangent_scale = 3.00;
  double multigon_label_inset = 0.24;
  double edge_waviness = 0.30;
  int relax_iters = 180;
  double relax_step = 0.07;
};

struct prong_pref_data
{
  std::vector<double> pref;
  std::vector<bool> has_pref;
  int k;
};

double pi() { return 3.14159265358979323846; }

vec2 add(const vec2& a, const vec2& b) { return {a.x+b.x,a.y+b.y}; }
vec2 sub(const vec2& a, const vec2& b) { return {a.x-b.x,a.y-b.y}; }
vec2 mul(const vec2& a, const double s) { return {a.x*s,a.y*s}; }
double dot2(const vec2& a, const vec2& b) { return a.x*b.x + a.y*b.y; }

vec2 dir_from_angle(const double a) { return {std::cos(a),std::sin(a)}; }

vec2 perp_ccw(const vec2& v) { return {-v.y,v.x}; }

double cubic_coord(const double p0, const double p1,
                  const double p2, const double p3, const double t)
{
  const double u = 1.0 - t;
  return u*u*u*p0 + 3.0*u*u*t*p1 + 3.0*u*t*t*p2 + t*t*t*p3;
}

vec2 cubic_point(const vec2& p0, const vec2& p1,
                 const vec2& p2, const vec2& p3, const double t)
{
  return {
    cubic_coord(p0.x,p1.x,p2.x,p3.x,t),
    cubic_coord(p0.y,p1.y,p2.y,p3.y,t)
  };
}

vec2 cubic_tangent(const vec2& p0, const vec2& p1,
                   const vec2& p2, const vec2& p3, const double t)
{
  const double u = 1.0 - t;
  const vec2 a = mul(sub(p1,p0),3.0*u*u);
  const vec2 b = mul(sub(p2,p1),6.0*u*t);
  const vec2 c = mul(sub(p3,p2),3.0*t*t);
  return add(add(a,b),c);
}

double norm(const vec2& v)
{
  return std::sqrt(v.x*v.x + v.y*v.y);
}

vec2 normalize_or_zero(const vec2& v)
{
  const double n = norm(v);
  if (n <= 1e-12) return {0.0,0.0};
  return mul(v,1.0/n);
}

vec2 unit_or_default(const vec2& v, const vec2& dv)
{
  const double n = norm(v);
  if (n <= 1e-12) return dv;
  return mul(v,1.0/n);
}

double branch_tangent_length(const options& opt, const double slot_offset)
{
  return opt.curvature + opt.edge_waviness*std::fabs(slot_offset);
}

double local_theta_from_pref(const prong_pref_data& pd, const int sign)
{
  double s = 0.0;
  double c = 0.0;
  for (int p = 0; p < pd.k; ++p)
    {
      if (!pd.has_pref[p]) continue;
      double base = sign*2.0*pi()*p/pd.k;
      double d = pd.pref[p] - base;
      s += std::sin(d);
      c += std::cos(d);
    }
  if (std::fabs(s) <= 1e-12 && std::fabs(c) <= 1e-12) return 0.0;
  return std::atan2(s,c);
}

double local_fit_error(const prong_pref_data& pd, const int sign, const double theta)
{
  double e = 0.0;
  for (int p = 0; p < pd.k; ++p)
    {
      if (!pd.has_pref[p]) continue;
      double a = theta + sign*2.0*pi()*p/pd.k;
      e += 1.0 - std::cos(pd.pref[p] - a);
    }
  return e;
}

double orient2d(const vec2& a, const vec2& b, const vec2& c)
{
  return (b.x-a.x)*(c.y-a.y) - (b.y-a.y)*(c.x-a.x);
}

bool segments_properly_intersect(const vec2& a, const vec2& b,
                                 const vec2& c, const vec2& d)
{
  const double o1 = orient2d(a,b,c);
  const double o2 = orient2d(a,b,d);
  const double o3 = orient2d(c,d,a);
  const double o4 = orient2d(c,d,b);
  const double eps = 1e-10;
  return ((o1*o2 < -eps) && (o3*o4 < -eps));
}

int crossing_score(const std::vector<edge_record>& edges,
                   const std::vector<vec2>& ctr,
                   const std::vector<std::vector<double> >& prong_angle,
                   const double prong_radius)
{
  int score = 0;
  const int E = edges.size();
  for (int i = 0; i < E; ++i)
    {
      const edge_record& ei = edges[i];
      vec2 ai = add(ctr[ei.a.m],mul(dir_from_angle(prong_angle[ei.a.m][ei.a.p]),prong_radius));
      vec2 bi = add(ctr[ei.b.m],mul(dir_from_angle(prong_angle[ei.b.m][ei.b.p]),prong_radius));

      for (int j = i+1; j < E; ++j)
        {
          const edge_record& ej = edges[j];
          if (ei.a.m == ej.a.m || ei.a.m == ej.b.m ||
              ei.b.m == ej.a.m || ei.b.m == ej.b.m)
            continue;

          vec2 aj = add(ctr[ej.a.m],mul(dir_from_angle(prong_angle[ej.a.m][ej.a.p]),prong_radius));
          vec2 bj = add(ctr[ej.b.m],mul(dir_from_angle(prong_angle[ej.b.m][ej.b.p]),prong_radius));

          if (segments_properly_intersect(ai,bi,aj,bj)) ++score;
        }
    }
  return score;
}

void relax_centers(std::vector<vec2>& ctr,
                   const std::vector<edge_record>& edges,
                   const options& opt)
{
  const int M = ctr.size();
  if (M <= 1 || opt.relax_iters <= 0 || opt.relax_step <= 0) return;

  std::vector<std::vector<double> > w(M,std::vector<double>(M,0.0));
  for (const edge_record& er : edges)
    {
      int a = er.a.m;
      int b = er.b.m;
      if (a == b) continue;
      w[a][b] += 1.0;
      w[b][a] += 1.0;
    }

  const double ideal = 2.2;
  const double k_rep = 2.2;
  const double k_spring = 0.22;
  const double k_center = 0.06;
  const double max_step = 0.18;

  for (int it = 0; it < opt.relax_iters; ++it)
    {
      std::vector<vec2> force(M,{0.0,0.0});

      for (int i = 0; i < M; ++i)
        {
          for (int j = i+1; j < M; ++j)
            {
              vec2 d = sub(ctr[j],ctr[i]);
              double r = norm(d);
              if (r < 1e-6) { d = {1e-3,0}; r = 1e-3; }
              vec2 u = mul(d,1.0/r);

              double fr = k_rep/(r*r);
              force[i] = add(force[i],mul(u,-fr));
              force[j] = add(force[j],mul(u, fr));

              if (w[i][j] > 0)
                {
                  double fs = k_spring*w[i][j]*(r-ideal);
                  force[i] = add(force[i],mul(u, fs));
                  force[j] = add(force[j],mul(u,-fs));
                }
            }
        }

      for (int i = 0; i < M; ++i)
        {
          force[i] = add(force[i],mul(ctr[i],-k_center));
          vec2 d = mul(force[i],opt.relax_step);
          double dn = norm(d);
          if (dn > max_step) d = mul(d,max_step/dn);
          ctr[i] = add(ctr[i],d);
        }
    }
}

std::vector<prong_pref_data>
compute_prong_preferences(const traintrack& tt,
                         const std::vector<edge_record>& edges,
                         const std::vector<vec2>& ctr)
{
  std::vector<prong_pref_data> pref_data(tt.multigons());
  for (int m = 0; m < tt.multigons(); ++m)
    {
      const multigon& mm = tt.Multigon(m);
      const int k = mm.prongs();
      pref_data[m].k = k;
      pref_data[m].pref.assign(k,0.0);
      pref_data[m].has_pref.assign(k,false);
      for (int p = 0; p < k; ++p)
        {
          vec2 accum = {0.0,0.0};
          for (const edge_record& er : edges)
            {
              if (er.a.m == m && er.a.p == p)
                {
                  vec2 v = sub(ctr[er.b.m],ctr[m]);
                  accum = add(accum,normalize_or_zero(v));
                }
              if (er.b.m == m && er.b.p == p)
                {
                  vec2 v = sub(ctr[er.a.m],ctr[m]);
                  accum = add(accum,normalize_or_zero(v));
                }
            }

          if (norm(accum) > 1e-12)
            {
              pref_data[m].pref[p] = std::atan2(accum.y,accum.x);
              pref_data[m].has_pref[p] = true;
            }
        }
    }
  return pref_data;
}

std::vector<std::vector<double> >
compute_prong_angles(const std::vector<prong_pref_data>& pref_data,
                     const std::vector<edge_record>& edges,
                     const std::vector<vec2>& ctr,
                     const double prong_radius)
{
  const int M = pref_data.size();
  std::vector<int> prong_sign(M,1);
  std::vector<double> prong_theta(M,0.0);
  std::vector<std::vector<double> > prong_angle(M);

  for (int m = 0; m < M; ++m)
    {
      const prong_pref_data& pd = pref_data[m];
      const double tpos = local_theta_from_pref(pd, 1);
      const double tneg = local_theta_from_pref(pd,-1);
      const double epos = local_fit_error(pd, 1, tpos);
      const double eneg = local_fit_error(pd,-1, tneg);
      if (eneg < epos)
        {
          prong_sign[m] = -1;
          prong_theta[m] = tneg;
        }
      else
        {
          prong_sign[m] = 1;
          prong_theta[m] = tpos;
        }

      prong_angle[m].assign(pd.k,0.0);
      for (int p = 0; p < pd.k; ++p)
        {
          prong_angle[m][p] = prong_theta[m] + prong_sign[m]*2.0*pi()*p/pd.k;
        }
    }

  int best_score = crossing_score(edges,ctr,prong_angle,prong_radius);
  for (int pass = 0; pass < 3; ++pass)
    {
      bool improved = false;
      for (int m = 0; m < M; ++m)
        {
          if (pref_data[m].k <= 1) continue;

          const int old_sign = prong_sign[m];
          const double old_theta = prong_theta[m];
          const std::vector<double> old_angles = prong_angle[m];

          prong_sign[m] = -old_sign;
          prong_theta[m] = local_theta_from_pref(pref_data[m],prong_sign[m]);
          for (int p = 0; p < pref_data[m].k; ++p)
            {
              prong_angle[m][p] = prong_theta[m] + prong_sign[m]*2.0*pi()*p/pref_data[m].k;
            }

          const int trial_score = crossing_score(edges,ctr,prong_angle,prong_radius);
          if (trial_score < best_score)
            {
              best_score = trial_score;
              improved = true;
            }
          else
            {
              prong_sign[m] = old_sign;
              prong_theta[m] = old_theta;
              prong_angle[m] = old_angles;
            }
        }
      if (!improved) break;
    }

  return prong_angle;
}

bool all_multigon_labels_equal(const traintrack& tt)
{
  if (tt.multigons() <= 1) return true;
  const int l0 = tt.Multigon(0).label();
  for (int m = 1; m < tt.multigons(); ++m)
    {
      if (tt.Multigon(m).label() != l0) return false;
    }
  return true;
}

std::vector<std::vector<double> >
compute_prong_side_tangents(const traintrack& tt,
                            const std::vector<edge_record>& edges,
                            const options& opt)
{
  std::vector<std::vector<double> > prong_side_tangent(tt.multigons());
  std::vector<std::vector<double> > tangent_sum(tt.multigons());
  std::vector<std::vector<int> > tangent_count(tt.multigons());

  for (int m = 0; m < tt.multigons(); ++m)
    {
      const int k = tt.Multigon(m).prongs();
      prong_side_tangent[m].assign(k,opt.curvature);
      tangent_sum[m].assign(k,0.0);
      tangent_count[m].assign(k,0);
    }

  for (const edge_record& er : edges)
    {
      const multigon& ma = tt.Multigon(er.a.m);
      const multigon& mb = tt.Multigon(er.b.m);
      const double sa = (er.a.e - 0.5*(ma.edges(er.a.p)-1))*opt.slot_spacing;
      const double sb = (er.b.e - 0.5*(mb.edges(er.b.p)-1))*opt.slot_spacing;
      const double ka = branch_tangent_length(opt,sa);
      const double kb = branch_tangent_length(opt,sb);

      tangent_sum[er.a.m][er.a.p] += ka;
      tangent_count[er.a.m][er.a.p] += 1;
      tangent_sum[er.b.m][er.b.p] += kb;
      tangent_count[er.b.m][er.b.p] += 1;
    }

  for (int m = 0; m < tt.multigons(); ++m)
    {
      for (int p = 0; p < (int)prong_side_tangent[m].size(); ++p)
        {
          if (tangent_count[m][p] == 0) continue;
          const double avgk = tangent_sum[m][p]/tangent_count[m][p];
          prong_side_tangent[m][p] = std::min(0.42,std::max(0.16,0.08 + 0.18*avgk));
        }
    }

  return prong_side_tangent;
}

std::vector<cusp_record> enumerate_cusps(const traintrack& tt)
{
  std::vector<cusp_record> out;
  int cusp = 0;
  for (int m = 0; m < tt.multigons(); ++m)
    {
      const multigon& mm = tt.Multigon(m);
      for (int p = 0; p < mm.prongs(); ++p)
        {
          for (int e = 0; e < mm.edges(p)-1; ++e)
            {
              cusp_record c;
              c.m = m;
              c.p = p;
              c.e = e;
              c.cusp_index = cusp;
              c.fold_cw = 2*cusp;
              c.fold_ccw = 2*cusp + 1;
              out.push_back(c);
              ++cusp;
            }
        }
    }
  return out;
}

std::string fmt(const vec2& p)
{
  std::ostringstream s;
  s << std::fixed << std::setprecision(5) << "(" << p.x << "," << p.y << ")";
  return s.str();
}

void usage(std::ostream& out)
{
  out
    << "Usage: ttplot [options]\n"
    << "\n"
    << "Plot a train track from coding and write TikZ output.\n"
    << "\n"
    << "Options:\n"
    << "  --coding \"...\"        Coding string (compact or integer-list form)\n"
    << "  --coding-file FILE    Read coding string from file\n"
    << "  --output FILE         Output file (default: ttplot.tex)\n"
    << "  --snippet             Emit TikZ picture only (no standalone preamble)\n"
    << "  --no-punctures        Do not draw puncture dots\n"
    << "  --first-monogon-fill COLOR  Fill first monogon with TikZ color (default: gray!12, use none to disable)\n"
    << "  --labels MODE         Label mode: none|multigons|prongs|edges|folds|all or comma list\n"
    << "  --scale X             Global scale (default: 1.0)\n"
    << "  --curvature X         Bezier control strength (default: 1.3)\n"
    << "  --slot-spacing X      Tangential fanout spacing (default: 0.18)\n"
    << "  --label-offset X      Label offset from geometry (default: 0.22)\n"
    << "  --label-font NAME     Label font: tiny|scriptsize|footnotesize|small|normalsize (default: tiny)\n"
    << "  --fold-label-sep X    Fold-label distance from cusp (default: 0.60)\n"
    << "  --fold-tangent-scale X  Tangential spread between fold labels (default: 3.00)\n"
    << "  --multigon-label-inset X  Inset of multigon labels from boundary (default: 0.24)\n"
    << "  --edge-waviness X     Extra bend from multi-edge fanout (default: 0.30)\n"
    << "  --relax-iters N       Layout relaxation iterations (default: 180)\n"
    << "  --relax-step X        Layout relaxation step (default: 0.07)\n"
    << "  --help                Show this help\n"
    << "\n"
    << "Input coding formats:\n"
    << "  1) Compact printed blocks, e.g. \"11111 13111 23111 ...\"\n"
    << "  2) Integer list, e.g. \"0 1 0 0 1 0 3 0 0 1 ...\"\n"
    << "\n"
    << "To make PDF from standalone output:\n"
    << "  pdflatex ttplot.tex\n";
}

bool parse_args(const int argc, char** argv, options& opt)
{
  for (int i = 1; i < argc; ++i)
    {
      std::string a(argv[i]);
      if (a == "--help")
        {
          usage(std::cout);
          return false;
        }
      else if (a == "--coding" && i+1 < argc)
        {
          opt.coding = argv[++i];
        }
      else if (a == "--coding-file" && i+1 < argc)
        {
          opt.coding_file = argv[++i];
        }
      else if (a == "--output" && i+1 < argc)
        {
          opt.output_file = argv[++i];
        }
      else if (a == "--labels" && i+1 < argc)
        {
          opt.labels = argv[++i];
        }
      else if (a == "--snippet")
        {
          opt.snippet = true;
        }
      else if (a == "--no-punctures")
        {
          opt.plot_punctures = false;
        }
      else if (a == "--first-monogon-fill" && i+1 < argc)
        {
          opt.first_monogon_fill = argv[++i];
          if (opt.first_monogon_fill == "none") opt.first_monogon_fill.clear();
        }
      else if (a == "--scale" && i+1 < argc)
        {
          opt.scale = std::atof(argv[++i]);
        }
      else if (a == "--curvature" && i+1 < argc)
        {
          opt.curvature = std::atof(argv[++i]);
        }
      else if (a == "--slot-spacing" && i+1 < argc)
        {
          opt.slot_spacing = std::atof(argv[++i]);
        }
      else if (a == "--label-offset" && i+1 < argc)
        {
          opt.label_offset = std::atof(argv[++i]);
        }
      else if (a == "--label-font" && i+1 < argc)
        {
          opt.label_font = argv[++i];
        }
      else if (a == "--fold-label-sep" && i+1 < argc)
        {
          opt.fold_label_sep = std::atof(argv[++i]);
        }
      else if (a == "--fold-tangent-scale" && i+1 < argc)
        {
          opt.fold_tangent_scale = std::atof(argv[++i]);
        }
      else if (a == "--multigon-label-inset" && i+1 < argc)
        {
          opt.multigon_label_inset = std::atof(argv[++i]);
        }
      else if (a == "--edge-waviness" && i+1 < argc)
        {
          opt.edge_waviness = std::atof(argv[++i]);
        }
      else if (a == "--relax-iters" && i+1 < argc)
        {
          opt.relax_iters = std::atoi(argv[++i]);
        }
      else if (a == "--relax-step" && i+1 < argc)
        {
          opt.relax_step = std::atof(argv[++i]);
        }
      else
        {
          std::cerr << "Unknown or incomplete option: " << a << "\n";
          usage(std::cerr);
          return false;
        }
    }

  return true;
}

std::string read_all_stdin()
{
  std::ostringstream s;
  s << std::cin.rdbuf();
  return s.str();
}

std::string read_file_text(const std::string& path)
{
  std::ifstream in(path.c_str());
  if (!in)
    {
      std::cerr << "Could not open coding file: " << path << "\n";
      std::exit(1);
    }
  std::ostringstream s;
  s << in.rdbuf();
  return s.str();
}

bool all_digits(const std::string& t)
{
  if (t.empty()) return false;
  for (char c : t)
    {
      if (!(c >= '0' && c <= '9')) return false;
    }
  return true;
}

std::vector<std::string> split_tokens(const std::string& text)
{
  std::string cleaned;
  cleaned.reserve(text.size());
  for (char c : text)
    {
      if (c == ',' || c == ';' || c == '[' || c == ']' || c == '{' || c == '}')
        cleaned.push_back(' ');
      else
        cleaned.push_back(c);
    }

  std::istringstream in(cleaned);
  std::vector<std::string> tok;
  std::string t;
  while (in >> t) tok.push_back(t);
  return tok;
}

traintrack::intVec parse_coding(const std::string& text)
{
  const int labeled_block_len = 5;
  const int unlabeled_block_len = 4;

  std::vector<std::string> tok(split_tokens(text));
  if (tok.empty())
    {
      std::cerr << "Empty coding input.\n";
      std::exit(1);
    }

  int compact_len = 0;
  bool compact_blocks = true;
  for (const std::string& t : tok)
    {
      if (!all_digits(t))
        {
          compact_blocks = false;
          break;
        }
      int tl = (int)t.size();
      if (tl != labeled_block_len && tl != unlabeled_block_len)
        {
          compact_blocks = false;
          break;
        }
      if (compact_len == 0) compact_len = tl;
      if (compact_len != tl)
        {
          compact_blocks = false;
          break;
        }
    }

  std::vector<int> ints;
  if (compact_blocks)
    {
      for (const std::string& t : tok)
        {
          for (char c : t) ints.push_back((int)(c - '0'));
        }
    }
  else
    {
      for (const std::string& t : tok)
        {
          if (!all_digits(t))
            {
              std::cerr << "Non-numeric token in coding: " << t << "\n";
              std::exit(1);
            }
          ints.push_back(std::atoi(t.c_str()));
        }
    }

  int block_len = 0;
  if ((int)ints.size() % labeled_block_len == 0)
    block_len = labeled_block_len;
  else if ((int)ints.size() % unlabeled_block_len == 0)
    block_len = unlabeled_block_len;
  else
    {
      std::cerr << "Coding length " << ints.size();
      std::cerr << " is not divisible by 5 (labeled) or 4 (legacy unlabeled).\n";
      std::exit(1);
    }

  traintrack::intVec code;
  code.reserve((block_len == labeled_block_len)
                 ? ints.size()
                 : (ints.size()/unlabeled_block_len)*labeled_block_len);

  bool one_based_pr_edge = true;
  if (block_len == labeled_block_len)
    {
      for (int i = 0; i < (int)ints.size(); i += labeled_block_len)
        {
          if (ints[i+0] == 0 || ints[i+3] == 0)
            {
              one_based_pr_edge = false;
              break;
            }
        }
    }
  else
    {
      for (int i = 0; i < (int)ints.size(); i += unlabeled_block_len)
        {
          if (ints[i+0] == 0 || ints[i+2] == 0)
            {
              one_based_pr_edge = false;
              break;
            }
        }
    }

  if (block_len == labeled_block_len)
    {
      for (int i = 0; i < (int)ints.size(); i += labeled_block_len)
        {
          int pr = ints[i+0] - (one_based_pr_edge ? 1 : 0);
          int np = ints[i+1];
          int lb = ints[i+2] - 1;
          int ed = ints[i+3] - (one_based_pr_edge ? 1 : 0);
          int ne = ints[i+4];
          if (pr < 0 || ed < 0)
            {
              std::cerr << "Printed coding must be 1-based in prong/edge fields.\n";
              std::exit(1);
            }
          code.push_back(pr);
          code.push_back(np);
          code.push_back(lb);
          code.push_back(ed);
          code.push_back(ne);
        }
    }
  else
    {
      for (int i = 0; i < (int)ints.size(); i += unlabeled_block_len)
        {
          int pr = ints[i+0] - (one_based_pr_edge ? 1 : 0);
          int np = ints[i+1];
          int ed = ints[i+2] - (one_based_pr_edge ? 1 : 0);
          int ne = ints[i+3];
          if (pr < 0 || ed < 0)
            {
              std::cerr << "Printed coding must be 1-based in prong/edge fields.\n";
              std::exit(1);
            }
          code.push_back(pr);
          code.push_back(np);
          code.push_back(0);
          code.push_back(ed);
          code.push_back(ne);
        }
    }

  return code;
}

const edge* edge_ptr(const multigon::edgep& ep)
{
#ifdef TRAINTRACKS_NO_SHARED_PTR
  return ep;
#else
  return ep.get();
#endif
}

std::vector<edge_record> enumerate_edges(const traintrack& tt)
{
  const traintrack& ctt = tt;
  std::map<const multigon*,int> mindex;
  for (int m = 0; m < ctt.multigons(); ++m) mindex[&ctt.Multigon(m)] = m;

  std::map<const edge*,edge_record> rec;
  int next_id = 1;
  for (int m = 0; m < ctt.multigons(); ++m)
    {
      const multigon& mm = ctt.Multigon(m);
      for (int p = 0; p < mm.prongs(); ++p)
        {
          for (int e = 0; e < mm.edges(p); ++e)
            {
              const edge* ep = edge_ptr(mm.Edge(p,e));
              auto it = rec.find(ep);
              if (it == rec.end())
                {
                  int tp = 0, te = 0;
                  multigon* tmm = ep->target_multigon(&mm,tp,te);
                  int tm = mindex[tmm];
                  edge_record er;
                  er.a = {m,p,e};
                  er.b = {tm,tp,te};
                  er.id = next_id++;
                  rec[ep] = er;
                }
            }
        }
    }

  std::vector<edge_record> out;
  for (auto it = rec.begin(); it != rec.end(); ++it) out.push_back(it->second);
  return out;
}

bool want_label(const std::string& mode, const std::string& which)
{
  if (mode == "all") return true;

  std::string tokens = mode;
  for (size_t i = 0; i < tokens.size(); ++i)
    {
      if (tokens[i] == ',' || tokens[i] == ';' || tokens[i] == '|') tokens[i] = ' ';
    }

  std::istringstream in(tokens);
  std::string tok;
  while (in >> tok)
    {
      if (tok == which) return true;
    }
  return false;
}

bool want_fold_labels(const std::string& mode)
{
  if (want_label(mode,"folds") || want_label(mode,"cusps")) return true;
  return false;
}

std::string tikz_font_command(const std::string& font)
{
  if (font.empty()) return "\\tiny";
  std::string f = font;
  if (f[0] == '\\') f = f.substr(1);

  if (f == "tiny") return "\\tiny";
  if (f == "scriptsize") return "\\scriptsize";
  if (f == "footnotesize") return "\\footnotesize";
  if (f == "small") return "\\small";
  if (f == "normalsize") return "\\normalsize";
  return "\\tiny";
}

vec2 multigon_label_position(const std::vector<vec2>& ctr,
                             const std::vector<std::vector<double> >& prong_angle,
                             const int m,
                             const double prong_radius,
                             const double monogon_depth,
                             const double multigon_label_inset,
                             const double label_offset,
                             const int k)
{
  (void)prong_radius;
  (void)multigon_label_inset;
  (void)label_offset;
  if (k <= 1)
    {
      const vec2 d = dir_from_angle(prong_angle[m][0]);
      const vec2 punct = add(ctr[m],mul(d,-0.62*monogon_depth));
      return add(punct,mul(d,-0.06*monogon_depth));
    }
  return ctr[m];
}

} // namespace

int main(int argc, char** argv)
{
  options opt;
  if (!parse_args(argc,argv,opt)) return 0;

  if (opt.scale <= 0 || opt.curvature <= 0 || opt.slot_spacing <= 0 ||
      opt.label_offset < 0 || opt.edge_waviness < 0 ||
      opt.fold_label_sep < 0 || opt.fold_tangent_scale <= 0 ||
      opt.multigon_label_inset < 0 ||
      opt.relax_iters < 0 || opt.relax_step <= 0)
    {
      std::cerr << "Scale/curvature/slot-spacing/relax-step must be positive;"
                << " label offset, edge-waviness, fold-label-sep, multigon-label-inset,"
                << " and relax-iters must be nonnegative; fold-tangent-scale must be positive.\n";
      return 1;
    }

  std::string coding_text;
  if (!opt.coding_file.empty())
    coding_text = read_file_text(opt.coding_file);
  else if (!opt.coding.empty())
    coding_text = opt.coding;
  else
    coding_text = read_all_stdin();

  traintrack::intVec code(parse_coding(coding_text));
  traintrack tt(code);
  tt.check();

  std::vector<edge_record> edges = enumerate_edges(tt);

  const traintrack& ctt = tt;
  const int M = ctt.multigons();
  int first_monogon = -1;
  for (int m = 0; m < M; ++m)
    {
      if (ctt.Multigon(m).prongs() == 1)
        {
          first_monogon = m;
          break;
        }
    }
  const double center_radius = std::max(2.5, 1.0 + 0.9*M);
  const double prong_radius = 0.45;
  const double monogon_depth = 0.90;
  const double monogon_tail_tangent = 0.55;
  const double trigon_side_tangent_scale = 0.78;
  const double puncture_radius = 0.040;
  const double prong_label_base_offset = 0.04;

  std::vector<vec2> ctr(M);
  for (int m = 0; m < M; ++m)
    {
      const double a = 2.0*pi()*m/M;
      ctr[m] = mul(dir_from_angle(a),center_radius);
    }

  relax_centers(ctr,edges,opt);

  const std::vector<prong_pref_data> pref_data = compute_prong_preferences(ctt,edges,ctr);
  const std::vector<std::vector<double> > prong_angle =
    compute_prong_angles(pref_data,edges,ctr,prong_radius);
  const bool same_multigon_label = all_multigon_labels_equal(ctt);
  const std::vector<std::vector<double> > prong_side_tangent =
    compute_prong_side_tangents(ctt,edges,opt);
  const std::vector<cusp_record> cusps = enumerate_cusps(ctt);

  std::ofstream out(opt.output_file.c_str());
  if (!out)
    {
      std::cerr << "Could not open output file: " << opt.output_file << "\n";
      return 1;
    }

  if (!opt.snippet)
    {
      out << "\\documentclass[tikz,border=4pt]{standalone}\n";
      out << "\\usepackage{tikz}\n";
      out << "\\begin{document}\n";
    }

  out << "\\begin{tikzpicture}[x=1cm,y=1cm,scale=" << opt.scale << ","
      << " line cap=round,line join=round,font=" << tikz_font_command(opt.label_font) << "]\n";

  out << "  \\tikzset{ttedge/.style={line width=0.9pt},"
      << " ttcore/.style={line width=0.7pt,draw=black!70},"
      << " ttlab/.style={inner sep=0.7pt}}\n";

  out << "  % Canonical prong tangency convention:\n";
  out << "  %   edges: +d direction at prong\n";
  out << "  %   side convention: fixed cw\n";

  constexpr int side_sign = -1;

  for (int m = 0; m < M; ++m)
    {
      const multigon& mm = ctt.Multigon(m);
      const int k = mm.prongs();

      if (k == 1)
        {
          const double a = prong_angle[m][0];
          vec2 d = dir_from_angle(a);
          vec2 t = perp_ccw(d);
          vec2 tip = add(ctr[m],mul(d,prong_radius));
          vec2 tail = add(ctr[m],mul(d,-monogon_depth));
          if (!opt.first_monogon_fill.empty() && m == first_monogon)
            {
              const double st = prong_side_tangent[m][0];
              vec2 c1f = add(tip,mul(d,side_sign*st));
              vec2 c2f = add(tail,mul(t,monogon_tail_tangent));
              vec2 c3f = add(tail,mul(t,-monogon_tail_tangent));
              vec2 c4f = add(tip,mul(d,side_sign*st));
              out << "  \\filldraw[fill=" << opt.first_monogon_fill
                  << ",draw=none] " << fmt(tip)
                  << " .. controls " << fmt(c1f) << " and " << fmt(c2f)
                  << " .. " << fmt(tail)
                  << " .. controls " << fmt(c3f) << " and " << fmt(c4f)
                  << " .. " << fmt(tip) << ";\n";
            }
          // Side arcs should approach the prong from inside the multigon,
          // while track edges leave outward along +d.
          const double st = prong_side_tangent[m][0];
          vec2 c1 = add(tip,mul(d,side_sign*st));
          vec2 c2 = add(tail,mul(t,monogon_tail_tangent));
          vec2 c3 = add(tail,mul(t,-monogon_tail_tangent));
          vec2 c4 = add(tip,mul(d,side_sign*st));
          out << "  \\draw[ttcore] " << fmt(tip)
              << " .. controls " << fmt(c1) << " and " << fmt(c2)
              << " .. " << fmt(tail)
              << " .. controls " << fmt(c3) << " and " << fmt(c4)
              << " .. " << fmt(tip) << ";\n";
        }
      else
        {
          for (int p = 0; p < k; ++p)
            {
              const int q = (p+1) % k;
              vec2 dp = dir_from_angle(prong_angle[m][p]);
              vec2 dq = dir_from_angle(prong_angle[m][q]);
              vec2 pp = add(ctr[m],mul(dp,prong_radius));
              vec2 qq = add(ctr[m],mul(dq,prong_radius));
              // Canonical tangent frame at each prong:
              //  - branch edges leave along +d
              //  - side arcs approach with fixed CW orientation
              const double stp = prong_side_tangent[m][p];
              const double stq = prong_side_tangent[m][q];
              const double side_scale = (k == 3 ? trigon_side_tangent_scale : 1.0);
              vec2 c1 = add(pp,mul(dp,side_sign*side_scale*stp));
              vec2 c2 = add(qq,mul(dq,side_sign*side_scale*stq));
              out << "  \\draw[ttcore] " << fmt(pp)
                  << " .. controls " << fmt(c1)
                  << " and " << fmt(c2)
                  << " .. " << fmt(qq) << ";\n";
            }
        }

      for (int p = 0; p < k; ++p)
        {
          const double a = prong_angle[m][p];
          vec2 d = dir_from_angle(a);

          if (want_label(opt.labels,"prongs"))
            {
              vec2 lp = add(ctr[m],mul(d,prong_radius+prong_label_base_offset));
              lp = add(lp,mul(d,0.35*opt.label_offset));
              out << "  \\node[ttlab] at " << fmt(lp) << " {$p_" << (p+1) << "$};\n";
            }
        }

      if (want_label(opt.labels,"multigons"))
        {
          vec2 ml = multigon_label_position(ctr,prong_angle,m,prong_radius,
                                            monogon_depth,opt.multigon_label_inset,
                                            opt.label_offset,k);
          out << "  \\node[ttlab] at " << fmt(ml) << " {$m_" << (m+1);
          if (!same_multigon_label) out << "\\!:\\!" << (mm.label()+1);
          out << "$};\n";
        }

      if (opt.plot_punctures && mm.punctured())
        {
          // In core traintracks code, punctures live on monogons and bigons
          // by default (k == 1 or k == 2), unless explicitly changed.
          if (k == 1)
            {
              const double a = prong_angle[m][0];
              vec2 d = dir_from_angle(a);
              vec2 punct = add(ctr[m],mul(d,-0.62*monogon_depth));
              out << "  \\fill[black] " << fmt(punct) << " circle (" << puncture_radius << ");\n";
            }
          else
            {
              out << "  \\fill[black] " << fmt(ctr[m]) << " circle (" << puncture_radius << ");\n";
            }
        }
    }

  for (const edge_record& er : edges)
    {
      const multigon& ma = ctt.Multigon(er.a.m);
      const multigon& mb = ctt.Multigon(er.b.m);

      const double aa = prong_angle[er.a.m][er.a.p];
      const double ab = prong_angle[er.b.m][er.b.p];

      vec2 da = dir_from_angle(aa);
      vec2 db = dir_from_angle(ab);
      const double sa = (er.a.e - 0.5*(ma.edges(er.a.p)-1))*opt.slot_spacing;
      const double sb = (er.b.e - 0.5*(mb.edges(er.b.p)-1))*opt.slot_spacing;

      vec2 pa = add(ctr[er.a.m],mul(da,prong_radius));
      vec2 pb = add(ctr[er.b.m],mul(db,prong_radius));

      const double ka = branch_tangent_length(opt,sa);
      const double kb = branch_tangent_length(opt,sb);
      vec2 c1 = add(pa,mul(da,ka));
      vec2 c2 = add(pb,mul(db,kb));

      vec2 ua = unit_or_default(sub(pa,ctr[er.a.m]),da);
      vec2 ub = unit_or_default(sub(pb,ctr[er.b.m]),db);
      const double ra0 = norm(sub(pa,ctr[er.a.m]));
      const double rb0 = norm(sub(pb,ctr[er.b.m]));
      const double ra = dot2(sub(c1,ctr[er.a.m]),ua);
      const double rb = dot2(sub(c2,ctr[er.b.m]),ub);
      const double min_out = 0.14;
      if (ra < ra0 + min_out) c1 = add(c1,mul(ua,(ra0 + min_out - ra)));
      if (rb < rb0 + min_out) c2 = add(c2,mul(ub,(rb0 + min_out - rb)));

      out << "  \\draw[ttedge] " << fmt(pa) << " .. controls " << fmt(c1)
          << " and " << fmt(c2) << " .. " << fmt(pb) << ";\n";

      if (want_label(opt.labels,"edges"))
        {
          vec2 mid = cubic_point(pa,c1,c2,pb,0.5);
          vec2 tan = cubic_tangent(pa,c1,c2,pb,0.5);
          vec2 nrm = unit_or_default(perp_ccw(tan),{0.0,1.0});
          vec2 lbl = add(mid,mul(nrm,opt.label_offset));
          out << "  \\node[ttlab] at " << fmt(lbl) << " {$e_" << er.id << "$};\n";
        }
    }

  if (want_fold_labels(opt.labels))
    {
      for (const cusp_record& c : cusps)
        {
          const double a = prong_angle[c.m][c.p];
          vec2 d = dir_from_angle(a);
          vec2 t = perp_ccw(d);

          const multigon& mm = ctt.Multigon(c.m);
          const double slotL = (c.e - 0.5*(mm.edges(c.p)-1))*opt.slot_spacing;
          const double slotR = ((c.e+1) - 0.5*(mm.edges(c.p)-1))*opt.slot_spacing;
          const double slotM = 0.5*(slotL + slotR);

          vec2 prong_point = add(ctr[c.m],mul(d,prong_radius));
          vec2 cusp_mid = add(prong_point,mul(t,opt.fold_tangent_scale*slotM));
          vec2 lbl = add(cusp_mid,mul(d,opt.fold_label_sep + opt.label_offset));

          out << "  \\node[ttlab] at " << fmt(lbl)
              << " {" << c.fold_cw << "/" << c.fold_ccw << "};\n";
        }
    }

  out << "\\end{tikzpicture}\n";
  if (!opt.snippet)
    {
      out << "\\end{document}\n";
    }

  out.close();

  std::cout << "Wrote TikZ output to " << opt.output_file << "\n";
  if (!opt.snippet)
    {
      std::cout << "Compile with: pdflatex " << opt.output_file << "\n";
    }

  return 0;
}
