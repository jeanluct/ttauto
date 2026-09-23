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

// Braid words, and the braid that restores proper embedding after a fold.
// See include/traintracks/braid.hpp.

#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <set>
#include <vector>

#include "traintracks/braid.hpp"
#include "traintracks/util.hpp"

namespace traintracks {

int braidword::exponent_sum() const
{
  int s = 0;
  for (int k = 0; k < (int)w.size(); ++k) s += (w[k] > 0 ? 1 : -1);
  return s;
}

std::vector<int> braidword::permutation() const
{
  // slot[k] holds the puncture that started at position slot[k].
  std::vector<int> slot(n);
  for (int k = 0; k < n; ++k) slot[k] = k+1;
  for (int k = 0; k < (int)w.size(); ++k)
    {
      const int g = std::abs(w[k]);
      std::swap(slot[g-1],slot[g]);
    }
  std::vector<int> perm(n);
  for (int k = 0; k < n; ++k) perm[slot[k]-1] = k+1;
  return perm;
}

braidword& braidword::operator*=(const braidword& b)
{
  if (n == 0) n = b.n;
  if (b.n != 0 && b.n != n)
    {
      std::cerr << "Braids on " << n << " and " << b.n
                << " strings in traintracks::braidword::operator*=.\n";
      std::exit(1);
    }
  w.insert(w.end(),b.w.begin(),b.w.end());
  return *this;
}

braidword braidword::inverse() const
{
  std::vector<int> v;
  v.reserve(w.size());
  for (int k = (int)w.size()-1; k >= 0; --k) v.push_back(-w[k]);
  return braidword(n,v);
}

braidword& braidword::reduce()
{
  std::vector<int> v;
  v.reserve(w.size());
  for (int k = 0; k < (int)w.size(); ++k)
    {
      if (!v.empty() && v.back() == -w[k]) v.pop_back(); else v.push_back(w[k]);
    }
  w.swap(v);
  return *this;
}

namespace {

// One generator acting on Dynnikov coordinates, as in Dynnikov's update
// rules (the same ones braidlab uses).
inline double dpos(const double x) { return (x > 0 ? x : 0); }
inline double dneg(const double x) { return (x < 0 ? x : 0); }

void dynnikov_step(const int n, const std::vector<int>& w,
                   std::vector<double>& a, std::vector<double>& b)
{
  std::vector<double> at(a), bt(b);
  for (std::size_t g = 0; g < w.size(); ++g)
    {
      const int i = std::abs(w[g]);
      if (w[g] > 0)
        {
          if (i == 1)
            { bt[1] = a[1] + dpos(b[1]); at[1] = -b[1] + dpos(bt[1]); }
          else if (i == n-1)
            { bt[n-2] = a[n-2] + dneg(b[n-2]);
              at[n-2] = -b[n-2] + dneg(bt[n-2]); }
          else
            { const double c = a[i-1] - a[i] - dpos(b[i]) + dneg(b[i-1]);
              at[i-1] = a[i-1] - dpos(b[i-1]) - dpos(dpos(b[i])+c);
              bt[i-1] = b[i] + dneg(c);
              at[i] = a[i] - dneg(b[i]) - dneg(dneg(b[i-1])-c);
              bt[i] = b[i-1] - dneg(c); }
        }
      else
        {
          if (i == 1)
            { bt[1] = -a[1] + dpos(b[1]); at[1] = b[1] - dpos(bt[1]); }
          else if (i == n-1)
            { bt[n-2] = -a[n-2] + dneg(b[n-2]);
              at[n-2] = b[n-2] - dneg(bt[n-2]); }
          else
            { const double d = a[i-1] - a[i] + dpos(b[i]) - dneg(b[i-1]);
              at[i-1] = a[i-1] + dpos(b[i-1]) + dpos(dpos(b[i])-d);
              bt[i-1] = b[i] - dpos(d);
              at[i] = a[i] + dneg(b[i]) + dneg(dneg(b[i-1])+d);
              bt[i] = b[i-1] + dpos(d); }
        }
      for (int k = 1; k <= n-2; ++k) { a[k] = at[k]; b[k] = bt[k]; }
    }
}

} // namespace

double braidword::growth() const
{
  if (n < 3 || w.empty()) return 1;

  double best = 0;
  // A few starting loops, since a particular one can miss the largest
  // piece of a reducible braid.
  for (int trial = 0; trial < 4; ++trial)
    {
      std::vector<double> a(n-1,0.0), b(n-1,0.0);
      for (int k = 1; k <= n-2; ++k)
        { a[k] = (trial == 0 ? 0 : (k % 3) - 1); b[k] = -1 - (trial % 2)*k; }

      double lam = 1, prev = 0;
      for (int it = 0; it < 400; ++it)
        {
          dynnikov_step(n,w,a,b);
          double nrm = 0;
          for (int k = 1; k <= n-2; ++k) nrm += a[k]*a[k] + b[k]*b[k];
          nrm = std::sqrt(nrm);
          if (!(nrm > 1e-300) || !std::isfinite(nrm)) break;
          for (int k = 1; k <= n-2; ++k) { a[k] /= nrm; b[k] /= nrm; }
          if (it > 200)
            {
              lam = nrm;
              if (std::fabs(lam-prev) < 1e-12*std::max(1.0,lam)) break;
              prev = lam;
            }
        }
      if (lam > best) best = lam;
    }
  return best;
}

braidword braidword::block_swap(const int n, const int i, const int a,
                                const int b, const int sign)
{
  if (a < 0 || b < 0 || i < 1 || i+a+b-1 > n)
    {
      std::cerr << "Block [" << i << "," << a << "," << b << "] does not fit "
                << n << " punctures in traintracks::braidword::block_swap.\n";
      std::exit(1);
    }
  std::vector<int> v;
  v.reserve(a*b);
  // Slide the left block right one puncture at a time, starting with its
  // rightmost.  With a == 1 this is sigma_i sigma_{i+1} ... sigma_{i+b-1}.
  for (int p = 1; p <= a; ++p)
    for (int k = 0; k < b; ++k)
      v.push_back(sign > 0 ? (i+a-p+k) : -(i+a-p+k));
  return braidword(n,v);
}

braidword braidword::delta(const int n, const int sign)
{
  return block_swap(n,1,1,n-1,sign);
}

std::ostream& braidword::print(std::ostream& strm) const
{
  for (int k = 0; k < (int)w.size(); ++k)
    {
      if (k) strm << " ";
      strm << w[k];
    }
  return strm;
}

std::ostream& braidword::printMathematicaForm(std::ostream& strm) const
{
  strm << "{";
  for (int k = 0; k < (int)w.size(); ++k)
    {
      if (k) strm << ",";
      strm << w[k];
    }
  return strm << "}";
}

namespace {

// Punctures of the subtree that hangs on the far end of edge `moved`, away
// from the multigon whose cusp is being folded.  The multigons and main
// edges form a tree, so this is one component of the tree cut at `moved`.
std::set<int> moved_subtree(const ttnumbering& num, const int moved,
                            const int cusp_prong)
{
  const int np = num.nprongs();
  int nmg = 0;
  for (int q = 0; q < np; ++q) nmg = std::max(nmg,num.prong[q].multigon+1);

  const int mg_cusp = num.prong[cusp_prong].multigon;
  const int mg_tail = num.prong[num.edge_tail[moved]].multigon;
  const int mg_head = num.prong[num.edge_head[moved]].multigon;

  std::vector<bool> seen(nmg,false);
  std::vector<int> todo(1,(mg_tail == mg_cusp ? mg_head : mg_tail));
  seen[todo.back()] = true;
  while (!todo.empty())
    {
      const int m = todo.back();
      todo.pop_back();
      for (int e = 0; e < num.nedges(); ++e)
        {
          if (e == moved) continue;
          const int a = num.prong[num.edge_tail[e]].multigon;
          const int b = num.prong[num.edge_head[e]].multigon;
          int o = -1;
          if (a == m) o = b; else if (b == m) o = a; else continue;
          if (!seen[o]) { seen[o] = true; todo.push_back(o); }
        }
    }

  std::set<int> s;
  for (int q = 0; q < np; ++q)
    if (num.prong[q].punctured && seen[num.prong[q].multigon]) s.insert(q);
  return s;
}

// Prong of the target multigon, the one at the far end of `onto`.
int target_prong(const ttnumbering& num, const int onto, const int cusp_prong)
{
  const int mg_cusp = num.prong[cusp_prong].multigon;
  const int pt = num.edge_tail[onto], ph = num.edge_head[onto];
  return (num.prong[pt].multigon == mg_cusp ? ph : pt);
}

} // namespace

fold_swap fold_block_swap(const ttnumbering& before_num,
                          const fold_map_data& fm,
                          const tt_embedding& before)
{
  fold_swap fs;
  fs.dir = fm.dir;

  const int np = before.npunctures();

  // A fold onto an unpunctured multigon leaves the track properly
  // embedded, so nothing moves.
  if (!fm.after.prong[fm.target_from].punctured)
    {
      fs.next_cut_dart = fm.after.side_letter(
                           fm.prong_image[before.puncture_order[0]]);
      return fs;
    }

  const std::set<int> S = moved_subtree(before_num,fm.moved,fm.cusp_prong);
  const int T = target_prong(before_num,fm.onto,fm.cusp_prong);
  const int a = (int)S.size();

  // S and T are consecutive around the disc, S first when the fold runs
  // clockwise and T first when it runs anticlockwise.
  const int t0 = before.position_of[T] - 1;
  const int start = (fm.dir > 0 ? traintracks::mod(t0-a,np) : t0);
  const int len = a + 1;

  std::set<int> block, want(S);
  want.insert(T);
  for (int i = 0; i < len; ++i)
    block.insert(before.puncture_order[(start+i) % np]);
  if (block != want)
    {
      std::cerr << "The punctures a fold moves are not consecutive"
                << " in traintracks::fold_block_swap.\n";
      std::exit(1);
    }

  // Move the cut off the block if it falls inside it.
  fs.rotate = std::max(0,start+len-np);
  const tt_embedding& rot =
    (fs.rotate == 0 ? before
     : outer_embedding(before_num,
         before_num.side_letter(before.puncture_order[fs.rotate])));

  const int t = rot.position_of[T];
  fs.first = (fm.dir > 0 ? t-a : t);
  fs.nleft = (fm.dir > 0 ? a : 1);
  fs.nright = (fm.dir > 0 ? 1 : a);
  fs.trivial = false;

  // Where the swap leaves the punctures, and so where to cut next.
  std::vector<int> pos(np);
  for (int i = 0; i < np; ++i) pos[i] = i+1;
  for (int x = fs.first; x < fs.first+fs.nleft; ++x) pos[x-1] = x + fs.nright;
  for (int x = fs.first+fs.nleft; x < fs.first+fs.nleft+fs.nright; ++x)
    pos[x-1] = x - fs.nleft;
  int at1 = -1;
  for (int i = 0; i < np; ++i) if (pos[i] == 1) at1 = i;
  fs.next_cut_dart =
    fm.after.side_letter(fm.prong_image[rot.puncture_order[at1]]);

  return fs;
}

namespace {

// Handedness of the block swap, and of the rotation that moves the cut.
// Which way the punctures go is already fixed by which of S and T comes
// first; this is the remaining bit, whether the moving block passes in
// front of or behind the other, and it follows the direction the fold
// carried the moved edge around the puncture.  Calibrated against Toby
// Hall's Trains: see testsuite/ttauto/test_braid_extraction.cpp.  The
// mirror convention, sign +fm.dir with rotation_sign -1, gives the mirror
// image of every braid, which is equally correct for the mirror image of
// every track.
inline int crossing_sign(const int dir) { return -dir; }
const int rotation_sign = 1;

} // namespace

braidword rotation_braid(const int n, const int k)
{
  braidword b(n);
  const braidword d = braidword::delta(n,rotation_sign);
  for (int i = 0; i < k; ++i) b *= d;
  return b;
}

braidword fold_braid(const ttnumbering& before_num, const fold_map_data& fm,
                     const tt_embedding& before)
{
  const fold_swap fs = fold_block_swap(before_num,fm,before);
  const int np = before.npunctures();
  braidword b = rotation_braid(np,fs.rotate);
  if (!fs.trivial)
    b *= braidword::block_swap(np,fs.first,fs.nleft,fs.nright,
                               crossing_sign(fs.dir));
  return b;
}

} // namespace traintracks
