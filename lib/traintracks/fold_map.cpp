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

// traintrack::fold_with_map and fold_map_data: one fold as a train-track
// map in canonical numbering.  See include/traintracks/fold_map.hpp for the
// conventions.

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <vector>

#include "traintracks/fold_map.hpp"
#include "traintracks/traintrack.hpp"
#include "traintracks/util.hpp"

namespace traintracks {

namespace {

// Geometric identity of a prong that survives normalise(): the type of its
// multigon, the set of edge objects attached at the prong, and the set of
// edge objects attached to the whole multigon.  Edge objects persist through
// fold() and normalise(); multigon objects do not (swap exchanges their
// contents), and edge ending indices do not either.
struct prong_signature
{
  int nprongs;
  bool punctured;
  std::vector<const edge*> at_prong;
  std::vector<const edge*> at_multigon;

  bool operator<(const prong_signature& o) const
  {
    if (nprongs != o.nprongs) return nprongs < o.nprongs;
    if (punctured != o.punctured) return punctured < o.punctured;
    if (at_prong != o.at_prong) return at_prong < o.at_prong;
    return at_multigon < o.at_multigon;
  }
};

void fail(const char* msg)
{
  std::cerr << msg << " in traintracks::traintrack::fold_with_map.\n";
  std::exit(1);
}

} // anonymous namespace


bool traintrack::fold_with_map(const int f, fold_map_data& fm)
{
  require_normalised("traintrack::fold_with_map");

  if (f < 0 || f >= foldings()) fail("Illegal folding index");

  fm = fold_map_data();
  fm.before = numbering();
  const ttnumbering& B = fm.before;
  fm.nmain = B.nedges();
  fm.nsides = B.nprongs();

  // Locate the cusp and the two edges from the numbering, as fold(int)
  // does.
  int mc0 = -1, pc = -1, ec = -1;
  B.fold_cusp(f,mc0,pc,ec);
  multigon* mmc = &Multigon(mc0);
  const int dir = 1 - 2*(f % 2);
  const int s0 = (dir == 1 ? ec : ec+1);   // slot of the moved edge
  const int s1 = (dir == 1 ? ec+1 : ec);   // slot of the edge folded onto
  const edge* E0 = mmc->Edge(pc,s0).get();
  const edge* E1 = mmc->Edge(pc,s1).get();

  int t_pr = -1, t_pre = -1;
  multigon* t_mm = E1->target_multigon(mmc,t_pr,t_pre);
  const int k = t_mm->prongs();
  const int t2_pr = traintracks::mod(t_pr+dir,k);
  const int mc = multigon_index(mmc);
  const int mt = multigon_index(t_mm);

  // Old numbers of the two edges and three prongs involved.
  std::map<const edge*,int> old_number;
  for (int e = 0; e < B.nedges(); ++e) old_number[B.edge_ptr[e]] = e;
  const int e0 = old_number.at(E0);
  const int e1 = old_number.at(E1);
  const int q_cusp = B.prong_number[mc][pc];
  const int q_t = B.prong_number[mt][t_pr];
  const int q_t2 = B.prong_number[mt][t2_pr];

  // Which end of each edge sits at the cusp, in the old orientation.
  const bool e0_head_at_cusp = (B.edge_head[e0] == q_cusp);
  const bool e1_head_at_cusp = (B.edge_head[e1] == q_cusp);
  if (!e0_head_at_cusp && B.edge_tail[e0] != q_cusp) fail("Moved edge not at cusp");
  if (!e1_head_at_cusp && B.edge_tail[e1] != q_cusp) fail("Onto edge not at cusp");

  // Predicted post-fold signature of every old prong: E0 leaves the cusp
  // prong and joins the target prong t2.
  std::vector<prong_signature> sig_old(B.nprongs());
  {
    std::vector<std::vector<const edge*> > mult_edges(multigons());
    for (int m = 0; m < multigons(); ++m)
      for (int p = 0; p < Multigon(m).prongs(); ++p)
        for (int e = 0; e < Multigon(m).edges(p); ++e)
          mult_edges[m].push_back(Multigon(m).Edge(p,e).get());
    // Move E0 from mc to mt at the multigon level.
    mult_edges[mc].erase(std::find(mult_edges[mc].begin(),mult_edges[mc].end(),E0));
    mult_edges[mt].push_back(E0);
    for (int m = 0; m < multigons(); ++m)
      std::sort(mult_edges[m].begin(),mult_edges[m].end());

    for (int q = 0; q < B.nprongs(); ++q)
      {
        const int m = B.prong[q].multigon, p = B.prong[q].prong;
        prong_signature& sg = sig_old[q];
        sg.nprongs = B.prong[q].nprongs;
        sg.punctured = B.prong[q].punctured;
        for (int e = 0; e < Multigon(m).edges(p); ++e)
          {
            const edge* E = Multigon(m).Edge(p,e).get();
            if (E == E0 && q == q_cusp) continue;
            sg.at_prong.push_back(E);
          }
        if (q == q_t2) sg.at_prong.push_back(E0);
        std::sort(sg.at_prong.begin(),sg.at_prong.end());
        sg.at_multigon = mult_edges[m];
      }
  }

  // The fold itself (normalises the track).
  if (!fold(*mmc,pc,ec,dir)) return false;

  fm.after = numbering();
  const ttnumbering& A = fm.after;
  if (A.nedges() != fm.nmain || A.nprongs() != fm.nsides)
    fail("Edge or prong count changed");

  // Match old prongs to new prongs by signature.
  std::map<prong_signature,int> new_by_sig;
  {
    std::vector<std::vector<const edge*> > mult_edges(multigons());
    for (int m = 0; m < multigons(); ++m)
      {
        for (int p = 0; p < Multigon(m).prongs(); ++p)
          for (int e = 0; e < Multigon(m).edges(p); ++e)
            mult_edges[m].push_back(Multigon(m).Edge(p,e).get());
        std::sort(mult_edges[m].begin(),mult_edges[m].end());
      }
    for (int q = 0; q < A.nprongs(); ++q)
      {
        const int m = A.prong[q].multigon, p = A.prong[q].prong;
        prong_signature sg;
        sg.nprongs = A.prong[q].nprongs;
        sg.punctured = A.prong[q].punctured;
        for (int e = 0; e < Multigon(m).edges(p); ++e)
          sg.at_prong.push_back(Multigon(m).Edge(p,e).get());
        std::sort(sg.at_prong.begin(),sg.at_prong.end());
        sg.at_multigon = mult_edges[m];
        if (!new_by_sig.insert(std::make_pair(sg,q)).second)
          fail("Two prongs with the same signature after fold");
      }
  }
  fm.prong_image.assign(B.nprongs(),-1);
  std::set<int> used;
  for (int q = 0; q < B.nprongs(); ++q)
    {
      auto it = new_by_sig.find(sig_old[q]);
      if (it == new_by_sig.end()) fail("Cannot match prong after fold");
      fm.prong_image[q] = it->second;
      if (!used.insert(it->second).second) fail("Prong matched twice after fold");
    }

  // Match old edges to new edges by object identity, with orientation.
  std::map<const edge*,int> new_number;
  for (int e = 0; e < A.nedges(); ++e) new_number[A.edge_ptr[e]] = e;
  fm.edge_image.assign(B.nedges(),0);
  for (int e = 0; e < B.nedges(); ++e)
    {
      auto it = new_number.find(B.edge_ptr[e]);
      if (it == new_number.end()) fail("Edge object lost in fold");
      const int en = it->second;
      // Geometric tail of the old edge after the fold: the moved end of
      // E0 now sits at the target prong t2.
      int q_tail = B.edge_tail[e];
      int q_head = B.edge_head[e];
      if (e == e0) { if (e0_head_at_cusp) q_head = q_t2; else q_tail = q_t2; }
      const int nt = fm.prong_image[q_tail], nh = fm.prong_image[q_head];
      if (A.edge_tail[en] == nt && A.edge_head[en] == nh)
        fm.edge_image[e] = A.main_letter(en);
      else if (A.edge_tail[en] == nh && A.edge_head[en] == nt)
        fm.edge_image[e] = -A.main_letter(en);
      else
        fail("Edge endpoints do not match after fold");
    }

  // The moved edge's image path.
  fm.moved = e0;
  fm.onto = e1;
  fm.dir = dir;
  fm.cusp_prong = q_cusp;
  fm.target_from = fm.prong_image[q_t];
  fm.target_to = fm.prong_image[q_t2];

  // Side letters for the two traversal directions between the target
  // prongs t and t2 = t + dir (mod k).  Side q runs from prong q to the
  // next prong in cycle_edges order.
  const int side_t_to_t2 = (dir == 1 ? A.side_letter(fm.target_from)
                                     : -A.side_letter(fm.target_to));
  const int side_t2_to_t = -side_t_to_t2;

  auto relabel = [&](const int old_letter) {
    const int e = std::abs(old_letter) - 1;
    return (old_letter > 0 ? fm.edge_image[e] : -fm.edge_image[e]);
  };
  // `onto` traversed from the target prong t to the cusp, in old letters.
  const int e1_t_to_cusp = (e1_head_at_cusp ? B.main_letter(e1) : -B.main_letter(e1));

  fm.moved_word.clear();
  if (e0_head_at_cusp)
    {
      // Old `moved` runs from its far end to the cusp: new copy, then the
      // side from t2 to t, then `onto` from t to the cusp.
      fm.moved_first = true;
      fm.side = side_t2_to_t;
      fm.moved_word.push_back(relabel(B.main_letter(e0)));
      fm.moved_word.push_back(fm.side);
      fm.moved_word.push_back(relabel(e1_t_to_cusp));
    }
  else
    {
      // Old `moved` runs from the cusp to its far end: `onto` from the
      // cusp to t, the side from t to t2, then the new copy.
      fm.moved_first = false;
      fm.side = side_t_to_t2;
      fm.moved_word.push_back(relabel(-e1_t_to_cusp));
      fm.moved_word.push_back(fm.side);
      fm.moved_word.push_back(relabel(B.main_letter(e0)));
    }

  // Continuity of the three-letter path in the new numbering.
  for (int i = 0; i + 1 < 3; ++i)
    {
      if (A.head_of(fm.moved_word[i]) != A.tail_of(fm.moved_word[i+1]))
        fail("Image of the moved edge is not a continuous path");
    }

  return true;
}


jlt::freeauto<int> identity_traintrack_map(const int nmain, const int nsides)
{
  return jlt::freeauto<int>(nmain + nsides);
}


jlt::freeauto<int> fold_map_data::to_freeauto() const
{
  jlt::freeauto<int> AM(nmain + nsides);

  for (int e = 0; e < nmain; ++e)
    {
      if (e == moved)
        {
          jlt::freeword<int> w;
          for (auto x : moved_word) w.push_back(x);
          AM[e+1] = w;
        }
      else
        {
          AM[e+1] = jlt::freeword<int>({edge_image[e]});
        }
    }
  for (int q = 0; q < nsides; ++q)
    {
      AM[nmain+1+q] = jlt::freeword<int>({nmain+1+prong_image[q]});
    }

  return AM;
}


std::ostream& fold_map_data::print(std::ostream& strm) const
{
  strm << "fold: edge " << moved+1 << " onto edge " << onto+1
       << " (dir " << dir << ") at cusp prong " << cusp_prong
       << "; target prongs " << target_from << " -> " << target_to
       << "; side letter " << side << "\n";
  strm << "  image of edge " << moved+1 << ":";
  for (auto x : moved_word) strm << " " << x;
  strm << "\n  edge images:";
  for (int e = 0; e < nmain; ++e)
    if (e != moved) strm << " " << e+1 << "->" << edge_image[e];
  strm << "\n  prong images:";
  for (int q = 0; q < nsides; ++q) strm << " " << q << "->" << prong_image[q];
  return strm << "\n";
}

} // namespace traintracks
