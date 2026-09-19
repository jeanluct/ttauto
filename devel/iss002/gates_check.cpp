// gates_check.cpp -- scratch diagnostic for issue #2 (not part of the build).
//
// For a closed path in the n=6, stratum-5, 90-vertex automaton, this
// program rebuilds the word-level train-track map with geometrically
// consistent conventions (edges oriented, the polygon side actually
// traversed by each fold, prongs tracked through normalisation and
// identified with the graph's stored vertex the way the automaton does),
// then runs the Bestvina-Handel gate connectivity test (BH Prop. 3.3.2):
// vertices = unpunctured multigons and prongs of punctured multigons,
// directions = main-edge ends plus the two peripheral-loop directions,
// gates = classes under D^k, joins = turns actually taken by iterated
// images.  A disconnected gate graph means the matrix criterion has
// accepted a reducible map.
//
// Build (after the normal CMake build, from the repo root):
//   g++ -std=c++17 -O1 -Iinclude -Iextern/jlt -Iextern/jlt/extern/CSparse/Include \
//       devel/iss002/gates_check.cpp lib/libttauto.a build/libcsparse.a -o /tmp/gates_check
// Run: with no arguments it analyses the bad cycle {28,45,42,70,87,84,28}.
// Otherwise pass the 0-based vertex cycle, then "M", then the n*n transition
// matrix entries row by row (as in pA_n=6_5_1_inv.m) to select the branch
// sequence.  Written by Claude with J-L Thiffeault, 2026-09-19.
#include <algorithm>
#include <cassert>
#include <functional>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

#include <jlt/freeauto.hpp>
#include <jlt/freeword.hpp>
#include <jlt/mathmatrix.hpp>

#include "traintracks/build.hpp"
#include "traintracks/map.hpp"
#include "traintracks/traintrack.hpp"
#include "traintracks/util.hpp"
#include "ttauto/folding_path.hpp"
#include "ttauto/ttfoldgraph.hpp"

using namespace traintracks;
using ttauto::ttfoldgraph;
typedef std::pair<int,int> MP;   // (multigon, prong)

static const multigon& M(const traintrack& t, int m) { return t.Multigon(m); }
static int midx(const traintrack& tt, const multigon* mm)
{
  for (int m = 0; m < tt.multigons(); ++m) if (&M(tt,m) == mm) return m;
  std::cerr << "multigon not found\n"; std::exit(1);
}
static int pidx(const traintrack& tt, int m, int p)
{
  int ix = 0; for (int j = 0; j < m; ++j) ix += M(tt,j).prongs();
  return ix + p;
}
static MP from_pidx(const traintrack& tt, int ix)
{
  for (int m = 0; m < tt.multigons(); ++m)
    { int k = M(tt,m).prongs(); if (ix < k) return MP(m,ix); ix -= k; }
  std::cerr << "prong index out of range\n"; std::exit(1);
}
static std::string mpstr(MP v) { return "(" + std::to_string(v.first) + "," + std::to_string(v.second) + ")"; }
static MP ending(const traintrack& tt, const edge* E, int en)
{ return MP(midx(tt,E->ending_multigon(en)),E->ending_prong(en)); }

// Map from edge pointer to 0-based label (coding/weight order) and back.
static void edge_labels(traintrack& tt, std::map<const edge*,int>& lab, std::vector<const edge*>& ptr)
{
  const int n = tt.edges();
  traintrack::dblVec wv(n); for (int e = 0; e < n; ++e) wv[e] = e+1;
  tt.weights(wv.begin());
  lab.clear(); ptr.assign(n,0);
  for (int m = 0; m < tt.multigons(); ++m)
    for (int p = 0; p < M(tt,m).prongs(); ++p)
      for (int e = 0; e < M(tt,m).edges(p); ++e)
        {
          const edge* E = M(tt,m).Edge(p,e).get();
          int l = (int)(E->weight()+0.5) - 1;
          assert(l >= 0 && l < n);
          lab[E] = l; ptr[l] = E;
        }
  assert((int)lab.size() == n);
}

int main(int argc, char** argv)
{
  std::cout << std::unitbuf;
  const int n_punct = 6, trk = 4, sgidx = 0;
  jlt::vector<traintrack> ttv = build_traintrack_list(n_punct);
  ttfoldgraph<traintrack> full(ttv[trk]);
  std::list<ttfoldgraph<traintrack> > sgs = ttauto::subgraphs(full);
  ttauto::prune_multihumps(sgs);
  auto it = sgs.begin(); std::advance(it,sgidx);
  const ttfoldgraph<traintrack>& ttg = *it;
  assert(ttg.vertices() == 90);

  // Cycle (0-based vertices) from argv, else the bad cycle.  Branches are found
  // by searching all branch sequences realising the vertex sequence; if a
  // matrix is given (row-major, n*n entries after the cycle, separated by 'M'),
  // pick the sequence whose composite matrix matches.
  std::vector<int> cyc; std::vector<int> want;
  { bool inM = false; for (int a = 1; a < argc; ++a) { std::string t(argv[a]); if (t == "M") { inM = true; continue; } (inM ? want : cyc).push_back(std::stoi(t)); } }
  if (cyc.empty()) { cyc = {28,45,42,70,87,84,28}; }
  const int L = cyc.size()-1;
  std::vector<int> br(L,-1);
  {
    const int nn = ttg.edges();
    std::vector<std::vector<int> > cands;
    std::vector<int> cur(L);
    std::function<void(int)> rec = [&](int i) {
      if (i == L) { cands.push_back(cur); return; }
      for (int b = 0; b < ttg.foldings(cyc[i]); ++b) if (ttg.target_vertex(cyc[i],b) == cyc[i+1]) { cur[i] = b; rec(i+1); }
    };
    rec(0);
    std::cout << cands.size() << " branch sequence(s) realise the vertex cycle";
    std::vector<std::vector<int> > ok;
    for (auto& c : cands)
      {
        ttauto::folding_path<traintrack> p(ttg,cyc[0]); for (int i = 0; i < L; ++i) p.push_back(c[i]);
        jlt::mathmatrix<int> A = p.transition_matrix();
        bool match = want.empty();
        if (!want.empty()) { match = true; for (int r = 0; r < nn; ++r) for (int c2 = 0; c2 < nn; ++c2) if (A(r,c2) != want[r*nn+c2]) match = false; }
        if (match) ok.push_back(c);
      }
    std::cout << ", " << ok.size() << " match the requested matrix\n";
    assert(!ok.empty());
    br = ok[0];
    std::cout << "branches:"; for (int b : br) std::cout << " " << b; std::cout << "\n";
  }

  const int n = ttg.edges();
  const int nsides = ttg.traintrack(cyc[0]).total_prongs();
  const int ngen = n + nsides;
  const jlt::mathmatrix<int> id = jlt::identity_matrix<int>(n);


  const traintrack tt0(ttg.traintrack(cyc[0]));
  std::vector<std::pair<MP,MP> > ends0(n);
  {
    traintrack T(tt0); std::map<const edge*,int> l; std::vector<const edge*> p; edge_labels(T,l,p);
    for (int e = 0; e < n; ++e) ends0[e] = std::make_pair(ending(T,p[e],0),ending(T,p[e],1));
  }

  // Fingerprint of a (multigon,prong): type + labels attached in slot order.
  auto fingerprint = [&](const traintrack& T, const std::map<const edge*,int>& lab, int m, int p) {
    std::vector<int> fp; fp.push_back(M(T,m).punctured()); const int km = M(T,m).prongs(); fp.push_back(km);
    for (int q = 0; q < km; ++q)
      {
        const int pq = (p+q) % km; fp.push_back(-1);
        for (int e = 0; e < M(T,m).edges(pq); ++e) fp.push_back(lab.at(M(T,m).Edge(pq,e).get()));
      }
    return fp;
  };

  jlt::freeauto<int> AMtot(ngen);

  for (int i = 0; i < L; ++i)
    {
      traintrack tt(ttg.traintrack(cyc[i]));   // fresh copy of the graph's stored object
      int f = -1, b = 0;
      for (int ff = 0; ff < tt.foldings(); ++ff)
        {
          traintrack t2(tt);
          jlt::freeauto<int> AM = t2.fold_traintrack_map(ff);
          if (transition_matrix_from_map(tt,AM) == id) continue;
          if (b == br[i]) { f = ff; break; }
          ++b;
        }
      assert(f >= 0);

      std::map<const edge*,int> labA; std::vector<const edge*> ptrA;
      edge_labels(tt,labA,ptrA);

      multigon* mmc = 0; int pc = -1, ec = -1;
      tt.fold_cusp_location(f,mmc,pc,ec);
      const int dir = 1 - 2*(f % 2);
      const int s0 = (dir == 1 ? ec : ec+1), s1 = (dir == 1 ? ec+1 : ec);
      const edge* E0 = mmc->Edge(pc,s0).get();
      const edge* E1 = mmc->Edge(pc,s1).get();
      int t_pr = -1, t_pre = -1;
      multigon* t_mm = E1->target_multigon(mmc,t_pr,t_pre);
      const int k = t_mm->prongs();
      const int t2_pr = traintracks::mod(t_pr+dir,k);
      const int mc = midx(tt,mmc), mt = midx(tt,t_mm);
      int en0 = -1, en1 = -1;
      for (int en = 0; en < 2; ++en)
        {
          if (E0->ending_multigon(en) == mmc && E0->ending_prong(en) == pc) en0 = en;
          if (E1->ending_multigon(en) == mmc && E1->ending_prong(en) == pc) en1 = en;
        }
      assert(en0 >= 0 && en1 >= 0);
      assert(ending(tt,E1,1-en1) == MP(mt,t_pr));

      std::cout << "step " << i << ": vertex " << cyc[i]+1 << " -> " << cyc[i+1]+1
                << ", fold f=" << f << " dir=" << dir
                << ": edge " << labA[E0]+1 << " folded onto edge " << labA[E1]+1
                << " at " << (M(tt,mc).punctured() ? "punctured " : "unpunctured ") << M(tt,mc).prongs() << "-gon " << mc
                << " (" << M(tt,mc).edges(pc) << " edges at prong); target " << (M(tt,mt).punctured() ? "punctured " : "unpunctured ") << k << "-gon " << mt
                << ", side prong " << t_pr << " -> " << t2_pr
                << "; code's inf. gen. = multigon " << from_pidx(tt,tt.fold_infinitesimal_generator(f,n)-n-1).first << "\n";

      // Predicted post-fold slot contents (edge objects) for every pre-fold prong.
      // Pre-fold endings of every edge as (m,p), for orientation bookkeeping.
      std::vector<std::vector<const edge*> > pred(nsides);
      std::vector<std::pair<MP,MP> > preEnds(n);
      for (int e = 0; e < n; ++e) preEnds[e] = std::make_pair(ending(tt,ptrA[e],0),ending(tt,ptrA[e],1));
      for (int m = 0; m < tt.multigons(); ++m)
        for (int p = 0; p < M(tt,m).prongs(); ++p)
          {
            std::vector<const edge*> c;
            for (int e = 0; e < M(tt,m).edges(p); ++e)
              {
                const edge* E = M(tt,m).Edge(p,e).get();
                if (E == E0 && m == mc && p == pc) continue;
                c.push_back(E);
              }
            if (m == mt && p == t2_pr) { if (dir == 1) c.insert(c.begin(),E0); else c.push_back(E0); }
            pred[pidx(tt,m,p)] = c;
          }
      const int nmult = tt.multigons();
      std::vector<int> preprongs(nmult); for (int m = 0; m < nmult; ++m) preprongs[m] = M(tt,m).prongs();
      auto side_step = [&](int from, int step) {
        return step == 1 ? std::make_pair(pidx(tt,mt,from),1) : std::make_pair(pidx(tt,mt,traintracks::mod(from-1,k)),-1);
      };
      const std::pair<int,int> S_t_to_t2 = side_step(t_pr,dir);
      const std::pair<int,int> S_t2_to_t = side_step(t2_pr,-dir);

      bool ok = tt.fold(f); assert(ok);
      assert(tt == ttg.traintrack(cyc[i+1]));

      // Match pre-fold prongs to post-fold prongs by slot contents (edge objects).
      std::vector<int> pi(nsides,-1);
      {
        std::vector<int> mmatch(nmult,-1);
        for (int m = 0; m < nmult; ++m)
          {
            std::set<const edge*> want;
            for (int p = 0; p < preprongs[m]; ++p) for (auto E : pred[pidx(ttg.traintrack(cyc[i]),m,p)]) want.insert(E);
            for (int m2 = 0; m2 < tt.multigons() && mmatch[m] < 0; ++m2)
              {
                if (M(tt,m2).prongs() != preprongs[m]) continue;
                std::set<const edge*> have;
                for (int p = 0; p < M(tt,m2).prongs(); ++p) for (int e = 0; e < M(tt,m2).edges(p); ++e) have.insert(M(tt,m2).Edge(p,e).get());
                if (have == want) mmatch[m] = m2;
              }
            if (mmatch[m] < 0) { std::cout << "   ERROR: cannot match multigon " << m << " after fold\n"; std::exit(1); }
            // rotation
            const int km = preprongs[m]; int rot = -1;
            for (int r = 0; r < km && rot < 0; ++r)
              {
                bool okr = true;
                for (int p = 0; p < km && okr; ++p)
                  {
                    const std::vector<const edge*>& c = pred[pidx(ttg.traintrack(cyc[i]),m,p)];
                    const int p2 = (p + r) % km;
                    if (M(tt,mmatch[m]).edges(p2) != (int)c.size()) { okr = false; break; }
                    for (int e = 0; e < (int)c.size(); ++e) if (M(tt,mmatch[m]).Edge(p2,e).get() != c[e]) okr = false;
                  }
                if (okr) rot = r;
              }
            if (rot < 0) { std::cout << "   ERROR: cannot match prongs of multigon " << m << " after fold\n"; std::exit(1); }
            for (int p = 0; p < km; ++p) pi[pidx(ttg.traintrack(cyc[i]),m,p)] = pidx(tt,mmatch[m],(p+rot)%km);
          }
        std::set<int> img(pi.begin(),pi.end()); assert((int)img.size() == nsides);
      }

      std::map<const edge*,int> labB; std::vector<const edge*> ptrB;
      edge_labels(tt,labB,ptrB);

      // Match the folded copy tt with the graph's stored object Gn at the next vertex,
      // using edge labels as the identification (this is what the automaton does).
      traintrack Gn(ttg.traintrack(cyc[i+1]));
      std::map<const edge*,int> labG; std::vector<const edge*> ptrG; edge_labels(Gn,labG,ptrG);
      std::map<std::vector<int>,int> fpG;
      for (int m = 0; m < Gn.multigons(); ++m)
        for (int p = 0; p < M(Gn,m).prongs(); ++p)
          { auto fp = fingerprint(Gn,labG,m,p); assert(!fpG.count(fp)); fpG[fp] = pidx(Gn,m,p); }
      std::vector<int> F2G(nsides,-1);
      for (int m = 0; m < tt.multigons(); ++m)
        for (int p = 0; p < M(tt,m).prongs(); ++p)
          {
            auto fp = fingerprint(tt,labB,m,p);
            if (!fpG.count(fp)) { std::cout << "   ERROR: no matching prong in graph object for fingerprint of (" << m << "," << p << ")\n"; std::exit(1); }
            F2G[pidx(tt,m,p)] = fpG[fp];
          }
      { std::set<int> img(F2G.begin(),F2G.end()); assert((int)img.size() == nsides); }
      // rotation compatibility
      for (int m = 0; m < tt.multigons(); ++m)
        {
          const int km = M(tt,m).prongs();
          for (int p = 0; p < km; ++p)
            {
              MP a = from_pidx(Gn,F2G[pidx(tt,m,p)]), c = from_pidx(Gn,F2G[pidx(tt,m,(p+1)%km)]);
              if (!(a.first == c.first && c.second == (a.second+1) % km)) std::cout << "   WARNING: prong matching not a rotation at multigon " << m << "\n";
            }
        }
      // Edge orientation: pre-fold ending0 -> ending1, transported through the fold, vs Gn's ending order.
      std::vector<int> sgn(n,0);
      bool relabelled = false;
      const traintrack& Gi = ttg.traintrack(cyc[i]);
      for (int e = 0; e < n; ++e)
        {
          const edge* E = ptrA[e];
          MP A = preEnds[e].first, B = preEnds[e].second;
          if (E == E0) { if (en0 == 0) A = MP(mt,t2_pr); else B = MP(mt,t2_pr); }
          MP a0 = from_pidx(Gn,F2G[pi[pidx(Gi,A.first,A.second)]]);
          MP a1 = from_pidx(Gn,F2G[pi[pidx(Gi,B.first,B.second)]]);
          const int eg = labB[E];   // post-fold label of this edge object
          MP g0 = ending(Gn,ptrG[eg],0), g1 = ending(Gn,ptrG[eg],1);
          if (a0 == g0 && a1 == g1) sgn[eg] = 1;
          else if (a0 == g1 && a1 == g0) { sgn[eg] = -1; relabelled = true; }
          else { std::cout << "   ERROR: edge " << eg+1 << " endings do not match graph object\n"; std::exit(1); }
        }
      if (relabelled) { std::cout << "   note: ending order differs from graph object for edges:"; for (int e = 0; e < n; ++e) if (sgn[e] < 0) std::cout << " " << e+1; std::cout << "\n"; }

      auto sideletter = [&](std::pair<int,int> s) { return s.second * (n + 1 + F2G[pi[s.first]]); };
      auto mainletter = [&](int x) { return sgn[std::abs(x)-1] * x; };

      jlt::freeauto<int> AMc(ngen);
      for (int e = 0; e < n; ++e)
        {
          const edge* E = ptrA[e];
          if (E != E0) { AMc[e+1] = jlt::freeword<int>({mainletter(labB[E]+1)}); continue; }
          const int e0n = labB[E0]+1, e1n = labB[E1]+1;
          jlt::freeword<int> w;
          if (en0 == 1)
            { w.push_back(mainletter(e0n)); w.push_back(sideletter(S_t2_to_t)); w.push_back(mainletter(en1 == 1 ? e1n : -e1n)); }
          else
            { w.push_back(mainletter(en1 == 0 ? e1n : -e1n)); w.push_back(sideletter(S_t_to_t2)); w.push_back(mainletter(e0n)); }
          AMc[e+1] = w;
        }
      for (int s = 0; s < nsides; ++s) AMc[n+1+s] = jlt::freeword<int>({n+1+F2G[pi[s]]});
      std::cout << "   one-step map:";
      for (int g = 1; g <= n; ++g) { std::cout << " " << g << "->("; for (auto x : AMc.get_action(g)) std::cout << x << " "; std::cout << ")"; }
      std::cout << "\n";
      {
        jlt::mathmatrix<int> A = ttg.transition_matrix(cyc[i],br[i]).full(), B = transition_matrix_from_map(ttg.traintrack(cyc[i]),AMc);
        if (A != B) { std::cout << "   ONE-STEP MATRIX MISMATCH vs graph\n"; std::exit(1); }
      }
      AMtot = AMtot * AMc;
    }

  // Abelianisation check against the path's transition matrix.
  {
    ttauto::folding_path<traintrack> p(ttg,cyc[0]);
    for (int i = 0; i < L; ++i) p.push_back(br[i]);
    jlt::mathmatrix<int> A = p.transition_matrix(), B = transition_matrix_from_map(tt0,AMtot);
    if (A != B) { std::cout << "MATRIX MISMATCH\n"; std::exit(1); } else std::cout << "\nComposite abelianises to the path's transition matrix: OK\n";
  }

  std::cout << "\n=== Train track at vertex " << cyc[0]+1 << "\n";
  tt0.print(std::cout);
  for (int e = 0; e < n; ++e)
    std::cout << "  main edge " << e+1 << ": " << mpstr(ends0[e].first) << " -> " << mpstr(ends0[e].second) << "\n";
  for (int s = 0; s < nsides; ++s)
    {
      MP a = from_pidx(tt0,s);
      std::cout << "  side " << n+1+s << ": multigon " << a.first << " prong " << a.second << " -> "
                << (a.second+1) % M(tt0,a.first).prongs()
                << (M(tt0,a.first).punctured() ? " (peripheral)" : " (infinitesimal)") << "\n";
    }
  std::cout << "\n=== Composed map along the cycle (geometric words):\n" << AMtot;

  // ---- Orientation/continuity check of words in AMtot^k. ----
  auto side_from_to = [&](int S, MP& from, MP& to) {
    MP a = from_pidx(tt0,std::abs(S)-n-1);
    MP b(a.first,(a.second+1) % M(tt0,a.first).prongs());
    if (S > 0) { from = a; to = b; } else { from = b; to = a; }
  };
  std::map<int,MP> head;   // head of signed main letter
  for (int e = 1; e <= n; ++e) { head[e] = ends0[e-1].second; head[-e] = ends0[e-1].first; }
  int bad = 0;
  {
    jlt::freeauto<int> P(ngen);
    for (int kk = 1; kk <= 3; ++kk)
      {
        P = P * AMtot;
        for (int g = 1; g <= n; ++g)
          {
            jlt::freeword<int> w = P.get_action(g);
            std::vector<int> v(w.begin(),w.end());
            for (size_t j = 0; j+1 < v.size(); ++j)
              {
                MP h, t;
                if (std::abs(v[j]) <= n) h = head[v[j]]; else { MP f_, t_; side_from_to(v[j],f_,t_); h = t_; }
                if (std::abs(v[j+1]) <= n) t = head[-v[j+1]]; else { MP f_, t_; side_from_to(v[j+1],f_,t_); t = f_; }
                if (h != t) { if (bad < 10) std::cout << "  DISCONTINUITY in image of " << g << " (iterate " << kk << "): " << v[j] << " ends at " << mpstr(h) << ", " << v[j+1] << " starts at " << mpstr(t) << "\n"; ++bad; }
                if (std::abs(v[j]) <= n && std::abs(v[j+1]) <= n) { if (bad < 10) std::cout << "  two consecutive main letters " << v[j] << " " << v[j+1] << " in image of " << g << "\n"; ++bad; }
              }
          }
      }
  }
  std::cout << "\nContinuity check of image paths (3 iterates): " << (bad ? std::to_string(bad) + " problems" : "all consistent") << "\n";

  // ---- Bestvina-Handel analysis. ----
  auto bh_vertex = [&](MP mp) { return M(tt0,mp.first).punctured() ? mp : MP(mp.first,-1); };
  std::map<MP,std::vector<int> > dirs;
  std::map<int,MP> start;
  for (int g = -n; g <= n; ++g) if (g) { start[g] = bh_vertex(head[-g]); dirs[start[g]].push_back(g); }
  for (int s = 0; s < nsides; ++s)
    {
      MP a = from_pidx(tt0,s);
      if (!M(tt0,a.first).punctured()) continue;
      MP b(a.first,(a.second+1) % M(tt0,a.first).prongs());
      start[n+1+s] = a; dirs[a].push_back(n+1+s);
      start[-(n+1+s)] = b; dirs[b].push_back(-(n+1+s));
    }
  std::map<int,int> D;
  for (auto& kv : start)
    {
      int g = kv.first;
      jlt::freeword<int> w = AMtot.get_action(g);
      if (std::abs(g) > n) { assert(w.size() == 1); D[g] = *w.begin(); continue; }
      // first REAL letter: main letter, or peripheral side
      int d = 0;
      for (auto x : w) { if (std::abs(x) <= n || M(tt0,from_pidx(tt0,std::abs(x)-n-1).first).punctured()) { d = x; break; } }
      assert(d != 0); D[g] = d;
    }
  std::map<MP,MP> vimage;
  for (auto& kv : dirs)
    {
      std::set<MP> tgt; for (int g : kv.second) tgt.insert(start[D[g]]);
      if (tgt.size() != 1) std::cout << "  WARNING: D at vertex " << mpstr(kv.first) << " lands on " << tgt.size() << " vertices\n";
      vimage[kv.first] = *tgt.begin();
    }
  const int N = start.size();
  auto meet = [&](int a, int b) { int x = a, y = b; for (int kk = 0; kk < N*N+1; ++kk) { x = D[x]; y = D[y]; if (x == y) return true; } return false; };

  std::set<std::pair<int,int> > turns;
  auto add_turn = [&](int a, int b) { if (a > b) std::swap(a,b); turns.insert(std::make_pair(a,b)); };
  for (int g = 1; g <= n; ++g)
    {
      jlt::freeword<int> w = AMtot.get_action(g);
      std::vector<int> real;
      for (auto x : w)
        {
          if (std::abs(x) <= n) { real.push_back(x); continue; }
          if (M(tt0,from_pidx(tt0,std::abs(x)-n-1).first).punctured()) real.push_back(x);
        }
      for (size_t j = 0; j+1 < real.size(); ++j) add_turn(-real[j],real[j+1]);
    }
  std::cout << "Turns taken by one-step images:";
  for (auto t : turns) std::cout << " (" << t.first << "," << t.second << ")";
  std::cout << "\n";
  for (bool grew = true; grew;)
    {
      grew = false;
      std::vector<std::pair<int,int> > cur(turns.begin(),turns.end());
      for (auto t : cur)
        {
          int a = D[t.first], b = D[t.second]; if (a > b) std::swap(a,b);
          if (a == b) { std::cout << "  WARNING: turn (" << t.first << "," << t.second << ") collapses under D (map not efficient)\n"; continue; }
          if (turns.insert(std::make_pair(a,b)).second) grew = true;
        }
    }
  for (auto t : turns) if (start[t.first] != start[t.second]) std::cout << "  WARNING: turn " << t.first << "," << t.second << " at different vertices\n";

  std::cout << "\n=== Bestvina-Handel gate analysis at vertex " << cyc[0]+1 << "\n";
  std::cout << "(vertex (m,-1) = unpunctured multigon m; (m,p) = prong p of punctured multigon m;\n"
            << " directions: signed main edges 1.." << n << ", signed peripheral loops " << n+1 << "..)\n";
  bool all_connected = true;
  for (auto& kv : dirs)
    {
      const MP v = kv.first; std::vector<int> ds = kv.second; std::sort(ds.begin(),ds.end());
      std::vector<int> par(ds.size()); for (size_t j = 0; j < ds.size(); ++j) par[j] = j;
      std::function<int(int)> find = [&](int x) { return par[x] == x ? x : par[x] = find(par[x]); };
      for (size_t a = 0; a < ds.size(); ++a) for (size_t b = a+1; b < ds.size(); ++b) if (meet(ds[a],ds[b])) par[find(b)] = find(a);
      std::map<int,std::vector<int> > gates; std::map<int,int> gate_of;
      for (size_t j = 0; j < ds.size(); ++j) gates[find(j)].push_back(ds[j]);
      std::vector<std::vector<int> > gl; for (auto& g : gates) { for (int d : g.second) gate_of[d] = gl.size(); gl.push_back(g.second); }
      std::vector<int> gpar(gl.size()); for (size_t j = 0; j < gl.size(); ++j) gpar[j] = j;
      std::function<int(int)> gfind = [&](int x) { return gpar[x] == x ? x : gpar[x] = gfind(gpar[x]); };
      std::vector<std::pair<int,int> > joins;
      for (auto t : turns) if (start[t.first] == v)
        {
          int ga = gate_of[t.first], gb = gate_of[t.second];
          if (ga == gb) { std::cout << "  WARNING: realized turn (" << t.first << "," << t.second << ") inside a gate at " << mpstr(v) << "\n"; continue; }
          joins.push_back(t); gpar[gfind(gb)] = gfind(ga);
        }
      std::set<int> comps; for (size_t j = 0; j < gl.size(); ++j) comps.insert(gfind(j));
      const bool connected = comps.size() == 1;
      if (!connected) all_connected = false;
      std::cout << "Vertex " << mpstr(v) << (v.second < 0 ? " [unpunctured " + std::to_string(M(tt0,v.first).prongs()) + "-gon]" : " [punctured]")
                << " -> image vertex " << mpstr(vimage[v]) << "\n  gates:";
      for (auto& g : gl) { std::cout << " {"; for (size_t j = 0; j < g.size(); ++j) std::cout << (j ? " " : "") << g[j]; std::cout << "}"; }
      std::cout << "\n  joins (realized turns between gates):";
      for (auto t : joins) std::cout << " " << t.first << "~" << t.second;
      std::cout << "\n  gate graph " << (connected ? "CONNECTED" : "NOT CONNECTED (" + std::to_string(comps.size()) + " components)") << "\n";
    }
  std::cout << "\nOverall: " << (all_connected ? "all vertices connected -> pA criterion satisfied" : "some vertex disconnected -> reducible by BH criterion") << "\n";
  return 0;
}
