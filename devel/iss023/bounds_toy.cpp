// Pruning bounds for searches over products of fold matrices (issue #23).
//
// The fast version of bounds_toy.py, and independent of ttauto: it uses
// nothing but the standard library.  See bounds_toy.py for the model and
// pruning_bounds.tex for the mathematics.
//
// A generator is F = P + e_ab.  Words are multiplied on the left, as in
// ttauto, and every generator set is a one-vertex automaton, so every word
// is a closed path.  A word is accepted when its matrix is primitive with
// 1 < rho <= Lambda.  Run with criterion HS (Ham-Song norm only) to get the
// reference list, which is complete; any other criterion must reproduce its
// hash, or it has missed words.
//
// Spectral radii are bracketed, never estimated: for B >= 0 and x > 0,
// rho(B) >= min_i (Bx)_i / x_i always (Collatz-Wielandt), and
// rho(B) <= max_i (Bx)_i / x_i when B is irreducible.  Every prune uses a
// lower bracket, so abandoning a prefix is safe whatever the iteration
// count (it never discards a prefix that could still succeed); an
// acceptance uses the upper bracket, on a primitive matrix.  Iterating
// x <- (B + I) x tightens both.
//
// Build:  g++ -std=c++17 -O2 -o bounds_toy bounds_toy.cpp
// Run:    ./bounds_toy SET LAMBDA CRITERION [BUDGET]
//         SET is RL, rauzy3, perm3 or sym3; CRITERION is HS, C0, C1, C2,
//         C3 (the C2 decision computed fast) or BAD.
// Prints: accepted words, cyclic classes, visited, wasted, longest prefix,
// longest accepted word, and a hash of the accepted set for comparison.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <set>
#include <string>
#include <vector>

namespace {

const int NMAX = 8;
typedef std::array<std::array<long long,NMAX>,NMAX> Mat;
typedef std::array<int,NMAX> Perm;

int n = 0;

Mat ident() { Mat M{}; for (int i = 0; i < n; ++i) M[i][i] = 1; return M; }

Mat mul(const Mat& A, const Mat& B)
{
  Mat C{};
  for (int i = 0; i < n; ++i)
    for (int k = 0; k < n; ++k)
      if (A[i][k])
        for (int j = 0; j < n; ++j) C[i][j] += A[i][k]*B[k][j];
  return C;
}

Mat pmat(const Perm& p)    // P e_c = e_p[c]
{
  Mat P{};
  for (int c = 0; c < n; ++c) P[p[c]][c] = 1;
  return P;
}

Perm pcompose(const Perm& p, const Perm& q)   // p o q
{
  Perm r{};
  for (int c = 0; c < n; ++c) r[c] = p[q[c]];
  return r;
}

Perm pinverse(const Perm& p)
{
  Perm r{};
  for (int c = 0; c < n; ++c) r[p[c]] = c;
  return r;
}

Perm pident() { Perm p{}; for (int c = 0; c < n; ++c) p[c] = c; return p; }

// Collatz-Wielandt brackets for rho(B), from iterating x <- (B+I)x.
void bracket(const Mat& B, double& lo, double& hi, const int iters = 60)
{
  std::array<double,NMAX> x;
  for (int i = 0; i < n; ++i) x[i] = 1.0;
  lo = 0; hi = INFINITY;
  for (int it = 0; it < iters; ++it)
    {
      std::array<double,NMAX> y{};
      for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) y[i] += B[i][j]*x[j];
      double l = INFINITY, h = 0;
      for (int i = 0; i < n; ++i)
        {
          const double r = y[i]/x[i];
          l = std::min(l,r); h = std::max(h,r);
        }
      lo = std::max(lo,l);
      hi = std::min(hi,h);
      if (hi - lo < 1e-12*hi) break;
      double s = 0;
      for (int i = 0; i < n; ++i) { x[i] = x[i] + y[i]; s += x[i]; }
      for (int i = 0; i < n; ++i) x[i] /= s;
      for (int i = 0; i < n; ++i) if (x[i] < 1e-300) x[i] = 1e-300;
    }
}

// Is rho(B) > L?  Iterates like bracket(), but stops as soon as the lower
// end passes L (a lower bound, so the answer "yes" is always right), or
// the bracket has converged below it.
bool rho_exceeds(const Mat& B, const double L, const int iters = 60)
{
  std::array<double,NMAX> x;
  for (int i = 0; i < n; ++i) x[i] = 1.0;
  double hi = INFINITY;
  for (int it = 0; it < iters; ++it)
    {
      std::array<double,NMAX> y{};
      for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) y[i] += B[i][j]*x[j];
      double l = INFINITY, h = 0;
      for (int i = 0; i < n; ++i)
        {
          const double r = y[i]/x[i];
          l = std::min(l,r); h = std::max(h,r);
        }
      if (l > L) return true;
      hi = std::min(hi,h);
      if (hi <= L) return false;   // converged below: rho <= hi <= L
      double s = 0;
      for (int i = 0; i < n; ++i) { x[i] = x[i] + y[i]; s += x[i]; }
      for (int i = 0; i < n; ++i) x[i] /= s;
      for (int i = 0; i < n; ++i) if (x[i] < 1e-300) x[i] = 1e-300;
    }
  return false;
}

bool primitive(const Mat& M)
{
  Mat B{};
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) B[i][j] = M[i][j] > 0;
  Mat X = B;
  for (int k = 0; k <= n*n - 2*n + 2; ++k)
    {
      bool pos = true;
      for (int i = 0; i < n && pos; ++i)
        for (int j = 0; j < n && pos; ++j) if (!X[i][j]) pos = false;
      if (pos) return true;
      X = mul(X,B);
      for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) X[i][j] = X[i][j] > 0;
    }
  return false;
}

// adj[i][j] true means an edge j -> i.
bool strongly_connected(const std::array<std::array<bool,NMAX>,NMAX>& adj)
{
  for (int dir = 0; dir < 2; ++dir)
    {
      std::vector<int> stack(1,0);
      std::array<bool,NMAX> seen{}; seen[0] = true; int count = 1;
      while (!stack.empty())
        {
          const int u = stack.back(); stack.pop_back();
          for (int v = 0; v < n; ++v)
            {
              const bool e = dir ? adj[u][v] : adj[v][u];
              if (e && !seen[v]) { seen[v] = true; ++count; stack.push_back(v); }
            }
        }
      if (count < n) return false;
    }
  return true;
}

struct Gen
{
  std::string name;
  Perm p;
  int a, b;           // F = P + e_ab
  int ea, eb;         // F = P (I + e_{ea eb})
  Mat F;
};

Gen make_gen(const std::string& name, const Perm& p, const int a, const int b)
{
  Gen g; g.name = name; g.p = p; g.a = a; g.b = b;
  const Perm pi = pinverse(p);
  g.ea = pi[a]; g.eb = b;
  g.F = pmat(p); g.F[a][b] += 1;
  if (g.F[a][b] != 1) { std::cerr << "extra 1 lands on the permutation\n"; std::exit(1); }
  return g;
}

std::vector<Gen> gens;
std::vector<Perm> H;                       // group generated by the p's
std::vector<std::pair<int,int> > positions;  // conjugates of the e'
double lam = 0;
std::string criterion;

void setup(const std::string& set)
{
  gens.clear();
  if (set == "RL")
    {
      n = 2;
      gens.push_back(make_gen("L",pident(),1,0));
      gens.push_back(make_gen("R",pident(),0,1));
    }
  else if (set == "rauzy3")
    {
      n = 3;
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
          if (a != b)
            gens.push_back(make_gen("e" + std::to_string(a) + std::to_string(b),
                                    pident(),a,b));
    }
  else if (set == "perm3")
    {
      n = 3;
      const Perm cyc = {1,2,0};
      gens.push_back(make_gen("a",pident(),0,1));
      gens.push_back(make_gen("b",pident(),1,2));
      gens.push_back(make_gen("c",cyc,1,2));
    }
  else if (set == "sym3")
    {
      // As perm3, plus a shear followed by a transposition, so that the
      // folds' permutations generate all of S_3, as in ttauto's automata.
      n = 3;
      const Perm cyc = {1,2,0}, swp = {1,0,2};
      gens.push_back(make_gen("a",pident(),0,1));
      gens.push_back(make_gen("b",pident(),1,2));
      gens.push_back(make_gen("c",cyc,1,2));
      gens.push_back(make_gen("d",swp,2,0));
    }
  else if (set.size() == 4 && set.compare(0,3,"sym") == 0 && set[3] >= '4'
           && set[3] <= '8')
    {
      // symD, D = 4..8: shears e_{i,i+1} along a chain, plus a shear
      // followed by the D-cycle and one followed by a transposition, so
      // that the permutations generate all of S_D.
      n = set[3] - '0';
      for (int i = 0; i+1 < n; ++i)
        gens.push_back(make_gen("s" + std::to_string(i),pident(),i,i+1));
      Perm cyc{}, swp = pident();
      for (int c = 0; c < n; ++c) cyc[c] = (c+1) % n;
      std::swap(swp[0],swp[1]);
      gens.push_back(make_gen("c",cyc,1,n-1));
      gens.push_back(make_gen("t",swp,n-1,0));
    }
  else { std::cerr << "unknown set " << set << "\n"; std::exit(1); }

  // Group closure.
  H.assign(1,pident());
  for (std::size_t k = 0; k < H.size(); ++k)
    for (const Gen& g : gens)
      {
        const Perm r = pcompose(g.p,H[k]);
        if (std::find(H.begin(),H.end(),r) == H.end()) H.push_back(r);
      }
  std::set<std::pair<int,int> > pos;
  for (const Perm& q : H)
    {
      const Perm qi = pinverse(q);
      for (const Gen& g : gens) pos.insert(std::make_pair(qi[g.ea],qi[g.eb]));
    }
  positions.assign(pos.begin(),pos.end());
}

long long hs_bound() { return (long long)std::floor(std::pow(lam,(double)n)) + n - 1; }

bool c0_prune(const Mat& A)
{
  long long norm = 0, colmin = -1, rowmin = -1;
  for (int j = 0; j < n; ++j)
    {
      long long c = 0, r = 0;
      for (int i = 0; i < n; ++i) { c += A[i][j]; r += A[j][i]; }
      norm += c;
      if (colmin < 0 || c < colmin) colmin = c;
      if (rowmin < 0 || r < rowmin) rowmin = r;
    }
  return norm > hs_bound() || colmin > lam || rowmin > lam;
}

// The completion bound: prune when every (Q, S) that allows an irreducible
// product gives rho(Q (I + E_S) N) > Lambda.  Uses lower brackets only.
bool beta_exceeds(const Mat& N)
{
  const int m = (int)positions.size();
  for (int mask = 0; mask < (1 << m); ++mask)
    {
      // Minimality is not required for safety, only for speed: a
      // superset of an admissible S gives a larger rho.
      Mat X = ident();
      std::array<std::array<bool,NMAX>,NMAX> R{};
      for (int i = 0; i < n; ++i) R[i][i] = true;
      for (int t = 0; t < m; ++t)
        if (mask & (1 << t))
          { X[positions[t].first][positions[t].second] += 1;
            R[positions[t].first][positions[t].second] = true; }
      for (int k = 0; k < n; ++k)
        for (int i = 0; i < n; ++i)
          if (R[i][k]) for (int j = 0; j < n; ++j) if (R[k][j]) R[i][j] = true;
      X = mul(X,N);
      for (const Perm& q : H)
        {
          std::array<std::array<bool,NMAX>,NMAX> adj{};
          for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
              if (N[k][j] > 0)
                for (int i = 0; i < n; ++i) if (R[i][k]) adj[q[i]][j] = true;
          if (!strongly_connected(adj)) continue;
          double lo, hi;
          bracket(mul(pmat(q),X),lo,hi);
          if (lo <= lam*(1+1e-12)) return false;   // this (Q,S) is still possible
        }
    }
  return true;
}

// The same decision as beta_exceeds, computed without trying every pair.
//
// A nonnegative matrix has rho at least its largest diagonal entry.  The
// diagonal of Q'X is X_{sigma(k),k}, k = 0..n-1, with sigma = Q'^{-1}, and
// X = (I + E_S) N >= N.  So a pair whose diagonal has an entry above
// Lambda has rho above Lambda and can be skipped: only permutations sigma
// with every N_{sigma(k),k} <= Lambda need be tried (found by
// backtracking, column by column), and for each only the minimal
// admissible S whose diagonal stays <= Lambda (found by growing S one
// position at a time, and stopping as soon as it is admissible).  Each
// candidate left is then bracketed as before.  Everything skipped has rho
// above Lambda, so the decision is the same as beta_exceeds'.
long long fast_candidates = 0;   // spectral brackets computed, for the cost

// Effort budget per evaluation, in calls to grow() (0: unlimited).  When it
// runs out the evaluation gives up and keeps the prefix, which is always
// safe: the search just does more work than it had to.
long long budget = 0, spent = 0, gaveup = 0;
struct out_of_budget {};

// A cheap lower bound on rho(Q'X), from quantities that only grow as
// entries of X grow: the largest diagonal entry, the 2-cycles
// sqrt(Y_ij Y_ji), and the smallest row and column sums.  Since X only
// grows as S grows, once this exceeds Lambda no larger S can fit either.
bool cheap_exceeds(const Mat& X, const Perm& q)
{
  const double L = lam*(1+1e-12);
  Mat Y{};
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) Y[q[i]][j] = X[i][j];
  long long rmin = -1, cmin = -1;
  for (int i = 0; i < n; ++i)
    {
      if (Y[i][i] > L) return true;
      long long r = 0, c = 0;
      for (int j = 0; j < n; ++j)
        {
          r += Y[i][j]; c += Y[j][i];
          if (j > i && (double)Y[i][j]*Y[j][i] > L*L) return true;
        }
      if (rmin < 0 || r < rmin) rmin = r;
      if (cmin < 0 || c < cmin) cmin = c;
    }
  return rmin > L || cmin > L;
}

// Reachability in D(S) as bitmasks: bit k of reach[i] is set when k
// reaches i (every vertex reaching itself).  Adding the edge b -> a of
// e_ab: everything that reaches b now reaches everything a reaches.
typedef std::array<unsigned,NMAX> Reach;

Reach reach_identity()
{
  Reach r{};
  for (int i = 0; i < n; ++i) r[i] = 1u << i;
  return r;
}

void reach_add(Reach& r, const int a, const int b)
{
  const unsigned from = r[b];
  for (int i = 0; i < n; ++i) if (r[i] & (1u << a)) r[i] |= from;
}

// Is G_{Q',S} strongly connected?  It has an edge j -> q(i) whenever
// N_kj > 0 and k reaches i; ncol[j] has bit k set when N_kj > 0.
bool admissible(const std::array<unsigned,NMAX>& ncol, const Perm& q,
                const Reach& r)
{
  std::array<unsigned,NMAX> out{}, in{};
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      if (ncol[j] & r[i]) { out[j] |= 1u << q[i]; in[q[i]] |= 1u << j; }
  const unsigned all = (1u << n) - 1;
  for (int dir = 0; dir < 2; ++dir)
    {
      const std::array<unsigned,NMAX>& e = dir ? in : out;
      unsigned seen = 1u, frontier = 1u;
      while (frontier)
        {
          unsigned next = 0;
          for (int u = 0; u < n; ++u) if (frontier & (1u << u)) next |= e[u];
          frontier = next & ~seen;
          seen |= next;
        }
      if (seen != all) return false;
    }
  return true;
}

// Grow S from position index `from`, with X = (I + E_S) N already formed
// and r the reachability of D(S).  Returns true if some minimal
// admissible S found this way has rho <= Lambda.
bool grow(const Mat& N, const std::array<unsigned,NMAX>& ncol, const Perm& q,
          const Perm& sigma, Mat& X, const Reach& r,
          std::vector<std::pair<int,int> >& S, const std::size_t from)
{
  if (budget > 0 && ++spent > budget) throw out_of_budget();
  if (cheap_exceeds(X,q)) return false;
  if (admissible(ncol,q,r))
    {
      // Minimal only: if some element can be dropped and the rest is still
      // admissible, that smaller set is tried too and gives a smaller rho.
      for (std::size_t drop = 0; drop < S.size(); ++drop)
        {
          Reach t = reach_identity();
          for (std::size_t u = 0; u < S.size(); ++u)
            if (u != drop) reach_add(t,S[u].first,S[u].second);
          if (admissible(ncol,q,t)) return false;
        }
      ++fast_candidates;
      return !rho_exceeds(mul(pmat(q),X),lam*(1+1e-12));
    }
  for (std::size_t t = from; t < positions.size(); ++t)
    {
      const int a = positions[t].first, b = positions[t].second;
      // If b already reaches a, adding e_ab changes no reachability, now
      // or for any larger S: no minimal admissible set contains it.
      if (r[a] & (1u << b)) continue;
      // Adding e_ab adds row b of N to row a of X.  Row a meets the
      // diagonal of Q'X in column k with sigma(k) = a.
      int kk = -1;
      for (int k = 0; k < n; ++k) if (sigma[k] == a) kk = k;
      if (X[a][kk] + N[b][kk] > lam*(1+1e-12)) continue;
      for (int j = 0; j < n; ++j) X[a][j] += N[b][j];
      Reach r2 = r;
      reach_add(r2,a,b);
      S.push_back(positions[t]);
      const bool ok = grow(N,ncol,q,sigma,X,r2,S,t+1);
      S.pop_back();
      for (int j = 0; j < n; ++j) X[a][j] -= N[b][j];
      if (ok) return true;
    }
  return false;
}

// Backtrack over sigma, column by column, keeping N_{sigma(k),k} <= Lambda.
bool some_pair_fits(const Mat& N, Perm& sigma, std::array<bool,NMAX>& used,
                    const int k)
{
  if (k == n)
    {
      const Perm q = pinverse(sigma);
      if (std::find(H.begin(),H.end(),q) == H.end()) return false;
      Mat X = N;
      std::array<unsigned,NMAX> ncol{};
      for (int j = 0; j < n; ++j)
        for (int k = 0; k < n; ++k) if (N[k][j] > 0) ncol[j] |= 1u << k;
      std::vector<std::pair<int,int> > S;
      return grow(N,ncol,q,sigma,X,reach_identity(),S,0);
    }
  for (int row = 0; row < n; ++row)
    {
      if (used[row] || N[row][k] > lam*(1+1e-12)) continue;
      used[row] = true; sigma[k] = row;
      const bool ok = some_pair_fits(N,sigma,used,k+1);
      used[row] = false;
      if (ok) return true;
    }
  return false;
}

bool beta_exceeds_fast(const Mat& N)
{
  Perm sigma{};
  std::array<bool,NMAX> used{};
  spent = 0;
  try { return !some_pair_fits(N,sigma,used,0); }
  catch (const out_of_budget&) { ++gaveup; return false; }   // keep: safe
}

bool prune(const Mat& A, const Mat& N)
{
  if (criterion == "HS")
    {
      long long norm = 0;
      for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) norm += A[i][j];
      return norm > hs_bound();
    }
  if (c0_prune(A)) return true;
  if (criterion == "C1" || criterion == "BAD")
    {
      double lo, hi; bracket(A,lo,hi);
      if (lo > lam*(1+1e-12)) return true;
    }
  if (criterion == "C2" && beta_exceeds(N)) return true;
  if (criterion == "C3" && beta_exceeds_fast(N)) return true;
  return false;
}

long long visited = 0, wasted = 0;
std::size_t maxlen = 0, maxacc = 0;
std::set<std::vector<int> > accepted;

bool walk(std::vector<int>& word, const Mat& A, const Perm& Q, const Mat& N)
{
  ++visited;
  maxlen = std::max(maxlen,word.size());
  bool found = false;
  if (!word.empty() && primitive(A))
    {
      double lo, hi; bracket(A,lo,hi,400);
      if (hi <= lam*(1+1e-12) && lo > 1 + 1e-9)
        { accepted.insert(word); maxacc = std::max(maxacc,word.size()); found = true; }
      else if (lo <= lam*(1+1e-12) && hi > lam*(1+1e-12))
        { std::cerr << "ambiguous bracket at the window edge\n"; std::exit(2); }
    }
  const Perm Qi = pinverse(Q);
  for (int gi = 0; gi < (int)gens.size(); ++gi)
    {
      const Gen& g = gens[gi];
      const Mat A2 = mul(g.F,A);
      // F Q N = (P Q) (I + Q^{-1} e' Q) N.
      Mat E = ident(); E[Qi[g.ea]][Qi[g.eb]] += 1;
      const Mat N2 = mul(E,N);
      const Perm Q2 = pcompose(g.p,Q);
      if (prune(A2,N2)) continue;
      word.push_back(gi);
      if (walk(word,A2,Q2,N2)) found = true;
      word.pop_back();
    }
  if (!found) ++wasted;
  return found;
}

} // namespace

int main(int argc, char** argv)
{
  if (argc < 4) { std::cerr << "usage: bounds_toy SET LAMBDA CRITERION\n"; return 1; }
  setup(argv[1]);
  lam = std::atof(argv[2]);
  criterion = argv[3];
  if (argc > 4) budget = std::atoll(argv[4]);
  std::vector<int> word;
  walk(word,ident(),pident(),ident());

  std::set<std::vector<int> > classes;
  std::uint64_t h = 1469598103934665603ULL;
  for (const auto& w : accepted)
    {
      std::vector<int> best = w;
      for (std::size_t i = 1; i < w.size(); ++i)
        {
          std::vector<int> r(w.begin()+i,w.end()); r.insert(r.end(),w.begin(),w.begin()+i);
          best = std::min(best,r);
        }
      classes.insert(best);
      for (int x : w) { h ^= (std::uint64_t)(x+1); h *= 1099511628211ULL; }
      h ^= 0xff; h *= 1099511628211ULL;
    }
  std::cout << argv[1] << " Lambda=" << lam << " " << criterion
            << " accepted=" << accepted.size() << " classes=" << classes.size()
            << " visited=" << visited << " wasted=" << wasted
            << " longest_prefix=" << maxlen << " longest_accepted=" << maxacc
            << " hash=" << std::hex << h << std::dec
            << " brackets=" << fast_candidates << " gaveup=" << gaveup << "\n";
  return 0;
}
