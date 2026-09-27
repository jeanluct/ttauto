#!/usr/bin/env python3
"""Pruning bounds for searches over products of fold matrices (issue #23).

A matrix-only model of the ttauto search, independent of ttauto.  A
generator is a fold matrix F = P + e_ab: a permutation matrix P plus a
single 1.  Words are multiplied on the left, as in ttauto
(TMl.push(F * TMl.top())), so the word F_k ... F_1 has matrix
M = F_k ... F_1.  Each generator set is a one-vertex automaton, so every
word is a closed path.  A word is accepted when M is primitive and
1 < rho(M) <= Lambda.

The search is depth-first, and a criterion decides whether a prefix can
be abandoned.  Every criterion must be safe: it may only abandon a prefix
if no accepted word begins with it (the prefix itself included); an unsafe
one makes the search miss words.  The reference list is what the search
finds using the Ham-Song norm alone, which is guaranteed complete (Ham &
Song 2007, Lemma 3.1); every other criterion must find exactly that
list.

Criteria:
  C0  what ttauto does now: Ham-Song norm, minimum row sum, minimum
      column sum.
  C1  C0, plus rho(prefix) > Lambda.  Safe when every P is the identity,
      since then F >= I and rho can only grow.
  C2  C0, plus the completion bound beta(prefix) > Lambda (see
      pruning_bounds.tex).
  BAD C0, plus rho(prefix) > Lambda applied to a set with permutations,
      where it is not justified.  Faulty on purpose: it should miss words,
      which shows the comparison notices an unsafe criterion.

Usage: python3 bounds_toy.py            (all generator sets, a few Lambda)
"""

import itertools
import math
import sys
from collections import Counter

import numpy as np

sys.setrecursionlimit(100000)


# ---------------------------------------------------------------------------
# Matrices

def perm_matrix(p):
    """Permutation matrix with P e_c = e_p[c], i.e. P[p[c], c] = 1."""
    n = len(p)
    P = np.zeros((n, n), dtype=np.int64)
    for c in range(n):
        P[p[c], c] = 1
    return P


def unit(n, a, b):
    E = np.zeros((n, n), dtype=np.int64)
    E[a, b] = 1
    return E


def rho(M):
    return float(max(abs(np.linalg.eigvals(M.astype(float)))))


def primitive(M):
    """Wielandt: a primitive n x n matrix has M^k > 0 for k = n^2-2n+2."""
    n = M.shape[0]
    B = (M > 0).astype(np.int64)
    X = B.copy()
    for _ in range(n * n - 2 * n + 2):
        if (X > 0).all():
            return True
        X = ((X @ B) > 0).astype(np.int64)
    return bool((X > 0).all())


def strongly_connected(adj):
    """adj[i][j] true means an edge j -> i (matrix convention X_ij > 0)."""
    n = len(adj)
    for start in (0,):
        for direction in (0, 1):
            seen = {start}
            stack = [start]
            while stack:
                u = stack.pop()
                for v in range(n):
                    e = adj[v][u] if direction == 0 else adj[u][v]
                    if e and v not in seen:
                        seen.add(v)
                        stack.append(v)
            if len(seen) < n:
                return False
    return True


# ---------------------------------------------------------------------------
# Generator sets

class Gen:
    """F = P + e_ab, stored also in the factored form F = P (I + e_{a'b})
    with a' = p^{-1}(a), which is what the frame rewriting uses."""

    def __init__(self, name, p, a, b):
        self.name = name
        self.p = tuple(p)
        n = len(p)
        self.P = perm_matrix(p)
        self.F = self.P + unit(n, a, b)
        pinv = [0] * n
        for c in range(n):
            pinv[p[c]] = c
        self.ep = (pinv[a], b)         # position of e' in I + e'
        assert (self.P @ (np.eye(n, dtype=np.int64) + unit(n, *self.ep))
                == self.F).all()
        assert self.F.max() <= 1 and self.F.sum() == n + 1, \
            "the extra 1 must not land on the permutation"

    @property
    def identity_perm(self):
        return self.p == tuple(range(len(self.p)))


def gens_RL():
    # Dimension 2, exactly the n = 3 automaton: L = I + e_10, R = I + e_01.
    return [Gen("L", (0, 1), 1, 0), Gen("R", (0, 1), 0, 1)]


def gens_rauzy3():
    # Dimension 3, all six I + e_ab.
    return [Gen("e%d%d" % (a, b), (0, 1, 2), a, b)
            for a in range(3) for b in range(3) if a != b]


def gens_perm3():
    # Dimension 3 with permutations: two shears, and a shear followed by
    # the cyclic permutation 0 -> 1 -> 2 -> 0.
    cyc = (1, 2, 0)
    return [Gen("e01", (0, 1, 2), 0, 1),
            Gen("e12", (0, 1, 2), 1, 2),
            Gen("c.e20", cyc, 1, 2)]    # P + e_{12}; e' lands at (0,2)


def group_closure(gens):
    n = len(gens[0].p)
    ident = tuple(range(n))
    H = {ident}
    frontier = [ident]
    while frontier:
        q = frontier.pop()
        for g in gens:
            r = tuple(g.p[q[c]] for c in range(n))   # p o q
            if r not in H:
                H.add(r)
                frontier.append(r)
    return sorted(H)


# ---------------------------------------------------------------------------
# Bounds

def hs_norm_bound(lam, n):
    return math.floor(lam ** n) + n - 1


def c0_prune(A, lam, n):
    if A.sum() > hs_norm_bound(lam, n):
        return True
    if A.sum(axis=0).min() > lam:
        return True
    if A.sum(axis=1).min() > lam:
        return True
    return False


def closure(S, n):
    """Reflexive-transitive closure of the edges b -> a of the e_ab in S,
    as a boolean matrix R with R[i][k] true when k ->* i."""
    R = [[i == k for k in range(n)] for i in range(n)]
    for (a, b) in S:
        R[a][b] = True
    for m in range(n):
        for i in range(n):
            if R[i][m]:
                for k in range(n):
                    if R[m][k]:
                        R[i][k] = True
    return R


def beta(Q, N, positions, H, n):
    """The completion bound: a lower bound on rho(M) for every completion
    M = Q' (I+e~)...(I+e~) N of the prefix Q N whose matrix is
    irreducible.  Q' ranges over the group H (all total permutations a
    completion can produce), S over subsets of the positions the
    conjugated units can take.  Returns +inf when no (Q', S) allows an
    irreducible product."""
    best = math.inf
    Nsupp = N > 0
    I = np.eye(n, dtype=np.int64)
    for r in range(len(positions) + 1):
        for S in itertools.combinations(positions, r):
            R = closure(S, n)
            ES = sum((unit(n, a, b) for (a, b) in S), np.zeros((n, n), dtype=np.int64))
            X = (I + ES) @ N
            for q in H:
                # Edge j -> q(i) whenever j ->N k ->*S i.
                adj = [[False] * n for _ in range(n)]
                for j in range(n):
                    for k in range(n):
                        if Nsupp[k, j]:
                            for i in range(n):
                                if R[i][k]:
                                    adj[q[i]][j] = True
                if not strongly_connected(adj):
                    continue
                best = min(best, rho(perm_matrix(q) @ X))
    return best


# ---------------------------------------------------------------------------
# Search

class Result:
    def __init__(self):
        self.accepted = set()
        self.visited = 0
        self.maxlen = 0
        self.wasted = 0       # visited prefixes whose subtree accepts nothing


def search(gens, lam, criterion):
    n = gens[0].F.shape[0]
    H = group_closure(gens)
    # Every position a conjugated unit q^{-1} e' q can take.
    positions = sorted({(qinv[a], qinv[b])
                        for q in H
                        for qinv in [tuple(sorted(range(n), key=lambda c: q[c]))]
                        for (a, b) in [g.ep for g in gens]})
    res = Result()
    eps = 1e-9 * lam

    def prune(A, Q, N):
        if criterion == "HS":
            return A.sum() > hs_norm_bound(lam, n)
        if c0_prune(A, lam, n):
            return True
        if criterion in ("C1", "BAD") and rho(A) > lam + eps:
            return True
        if criterion == "C2" and beta(Q, N, positions, H, n) > lam + eps:
            return True
        return False

    def walk(word, A, Q, N):
        res.visited += 1
        res.maxlen = max(res.maxlen, len(word))
        found = False
        if word:
            r = rho(A)
            if 1 + 1e-9 < r <= lam + eps and primitive(A):
                res.accepted.add(word)
                found = True
        for g in gens:
            A2 = g.F @ A
            # Frame rewriting: F Q N = (P Q)(I + Q^{-1} e' Q) N.
            q = Q
            qinv = tuple(sorted(range(n), key=lambda c: q[c]))
            a, b = g.ep
            e2 = unit(n, qinv[a], qinv[b])
            Q2 = tuple(g.p[q[c]] for c in range(n))
            N2 = (np.eye(n, dtype=np.int64) + e2) @ N
            assert (perm_matrix(Q2) @ N2 == A2).all()
            if prune(A2, Q2, N2):
                continue
            if walk(word + (g.name,), A2, Q2, N2):
                found = True
        if not found:
            res.wasted += 1
        return found

    ident = tuple(range(n))
    walk((), np.eye(n, dtype=np.int64), ident, np.eye(n, dtype=np.int64))
    return res


def cyclic_classes(words):
    return {min(w[i:] + w[:i] for i in range(len(w))) for w in words}


# ---------------------------------------------------------------------------

def report(label, gens, lams, criteria):
    print("\n== %s: dimension %d, %d generators, permutation group of order %d"
          % (label, gens[0].F.shape[0], len(gens), len(group_closure(gens))))
    rows = []
    for lam in lams:
        truth = search(gens, lam, "HS")
        line = {"lam": lam, "HS": truth}
        for c in criteria:
            line[c] = search(gens, lam, c)
        rows.append(line)
        print("Lambda %-6g accepted %5d words (%d cyclic classes), longest %d"
              % (lam, len(truth.accepted), len(cyclic_classes(truth.accepted)),
                 max((len(w) for w in truth.accepted), default=0)))
        for c in ["HS"] + criteria:
            r = line[c]
            ok = "finds all" if r.accepted == truth.accepted else \
                 "MISSES %d" % len(truth.accepted - r.accepted)
            print("   %-4s visited %8d  wasted %8d  longest prefix %5d  %s"
                  % (c, r.visited, r.wasted, r.maxlen, ok))
    return rows


if __name__ == "__main__":
    report("R, L (the n = 3 automaton)", gens_RL(), [3, 5, 8, 12],
           ["C0", "C1", "C2"])
    report("all I + e_ab", gens_rauzy3(), [2.2, 2.6, 3.0],
           ["C0", "C1", "C2"])
    report("shears and a cyclic permutation", gens_perm3(), [2.0, 2.5, 3.0],
           ["C0", "C2", "BAD"])
