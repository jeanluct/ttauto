# Independent n=3 check.  Hyperbolic conjugacy classes of PSL(2,Z) with
# trace <= T are the cyclic words in R, L containing both letters.  Each
# braid from ttauto is mapped to PSL(2,Z) by s1 -> R, s2 -> L^-1, brought
# to a conjugate with nonnegative entries, factored into R and L, and its
# cyclic word compared with the enumeration.
import itertools, math, sys
from collections import Counter

R = ((1,1),(0,1)); L = ((1,0),(1,1))
Ri = ((1,-1),(0,1)); Li = ((1,0),(-1,1))
def mul(A,B):
    return ((A[0][0]*B[0][0]+A[0][1]*B[1][0], A[0][0]*B[0][1]+A[0][1]*B[1][1]),
            (A[1][0]*B[0][0]+A[1][1]*B[1][0], A[1][0]*B[0][1]+A[1][1]*B[1][1]))
def tr(A): return A[0][0]+A[1][1]
def neg(A): return tuple(tuple(-x for x in r) for r in A)
def least_rotation(w):
    return min(w[i:]+w[:i] for i in range(len(w)))

def enumerate_classes(T):
    # Depth-first over words in R and L.  Appending R or L to a matrix
    # with nonnegative entries never lowers its trace, so a prefix whose
    # trace exceeds T can be abandoned.  A pure power R^k or L^k keeps
    # trace 2, so length is capped too: a word of length k containing
    # both letters has trace at least k+1.
    import sys
    sys.setrecursionlimit(10000)
    out = {}
    def walk(w, M):
        if tr(M) > T or len(w) > T-1: return
        if "R" in w and "L" in w:
            c = least_rotation(w)
            if c not in out: out[c] = tr(M)
        walk(w+"R", mul(M,R)); walk(w+"L", mul(M,L))
    walk("R", R); walk("L", L)
    return out

def positive_conjugate(A):
    # Breadth-first over conjugations by R, L and their inverses, within
    # a bound on the entries that is raised until a conjugate with
    # nonnegative entries turns up.  One exists for every hyperbolic
    # class with positive trace.
    from collections import deque
    if tr(A) < 0: A = neg(A)
    gens = [(R,Ri),(Ri,R),(L,Li),(Li,L)]
    size = lambda M: max(abs(x) for r in M for x in r)
    cap = size(A)
    while True:
        seen = {A}; q = deque([A])
        while q:
            M = q.popleft()
            if all(x >= 0 for r in M for x in r): return M
            for g, gi in gens:
                N = mul(mul(g,M),gi)
                if N not in seen and size(N) <= cap:
                    seen.add(N); q.append(N)
        cap *= 2

def factor(M):
    w = ""
    while M != ((1,0),(0,1)):
        (a,b),(c,d) = M
        if a >= c and b >= d: w += "R"; M = mul(Ri,M)
        elif c >= a and d >= b: w += "L"; M = mul(Li,M)
        else: raise RuntimeError("not in the positive monoid: %s" % (M,))
    return w

def braid_class(word, s2=Li):
    M = ((1,0),(0,1))
    for g in word:
        M = mul(M, {1:R, -1:Ri, 2:s2, -2:(L if s2 == Li else Li)}[g])
    if abs(tr(M)) <= 2: return "not hyperbolic, trace %d" % tr(M), abs(tr(M))
    return least_rotation(factor(positive_conjugate(M))), abs(tr(M))

T = int(sys.argv[1]); corrupt = len(sys.argv) > 3
expected = enumerate_classes(T)
found = Counter(); bad_trace = 0
for line in open(sys.argv[2]):
    if line.startswith("T="): continue
    lam, _, path, braid, ver = line.rstrip("\n").split("|")
    word = [int(x) for x in braid.split()]
    c, t = braid_class(word, L if corrupt else Li)
    found[c] += 1
    lam = float(lam)
    if abs(lam + 1/lam - t) > 1e-6*t: bad_trace += 1
missing = sorted(set(expected) - set(found), key=lambda c: (expected[c], c))
extra = sorted(set(found) - set(expected))
dup = {c: n for c, n in found.items() if n > 1}
per_trace = Counter(expected.values())
print("T=%d: %d classes expected, %d found, %d missing, %d extra, %d found more than once"
      % (T, len(expected), len(found), len(missing), len(extra), len(dup)))
print("dilatation disagreeing with the braid's trace: %d" % bad_trace)
print("classes per trace:", dict(sorted(per_trace.items())))
if missing: print("missing (first 10):", missing[:10])
if extra: print("extra (first 10):", extra[:10])
if dup: print("duplicated (first 10):", list(dup.items())[:10])
