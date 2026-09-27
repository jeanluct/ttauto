# n = 3 against PSL(2,Z)

Section 2 of issue #22, done 2026-09-26: a complete enumeration of the
pseudo-Anosov 3-braids to a dilatation bound, checked class for class
against an independent enumeration of the hyperbolic classes of
PSL(2,Z).  Everything matches, to trace 150.

## The two objects

**ttauto.**  `build_traintrack_list(3)` has one stratum, and its automaton
has one vertex (`1111 1112 1122 1111`, no cyclic symmetry) with two
branches, both loops.  Their transition matrices are exactly

    branch 0:  [[1,0],[1,1]] = L        braid  s2
    branch 1:  [[1,1],[0,1]] = R        braid  s1 s2^-1 s1^-1

so a closed path is a word in L and R, and its length is the length of
that word.

**PSL(2,Z).**  B3 modulo its centre is PSL(2,Z), by s1 -> R,
s2 -> L^-1.  A pA 3-braid maps to a hyperbolic class, and
lambda + 1/lambda = |trace|.  Hyperbolic classes with positive trace
correspond one to one with cyclic words in R and L containing both
letters, so they can be enumerated directly.  Nothing from class-number
theory is needed.

## Method

- `n3_search.cpp`: the automaton, `max_dilatation(lambda_T)` with
  lambda_T = (T + sqrt(T^2 - 4))/2, `check_norms()`, bad words off
  (the default), gates on, and `max_paths_to_save` large enough to keep
  every accepted path.  `path_length_exceeded()` is 0 in every run, so
  each list is complete for its window.
- `n3_psl.py`: enumerates the cyclic words with trace <= T, pruning on
  trace (appending R or L to a nonnegative matrix never lowers it) and on
  length (a word of length k with both letters has trace >= k+1).  For
  each ttauto path it takes the *braid*, maps it to PSL(2,Z), finds a
  conjugate with nonnegative entries by breadth-first search, factors it
  into R and L, and takes the least rotation.  The two sides share
  nothing but the braid word.
- The comparison checks: every class found, none extra, none twice, and
  lambda + 1/lambda equal to the braid's |trace| for every path.

Build and run from the repository root, with the CMake build done:

    g++ -std=c++17 -O2 -Itestsuite -Iinclude -Iextern/jlt \
      -Iextern/jlt/extern/CSparse/Include devel/iss022/n3_search.cpp \
      lib/libttauto.a build/libcsparse.a -lm -o /tmp/n3_search
    /tmp/n3_search 150 >| /tmp/t150.txt
    python3 devel/iss022/n3_psl.py 150 /tmp/t150.txt
    python3 devel/iss022/n3_growth.py

## Results

| T | classes | found by ttauto | missing | extra | twice | ttauto time |
|---|---|---|---|---|---|---|
| 10 | 23 | 23 | 0 | 0 | 0 | 0.02 s |
| 20 | 78 | 78 | 0 | 0 | 0 | 0.09 s |
| 40 | 254 | 254 | 0 | 0 | 0 | 0.45 s |
| 60 | 502 | 502 | 0 | 0 | 0 | 0.59 s |
| 100 | 1240 | 1240 | 0 | 0 | 0 | 2.4 s |
| 150 | 2541 | 2541 | 0 | 0 | 0 | 6.9 s |

The gate test rejected nothing.  Every braid was verified (its Dynnikov
growth equals the path's Perron root).

**The check is not vacuous.**  With s2 -> L instead of L^-1, at T = 10:
18 classes missing, 15 extra, 3 found twice, and 20 dilatations
disagreeing with the braid's trace.

**Per trace**, the counts are those of the class numbers, e.g. trace 4,
discriminant 12, has 2 classes; trace 6, discriminant 32 = 8 * 2^2, has
h+(32) + h+(8) = 2 + 1 = 3.  (Spot checks only; the enumeration does not
use them.)

## Deduplication: what ttauto's own count says

`pA_list()` groups by characteristic polynomial, and a 2x2 polynomial is
x^2 - t x + 1, fixed by the trace.  So ttauto reports one "class" per
trace: 148 at T = 150, against 2541 conjugacy classes.  The paths are
all there (`max_paths_to_save`), but any count taken from `pA_list()`
alone undercounts, badly.  This is pitfall 1 of #22, now measured.  At
n = 3 the least rotation of the branch word is already a complete
conjugacy invariant; for larger n it will not be.

## Growth, against the counting theorems

Primitive classes are the words that are not proper powers.

| T | lambda | classes | primitive | lambda^2/(2 log lambda) | ratio | li(lambda^2) | ratio |
|---|---|---|---|---|---|---|---|
| 10 | 9.90 | 23 | 22 | 21.4 | 1.029 | 29.7 | 0.741 |
| 20 | 19.95 | 78 | 74 | 66.5 | 1.113 | 85.1 | 0.870 |
| 40 | 39.97 | 254 | 245 | 216.6 | 1.131 | 261.2 | 0.938 |
| 80 | 79.99 | 844 | 824 | 730.1 | 1.129 | 846.0 | 0.974 |
| 150 | 149.99 | 2541 | 2505 | 2245.0 | 1.116 | 2539.3 | 0.986 |

The primitive count tends to li(lambda^2): the prime geodesic theorem for
the modular group, with x = lambda^2 the norm (Sarnak 1982).  The
Eskin–Mirzakhani form e^{hR}/(hR), with h = 2n - 4 = 2 and R = log lambda,
is its leading term, and the ratio to it sits near 1.12 because of the
1/log x correction that li carries.  So at n = 3 the certified counts show
the predicted exponent.

## What this settles, and what it does not

- At n = 3 the pipeline (automaton, primitivity, gate test, braid
  extraction, dilatation) is right, class for class, to trace 150.
- n = 3 is special: one vertex, no symmetry, no rejected paths, and a
  free monoid of paths.  None of the hard parts of larger n (several
  vertices, cyclic symmetry, gate rejections, strata) is exercised here.
- Path length is word length, so the shortest path of a class is the
  length of its cyclic word, the sum of the partial quotients of its
  continued-fraction period.  That is the n = 3 case of #22 section 4;
  Kawamuro–Kin's Agol cycles for 3-braids would be the comparison.
