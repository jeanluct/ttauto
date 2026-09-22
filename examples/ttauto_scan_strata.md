# ttauto Strata Scan (n=3..7)

> Minimum dilatations by stratum, from `ttauto_scan_strata.sh`.  `Main`
> and `non-Main` count the vertices of the main subautomaton and of the
> rest.  `Min. Dil Len.` is the folding-path length of the minimiser and
> `Shortest Len.` that of the shortest pseudo-Anosov, in bold where the
> two differ.  `Max Len.` is the length the search ran to.
>
> A closed folding path is accepted as pseudo-Anosov when its transition
> matrix is primitive and the Bestvina-Handel gates are connected at
> every vertex; see `doc/ttauto.tex`.  Since the search is bounded by
> path length, each value is an upper bound on the minimum of its
> stratum rather than a proof of it.
>
> Automaton sizes alone, without the dilatation search and reaching to
> nine punctures, are in `ttauto_strata_sizes.md`.

## 3 Punctures

| Stratum | Singularity Data | Main | non-Main | Min Dil. | Min. Dil Len. | Shortest Len. | Max Len. |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | 1. 1. 1. (1) | 1 | 0 | 2.61803 | 2 | 2 | 4 |

## 4 Punctures

| Stratum | Singularity Data | Main | non-Main | Min Dil. | Min. Dil Len. | Shortest Len. | Max Len. |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | 1. 1. 1. 1. (2) | 3 | 1 | 2.61803 | 3 | 3 | 4 |
| 2 | 1. 1. 1. 1. 3 (1) | 3 | 0 | 2.29663 | 3 | 3 | 4 |

## 5 Punctures

| Stratum | Singularity Data | Main | non-Main | Min Dil. | Min. Dil Len. | Shortest Len. | Max Len. |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | 1. 1. 1. 1. 1. (3) | 10 | 1 | 1.72208 | 2 | 2 | 4 |
| 2 | 1. 1. 1. 1. 1. 3 (2) | 16 | 5 | 1.72208 | 3 | 3 | 4 |
| 3 | 1. 1. 1. 1. 1. 4 (1) | 3 | 0 | 2.15372 | 3 | 3 | 4 |
| 4 | 1. 1. 1. 1. 1. 3 3 (1) | 9 | 0 | 2.01536 | 4 | 4 | 4 |

## 6 Punctures

| Stratum | Singularity Data | Main | non-Main | Min Dil. | Min. Dil Len. | Shortest Len. | Max Len. |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | 1. 1. 1. 1. 1. 1. (4) | 43 | 6 | 1.88320 | 3 | 3 | 4 |
| 2 | 1. 1. 1. 1. 1. 1. 3 (3) | 131 | 7 | 1.83929 | 4 | 4 | 4 |
| 3 | 1. 1. 1. 1. 1. 1. 4 (2) | 22 | 5 | 1.88320 | 5 | 5 | 5 |
| 4 | 1. 1. 1. 1. 1. 1. 5 (1) | 3 | 0 | 2.08102 | 3 | 3 | 4 |
| 5 | 1. 1. 1. 1. 1. 1. 3 3 (2) | 90 | 20 | 2.08102 | 7 | **4** | 8 |
| 6 | 1. 1. 1. 1. 1. 1. 3 4 (1) | 21 | 0 | 1.88320 | 4 | 4 | 4 |
| 7 | 1. 1. 1. 1. 1. 1. 3 3 3 (1) | 28 | 0 | 2.17113 | 8 | 8 | 8 |

> On strata 3 and 5 the search finds candidates of lower dilatation
> that are not pseudo-Anosov, and rejects them: on stratum 5 a 2.01536
> class, which is the five-puncture 3 3 (1) minimum with an idle
> puncture, and at length 8 a 2.54205 class extending it; on stratum 3
> classes at 1.61803 and 1.93185 whose characteristic polynomials are
> polynomials in x^2, so the matrix is irreducible but not primitive,
> with eigenvalues +-lambda and no Perron root.  Stratum 5 needs the
> longest search here, to length 8, and its minimiser is longer than its
> shortest pseudo-Anosov.
>
> Every row agrees with the per-stratum table of
> `devel/iss002/braids.tex` (Lanneau and Thiffeault, "On the minimum
> dilatation of braids on the punctured disc", Geom. Dedicata 152,
> 2011), which reaches the same values through the Lefschetz formula.
> The imprimitive classes above were never candidates there, since that
> method enumerates polynomials with a Perron root.

## 7 Punctures

| Stratum | Singularity Data | Main | non-Main | Min Dil. | Min. Dil Len. | Shortest Len. | Max Len. |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | 1. 1. 1. 1. 1. 1. 1. (5) | 192 | 12 | 1.55603 | 3 | 3 | 4 |
| 2 | 1. 1. 1. 1. 1. 1. 1. 3 (4) | 793 | 62 | 1.46557 | 3 | 3 | 4 |
| 3 | 1. 1. 1. 1. 1. 1. 1. 4 (3) | 183 | 7 | 1.46557 | 3 | 3 | 4 |
| 4 | 1. 1. 1. 1. 1. 1. 1. 5 (2) | 25 | 5 | 1.55603 | 3 | 3 | 4 |
| 5 | 1. 1. 1. 1. 1. 1. 1. 6 (1) | 3 | 0 | 2.04249 | 3 | 3 | 4 |
| 6 | 1. 1. 1. 1. 1. 1. 1. 3 3 (3) | 977 | 35 | 1.61094 | 5 | 5 | 5 |
| 7 | 1. 1. 1. 1. 1. 1. 1. 3 4 (2) | 231 | 45 | 2.47541 | 7 | **6** | 7 |
| 8 | 1. 1. 1. 1. 1. 1. 1. 3 5 (1) | 24 | 0 | 1.80979 | 4 | 4 | 4 |
| 9 | 1. 1. 1. 1. 1. 1. 1. 4 4 (1) | 12 | 0 | 1.75488 | 4 | 4 | 4 |
| 10 | 1. 1. 1. 1. 1. 1. 1. 3 3 3 (2) | 393 | 75 | 1.61094 | 4 | 4 | 4 |
| 11 | 1. 1. 1. 1. 1. 1. 1. 3 3 4 (1) | 108 | 0 | 2.04249 | 8 | 8 | 8 |
| 12 | 1. 1. 1. 1. 1. 1. 1. 3 3 3 3 (1) | 90 | 0 | 2.02598 | 10 | **8** | 10 |

> On stratum 7 the search rejects a 1.88320 class, the six-puncture (4)
> minimum with an idle puncture, and the minimiser shown is longer than
> the shortest pseudo-Anosov.  Stratum 12 needs the longest search, to
> length 10: its 2.02598 class, found by Lizi Guo, lies below anything
> reachable at length 8, so this stratum in particular should be read as
> an upper bound.  Its braid has now been read off the folding path
> (`examples/ttbraid`, issue #4): the closed path of length 10 with branch
> sequence `0 0 1 0 0 1 1 1 0 0` gives `-2 -1 3 4 3 4 5 6`, which Toby
> Hall's Trains independently calls pseudo-Anosov with dilatation
> 2.02598.  So the value no longer rests on ttauto alone.  The class at
> length 8 on the same stratum gives a braid Trains calls pseudo-Anosov
> with dilatation 2.21497, confirming the starred entry of
> `devel/iss002/braids.tex` as the best reachable at that length.
>
> Every row agrees with the seven-puncture table of
> `devel/iss002/braids.tex` except s_12, where that paper's starred
> 2.21497 comes from an automaton search of length 8 and is superseded
> by 2.02598.
