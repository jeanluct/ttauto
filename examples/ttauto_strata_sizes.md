# Automaton size by stratum (n = 3..9)

> Structural statistics for every stratum: the size of the main
> subautomaton and of the part outside it.  `ttauto::subgraphs` splits
> an automaton into invariant subgraphs by the Dulmage-Mendelsohn
> decomposition; the largest is the main one, and no arrow leads from
> any of the others back to it.
>
> `non-main` counts the vertices outside the main subautomaton and
> `subg` how many invariant subgraphs they form; `biggest` is the
> largest of those.  `pA` is how many pseudo-Anosov classes were found
> anywhere outside the main subautomaton, searching every such
> subgraph to the path length in the section heading, with the
> Bestvina-Handel gate test both on and off.  It is zero everywhere,
> which is what lets the matching column of the paper's tables be
> called `non-pseudo-Anosov` rather than merely `non-main`.
>
> The search length is 12 up to seven punctures and 8 beyond, because
> the cost grows by roughly a factor of five per unit of length: at
> n = 8 stratum 1 it is 4 s at length 8, 17 s at 9 and 88 s at 10.
>
> Minimum dilatations are not here.  They come from a far more
> expensive search and are in `ttauto_scan_strata.md`, for n <= 7 only.
>
> The number of strata is the sum of p(j) for j = 0 to n-3, with p the
> partition function: 1, 2, 4, 7, 12, 19, 30 for n = 3 to 9.

## 3 Punctures (search length <= 12)

| Stratum | Singularity Data | Main | non-main | subg | biggest | pA |
|---:|---|---:|---:|---:|---:|---:|
| 1 | 1. 1. 1. (1) | 1 | 0 | 0 | 0 | 0 |

Total: 1 main + 0 non-main = 1 train tracks.

## 4 Punctures (search length <= 12)

| Stratum | Singularity Data | Main | non-main | subg | biggest | pA |
|---:|---|---:|---:|---:|---:|---:|
| 1 | 1. 1. 1. 1. (2) | 3 | 1 | 1 | 1 | 0 |
| 2 | 1. 1. 1. 1. 3 (1) | 3 | 0 | 0 | 0 | 0 |

Total: 6 main + 1 non-main = 7 train tracks.

## 5 Punctures (search length <= 12)

| Stratum | Singularity Data | Main | non-main | subg | biggest | pA |
|---:|---|---:|---:|---:|---:|---:|
| 1 | 1. 1. 1. 1. 1. (3) | 10 | 1 | 1 | 1 | 0 |
| 2 | 1. 1. 1. 1. 1. 3 (2) | 16 | 5 | 3 | 2 | 0 |
| 3 | 1. 1. 1. 1. 1. 4 (1) | 3 | 0 | 0 | 0 | 0 |
| 4 | 1. 1. 1. 1. 1. 3 3 (1) | 9 | 0 | 0 | 0 | 0 |

Total: 38 main + 6 non-main = 44 train tracks.

## 6 Punctures (search length <= 12)

| Stratum | Singularity Data | Main | non-main | subg | biggest | pA |
|---:|---|---:|---:|---:|---:|---:|
| 1 | 1. 1. 1. 1. 1. 1. (4) | 43 | 6 | 4 | 2 | 0 |
| 2 | 1. 1. 1. 1. 1. 1. 3 (3) | 131 | 7 | 5 | 2 | 0 |
| 3 | 1. 1. 1. 1. 1. 1. 4 (2) | 22 | 5 | 3 | 2 | 0 |
| 4 | 1. 1. 1. 1. 1. 1. 5 (1) | 3 | 0 | 0 | 0 | 0 |
| 5 | 1. 1. 1. 1. 1. 1. 3 3 (2) | 90 | 20 | 7 | 5 | 0 |
| 6 | 1. 1. 1. 1. 1. 1. 3 4 (1) | 21 | 0 | 0 | 0 | 0 |
| 7 | 1. 1. 1. 1. 1. 1. 3 3 3 (1) | 28 | 0 | 0 | 0 | 0 |

Total: 338 main + 38 non-main = 376 train tracks.

## 7 Punctures (search length <= 12)

| Stratum | Singularity Data | Main | non-main | subg | biggest | pA |
|---:|---|---:|---:|---:|---:|---:|
| 1 | 1. 1. 1. 1. 1. 1. 1. (5) | 192 | 12 | 8 | 3 | 0 |
| 2 | 1. 1. 1. 1. 1. 1. 1. 3 (4) | 793 | 62 | 32 | 7 | 0 |
| 3 | 1. 1. 1. 1. 1. 1. 1. 4 (3) | 183 | 7 | 5 | 2 | 0 |
| 4 | 1. 1. 1. 1. 1. 1. 1. 5 (2) | 25 | 5 | 3 | 2 | 0 |
| 5 | 1. 1. 1. 1. 1. 1. 1. 6 (1) | 3 | 0 | 0 | 0 | 0 |
| 6 | 1. 1. 1. 1. 1. 1. 1. 3 3 (3) | 977 | 35 | 18 | 5 | 0 |
| 7 | 1. 1. 1. 1. 1. 1. 1. 3 4 (2) | 231 | 45 | 13 | 12 | 0 |
| 8 | 1. 1. 1. 1. 1. 1. 1. 3 5 (1) | 24 | 0 | 0 | 0 | 0 |
| 9 | 1. 1. 1. 1. 1. 1. 1. 4 4 (1) | 12 | 0 | 0 | 0 | 0 |
| 10 | 1. 1. 1. 1. 1. 1. 1. 3 3 3 (2) | 393 | 75 | 16 | 14 | 0 |
| 11 | 1. 1. 1. 1. 1. 1. 1. 3 3 4 (1) | 108 | 0 | 0 | 0 | 0 |
| 12 | 1. 1. 1. 1. 1. 1. 1. 3 3 3 3 (1) | 90 | 0 | 0 | 0 | 0 |

Total: 3031 main + 241 non-main = 3272 train tracks.

## 8 Punctures (search length <= 8)

| Stratum | Singularity Data | Main | non-main | subg | biggest | pA |
|---:|---|---:|---:|---:|---:|---:|
| 1 | 1. 1. 1. 1. 1. 1. 1. 1. (6) | 932 | 52 | 28 | 9 | 0 |
| 2 | 1. 1. 1. 1. 1. 1. 1. 1. 3 (5) | 5104 | 209 | 115 | 14 | 0 |
| 3 | 1. 1. 1. 1. 1. 1. 1. 1. 4 (4) | 1261 | 77 | 41 | 10 | 0 |
| 4 | 1. 1. 1. 1. 1. 1. 1. 1. 5 (3) | 246 | 7 | 5 | 2 | 0 |
| 5 | 1. 1. 1. 1. 1. 1. 1. 1. 6 (2) | 31 | 5 | 3 | 2 | 0 |
| 6 | 1. 1. 1. 1. 1. 1. 1. 1. 7 (1) | 3 | 0 | 0 | 0 | 0 |
| 7 | 1. 1. 1. 1. 1. 1. 1. 1. 3 3 (4) | 8814 | 0 | 0 | 0 | 0 |
| 8 | 1. 1. 1. 1. 1. 1. 1. 1. 3 4 (3) | 2848 | 77 | 37 | 12 | 0 |
| 9 | 1. 1. 1. 1. 1. 1. 1. 1. 3 5 (2) | 301 | 50 | 14 | 14 | 0 |
| 10 | 1. 1. 1. 1. 1. 1. 1. 1. 3 6 (1) | 27 | 0 | 0 | 0 | 0 |
| 11 | 1. 1. 1. 1. 1. 1. 1. 1. 4 4 (2) | 157 | 25 | 8 | 7 | 0 |
| 12 | 1. 1. 1. 1. 1. 1. 1. 1. 4 5 (1) | 27 | 0 | 0 | 0 | 0 |
| 13 | 1. 1. 1. 1. 1. 1. 1. 1. 3 3 3 (3) | 5936 | 154 | 59 | 14 | 0 |
| 14 | 1. 1. 1. 1. 1. 1. 1. 1. 3 3 4 (2) | 1690 | 275 | 54 | 24 | 0 |
| 15 | 1. 1. 1. 1. 1. 1. 1. 1. 3 3 5 (1) | 135 | 0 | 0 | 0 | 0 |
| 16 | 1. 1. 1. 1. 1. 1. 1. 1. 3 4 4 (1) | 135 | 0 | 0 | 0 | 0 |
| 17 | 1. 1. 1. 1. 1. 1. 1. 1. 3 3 3 3 (2) | 1725 | 275 | 45 | 28 | 0 |
| 18 | 1. 1. 1. 1. 1. 1. 1. 1. 3 3 3 4 (1) | 495 | 0 | 0 | 0 | 0 |
| 19 | 1. 1. 1. 1. 1. 1. 1. 1. 3 3 3 3 3 (1) | 297 | 0 | 0 | 0 | 0 |

Total: 30164 main + 1206 non-main = 31370 train tracks.

## 9 Punctures (search length <= 8)

| Stratum | Singularity Data | Main | non-main | subg | biggest | pA |
|---:|---|---:|---:|---:|---:|---:|
| 1 | 1. 1. 1. 1. 1. 1. 1. 1. 1. (7) | 4618 | 189 | 84 | 20 | 0 |
| 2 | 1. 1. 1. 1. 1. 1. 1. 1. 1. 3 (6) | 31634 | 1264 | 513 | 48 | 0 |
| 3 | 1. 1. 1. 1. 1. 1. 1. 1. 1. 4 (5) | 8691 | 279 | 150 | 18 | 0 |
| 4 | 1. 1. 1. 1. 1. 1. 1. 1. 1. 5 (4) | 1863 | 87 | 49 | 11 | 0 |
| 5 | 1. 1. 1. 1. 1. 1. 1. 1. 1. 6 (3) | 320 | 7 | 5 | 2 | 0 |
| 6 | 1. 1. 1. 1. 1. 1. 1. 1. 1. 7 (2) | 34 | 5 | 3 | 2 | 0 |
| 7 | 1. 1. 1. 1. 1. 1. 1. 1. 1. 8 (1) | 3 | 0 | 0 | 0 | 0 |
| 8 | 1. 1. 1. 1. 1. 1. 1. 1. 1. 3 3 (5) | 69219 | 2034 | 879 | 55 | 0 |
| 9 | 1. 1. 1. 1. 1. 1. 1. 1. 1. 3 4 (4) | 26239 | 1166 | 494 | 78 | 0 |
| 10 | 1. 1. 1. 1. 1. 1. 1. 1. 1. 3 5 (3) | 3976 | 84 | 40 | 14 | 0 |
| 11 | 1. 1. 1. 1. 1. 1. 1. 1. 1. 3 6 (2) | 380 | 55 | 15 | 16 | 0 |
| 12 | 1. 1. 1. 1. 1. 1. 1. 1. 1. 3 7 (1) | 30 | 0 | 0 | 0 | 0 |
| 13 | 1. 1. 1. 1. 1. 1. 1. 1. 1. 4 4 (3) | 1988 | 42 | 21 | 7 | 0 |
| 14 | 1. 1. 1. 1. 1. 1. 1. 1. 1. 4 5 (2) | 380 | 55 | 15 | 16 | 0 |
| 15 | 1. 1. 1. 1. 1. 1. 1. 1. 1. 4 6 (1) | 30 | 0 | 0 | 0 | 0 |
| 16 | 1. 1. 1. 1. 1. 1. 1. 1. 1. 5 5 (1) | 15 | 0 | 0 | 0 | 0 |
| 17 | 1. 1. 1. 1. 1. 1. 1. 1. 1. 3 3 3 (4) | 65411 | 2789 | 939 | 110 | 0 |
| 18 | 1. 1. 1. 1. 1. 1. 1. 1. 1. 3 3 4 (3) | 26734 | 546 | 201 | 24 | 0 |
| 19 | 1. 1. 1. 1. 1. 1. 1. 1. 1. 3 3 5 (2) | 2310 | 330 | 63 | 28 | 0 |
| 20 | 1. 1. 1. 1. 1. 1. 1. 1. 1. 3 3 6 (1) | 165 | 0 | 0 | 0 | 0 |
| 21 | 1. 1. 1. 1. 1. 1. 1. 1. 1. 3 4 4 (2) | 2310 | 330 | 63 | 24 | 0 |
| 22 | 1. 1. 1. 1. 1. 1. 1. 1. 1. 3 4 5 (1) | 330 | 0 | 0 | 0 | 0 |
| 23 | 1. 1. 1. 1. 1. 1. 1. 1. 1. 4 4 4 (1) | 55 | 0 | 0 | 0 | 0 |
| 24 | 1. 1. 1. 1. 1. 1. 1. 1. 1. 3 3 3 3 (3) | 32096 | 637 | 196 | 28 | 0 |
| 25 | 1. 1. 1. 1. 1. 1. 1. 1. 1. 3 3 3 4 (2) | 10120 | 1430 | 211 | 60 | 0 |
| 26 | 1. 1. 1. 1. 1. 1. 1. 1. 1. 3 3 3 5 (1) | 660 | 0 | 0 | 0 | 0 |
| 27 | 1. 1. 1. 1. 1. 1. 1. 1. 1. 3 3 4 4 (1) | 990 | 0 | 0 | 0 | 0 |
| 28 | 1. 1. 1. 1. 1. 1. 1. 1. 1. 3 3 3 3 3 (2) | 7150 | 1001 | 126 | 70 | 0 |
| 29 | 1. 1. 1. 1. 1. 1. 1. 1. 1. 3 3 3 3 4 (1) | 2145 | 0 | 0 | 0 | 0 |
| 30 | 1. 1. 1. 1. 1. 1. 1. 1. 1. 3 3 3 3 3 3 (1) | 1001 | 0 | 0 | 0 | 0 |

Total: 300897 main + 12330 non-main = 313227 train tracks.

## Totals

| n | strata | main | non-main | total | largest stratum | growth |
|---:|---:|---:|---:|---:|---:|---:|
| 3 | 1 | 1 | 0 | 1 | 1 | - |
| 4 | 2 | 6 | 1 | 7 | 4 | 7.0 |
| 5 | 4 | 38 | 6 | 44 | 21 | 6.3 |
| 6 | 7 | 338 | 38 | 376 | 138 | 8.5 |
| 7 | 12 | 3031 | 241 | 3272 | 1012 | 8.7 |
| 8 | 19 | 30164 | 1206 | 31370 | 8814 | 9.6 |
| 9 | 30 | 300897 | 12330 | 313227 | 71253 | 10.0 |
