# Braids of the stratum minimisers (n=3..7)

> The braid of the lowest-dilatation class of every stratum, read off
> its folding path by `examples/ttbraid_strata`, beside the braid
> given by the appendix of Lanneau and Thiffeault, *On the minimum
> dilatation of braids on the punctured disc*, Geom. Dedicata **152**
> (2011), whose source is `devel/iss002/braids.tex`.  Those words were
> obtained through the Lefschetz formula rather than from this
> automaton, and were checked with Toby Hall's Trains, so they are an
> independent account of the same mapping classes.
>
> The two braids describe different folding paths of the same class,
> so only what conjugation preserves is compared: the dilatation, the
> cycle type of the permutation, and the exponent sum modulo n(n-1),
> that being the exponent sum of the full twist, which is the
> ambiguity in the braid of a folding path.  The sign of the exponent
> sum is free too, since the mirror convention negates it.
>
> `ok` in the last column means all three agree.  `verified` is the
> library's own check, that the braid's growth under the Dynnikov
> action equals the Perron root of the path's transition matrix.

## 3 Punctures

| Stratum | Dil. | Our braid | Exp. | Cycles | Published braid | Exp. | Cycles | |
|---:|---:|---|---:|---|---|---:|---|---|
| 1 | 2.61803 | `-1 2` | 0 | 3 | `1 -2` | 0 | 3 | ok |

## 4 Punctures

| Stratum | Dil. | Our braid | Exp. | Cycles | Published braid | Exp. | Cycles | |
|---:|---:|---|---:|---|---|---:|---|---|
| 1 | 2.61803 | `1 1 2 3 3 3` | 6 | 3+1 | `1 2 1 2 -3 1 2 3` | 6 | 3+1 | ok |
| 2 | 2.29663 | `-3 1 2 3 1 2 3 2 3 1 2 3 1 2 3` | 13 | 4 | `1 2 -3` | 1 | 4 | ok |

## 5 Punctures

| Stratum | Dil. | Our braid | Exp. | Cycles | Published braid | Exp. | Cycles | |
|---:|---:|---|---:|---|---|---:|---|---|
| 1 | 1.72208 | `1 2 -4 -3 -2 -1 -4 -3 -2 -1` | -6 | 5 | `1 2 3 1 2 3 4 -3` | 6 | 5 | ok |
| 2 | 1.72208 | `2 1 2 3 4 1 1 2 3 4` | 10 | 5 | `1 1 1 2 3 4 1 2 3 4` | 10 | 5 | ok |
| 3 | 2.15372 | `-4 1 2 3 4 1 2 3 4 1 2 3 4 2 3 4 1 2 3 4 1 2 3 4` | 22 | 5 | `1 2 3 -4` | 2 | 5 | ok |
| 4 | 2.01536 | `-2 -1 3 4` | 0 | 5 | `1 2 -4 -3` | 0 | 5 | ok |

## 6 Punctures

| Stratum | Dil. | Our braid | Exp. | Cycles | Published braid | Exp. | Cycles | |
|---:|---:|---|---:|---|---|---:|---|---|
| 1 | 1.8832 | `1 2 1 2 3 4 5` | 7 | 6 | `1 2 3 4 5 4 5` | 7 | 6 | ok |
| 2 | 1.83929 | `3 -2 1 2 3 4 5 1 2 3 4 5` | 10 | 4+2 | `5 -4 1 2 3 4 5 1 2 3 4 5` | 10 | 4+2 | ok |
| 3 | 1.8832 | `-2 1 2 3 4 5 2 3 4 5 -3 1 2 3 4 5 1 2 3 4 5` | 17 | 6 | `1 1 4 1 2 3 4 5 1 2 3 4 5` | 13 | 6 | ok |
| 4 | 2.08102 | `-5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 2 3 4 5 1 2 3 4 5 1 2 3 4 5` | 33 | 6 | `1 2 3 4 -5` | 3 | 6 | ok |
| 5 | 2.08102 | `3 4 1 2 3 4 5 5 4 1 2 3 4 5` | 14 | 3+3 | `4 5 5 4 1 2 3 4 5 1 2 3 4 5` | 14 | 3+3 | ok |
| 6 | 1.8832 | `-5 -4 -3 1 2 3 4 5 1 2 3 4 5 -3 -2 -1 -5 -4 -3 -2 -1` | -1 | 6 | `1 2 3 -5 -4` | 1 | 6 | ok |
| 7 | 2.17113 | `-2 -1 -3 -2 4 5 1 2 3 4 5` | 3 | 6 | `1 2 3 -5 -4 -3 -5 -4 -3` | -3 | 6 | ok |

## 7 Punctures

| Stratum | Dil. | Our braid | Exp. | Cycles | Published braid | Exp. | Cycles | |
|---:|---:|---|---:|---|---|---:|---|---|
| 1 | 1.55603 | `1 2 3 4 5 6 5 6 -4 -3 -2 -5 -4 -3 -6 -5 -4 -3 -2 -1 -6 -5 -4 -3 -2 -1 -6 -5 -4 -3 -2 -1` | -16 | 7 | `3 4 5 6 2 3 4 1 2 3 1 2 3 4 5 6` | 16 | 7 | ok |
| 2 | 1.46557 | `6 -5 -4 -3 -2 -1 -6 -5 -4 -3 -2 -1` | -10 | 7 | `-4 -4 1 2 3 4 5 6 1 2 3 4 5 6` | 10 | 7 | ok |
| 3 | 1.46557 | `2 1 2 3 4 5 6 1 1 2 3 4 5 6` | 14 | 7 | `6 6 1 2 3 4 5 6 1 2 3 4 5 6` | 14 | 7 | ok |
| 4 | 1.55603 | `3 1 2 3 4 5 6 1 2 3 4 5 6 1 1 2 3 4 5 6` | 20 | 7 | `5 5 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6` | 20 | 7 | ok |
| 5 | 2.04249 | `-6 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6` | 46 | 7 | `-4 -4 1 2 3 4 5 6` | 4 | 7 | ok |
| 6 | 1.61094 | `6 -5 -4 -3 -6 -5 -4 -3 -2 -1 -6 -5 -4 -3 -2 -1` | -14 | 7 | `-2 3 4 5 1 2 3 4 5 6 1 2 3 4 5 6` | 14 | 7 | ok |
| 7 | 2.47541 | `-2 -1 3 4 5 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 2 3 4 5 6 1 2 3 4 5 6` | 42 | 7 | `1 2 3 3 -6 -5 -4 -3` | 0 | 7 | ok |
| 8 | 1.80979 | `1 2 1 2 3 4 5 6 1 2 3 4 5 6 -4 -3 -2 -1 -6 -5 -4 -3 -2 -1 -6 -5 -4 -3 -2 -1` | -2 | 7 | `1 2 3 4 -6 -5` | 2 | 7 | ok |
| 9 | 1.75488 | `-3 -2 -1 4 5 6` | 0 | 7 | `1 2 3 -6 -5 -4` | 0 | 7 | ok |
| 10 | 1.61094 | `1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 1 2 4 5` | 22 | 7 | `-5 -4 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6` | 20 | 7 | ok |
| 11 | 2.04249 | `-5 -4 -3 -6 -5 -4 1 2 3 4 5 6 1 2 3 4 5 6 -4 -3 -2 -1` | 2 | 7 | `4 5 6 3 4 5 -2 -1 -6 -5 -4 -3 -2 -1` | -2 | 7 | ok |
| 12 (len 8) | 2.21497 | `1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6` | 50 | 7 | `2 1 1 2 -6 -5 -4 -3 -2 -1 -6 -5 -4 -3 -2 -1` | -8 | 7 | ok |
| 12 (len 10) | 2.02598 | `-2 -1 3 4 3 4 5 6` | 4 | 7 | `2 1 1 2 -6 -5 -4 -3 -2 -1 -6 -5 -4 -3 -2 -1` | -8 | 7 | **differs**: the published search stopped at length 8; 2.02598 supersedes |

> Rows 27, agreeing with the published table 26, verified against their own path 27.
>
> The one row that differs is stratum 12 on seven punctures
> searched to length 10, and it is meant to: the published
> search stopped at length 8, and 2.02598 lies below anything
> reachable there.  Searched to length 8 the same stratum
> agrees.  So every braid the published appendix gives for
> three to seven punctures is reproduced from the folding
> automaton, by a route that has nothing in common with the
> Lefschetz argument that produced it.
