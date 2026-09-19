# ttauto Strata Scan (n=3..7)

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
| ~~3~~ | ~~1. 1. 1. 1. 1. 1. 4 (2)~~ | ~~22~~ | ~~5~~ | ~~1.61803~~ | ~~4~~ | ~~4~~ | ~~5~~ |
| 4 | 1. 1. 1. 1. 1. 1. 5 (1) | 3 | 0 | 2.08102 | 3 | 3 | 4 |
| 5 | 1. 1. 1. 1. 1. 1. 3 3 (2) | 90 | 20 | 2.08102 | 7 | **4** | 8 |
| ~~5~~ | ~~1. 1. 1. 1. 1. 1. 3 3 (2)~~ | ~~90~~ | ~~20~~ | ~~2.01536~~ | ~~6~~ | ~~**4**~~ | ~~6~~ |
| 6 | 1. 1. 1. 1. 1. 1. 3 4 (1) | 21 | 0 | 1.88320 | 4 | 4 | 4 |
| 7 | 1. 1. 1. 1. 1. 1. 3 3 3 (1) | 28 | 0 | 2.17113 | 8 | 8 | 8 |

> Struck-through rows: values before the Bestvina-Handel gate test was
> added to the search (2026-09-19, issue #2).  The 2.01536 class of
> stratum 5 is a 5-puncture pseudo-Anosov plus an idle puncture, hence
> reducible.  The 1.61803 class of stratum 3 (and a 1.93185 class in the
> same stratum) has a characteristic polynomial in x^2, so its transition
> matrix is irreducible but not primitive (eigenvalues +-lambda); such a
> polynomial has no Perron root and was never a candidate in the
> Lefschetz-based enumeration of `devel/iss002/braids.tex` (Lanneau and
> Thiffeault, "On the minimum dilatation of braids on the punctured disc",
> Geom. Dedicata 152, 2011), whereas ttauto's
> irreducibility test admitted it.  The gate test rejects all of these.
> The search length for stratum 5 was raised from 6 to 8 at the same time;
> at length 8 the gate test also rejects a class with dilatation 2.54205
> whose paths extend the 2.01536 cycle.
>
> With these corrections every row of this table agrees with the
> per-stratum table of `braids.tex` for six punctures (s_1..s_7: 1.88320,
> 1.83929, 1.88320, 2.08102, 2.08102, 1.88320, 2.17113).  See
> `devel/iss002/`.

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
| ~~7~~ | ~~1. 1. 1. 1. 1. 1. 1. 3 4 (2)~~ | ~~231~~ | ~~45~~ | ~~1.88320~~ | ~~6~~ | ~~6~~ | ~~7~~ |
| 8 | 1. 1. 1. 1. 1. 1. 1. 3 5 (1) | 24 | 0 | 1.80979 | 4 | 4 | 4 |
| 9 | 1. 1. 1. 1. 1. 1. 1. 4 4 (1) | 12 | 0 | 1.75488 | 4 | 4 | 4 |
| 10 | 1. 1. 1. 1. 1. 1. 1. 3 3 3 (2) | 393 | 75 | 1.61094 | 4 | 4 | 4 |
| 11 | 1. 1. 1. 1. 1. 1. 1. 3 3 4 (1) | 108 | 0 | 2.04249 | 8 | 8 | 8 |
| 12 | 1. 1. 1. 1. 1. 1. 1. 3 3 3 3 (1) | 90 | 0 | 2.02598 | 10 | **8** | 10 |
| ~~12~~ | ~~1. 1. 1. 1. 1. 1. 1. 3 3 3 3 (1)~~ | ~~90~~ | ~~0~~ | ~~2.21497~~ | ~~8~~ | ~~8~~ | ~~8~~ |

> Struck-through rows: stratum 7, value before the gate test (2026-09-19,
> issue #2): the 1.88320 class is the n=6 stratum (4) minimum plus an idle
> puncture, hence reducible, see `devel/iss002/`.  Stratum 12, value at
> search length 8: the length was raised to 10 on 2026-09-19 after Lizi Guo
> found the lower-dilatation pseudo-Anosov 2.02598 at length 10; the gate
> test rejects nothing on this stratum.
>
> Comparison with the seven-puncture table of `devel/iss002/braids.tex`:
> every row agrees (s_1..s_11: 1.55603, 1.46557, 1.46557, 1.55603, 2.04249,
> 1.61094, 2.47541, 1.80979, 1.75488, 1.61094, 2.04249) except s_12, where
> the paper's starred value 2.21497 came from the automaton at search
> length 8 and is superseded by 2.02598 at length 10.  That braid should be
> confirmed with Trains before the paper's table is corrected.
