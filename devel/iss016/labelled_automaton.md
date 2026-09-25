# What labelling does, measured

Behind item 1 of issue #16: the class counts claimed by
`examples/ttauto_labels.cpp`.  Measured 2026-09-25 on this machine with
the CMake Release build.  Class counts and branch counts are
deterministic; timings are not.

## Fixture

Five punctures, `build_traintrack_list(5)[3]`, stratum
`1. 1. 1. 1. 1. 3 3 (1)`, six edges.  One subgraph.  Multigon 6 is one of
the two trigons; `set_label(6,1)` distinguishes it from the other.

| | vertices |
|---|---:|
| unlabelled | 9 |
| labelled | 18 |

Forgetting the label projects the labelled fold graph onto the unlabelled
one.  The **fibre** over an unlabelled vertex is the set of labelled
vertices above it, that is, the distinct labellings of that track up to
isomorphism.

Here there are two labellings -- the label sits on one trigon or the
other -- and the fibre over a track `T` is the set of orbits of `Aut(T)`
on them: 1 if some automorphism of `T` exchanges the trigons, 2 if none
does.  Since the totals are 18 and 9 and no fibre can exceed 2, every
fibre is exactly 2, and so **no track in this automaton has an
automorphism exchanging its two trigons**.  That is derived from the
counts rather than assumed.

Note this is much smaller than the `n!` of a full `pure_braid()`
labelling: only one label is set, on one of two multigons.

## The comparison

Window `max_dilatation(3).check_norms()`, `badword_length(0)`, gate test
on, path-length cap 12.

**10 classes without the label, 6 with**, so 4 of the 10 have folding
paths that permute the two trigons.  Checked by characteristic
polynomial, the labelled set is a subset of the unlabelled one: no class
appears only when labelled, which is the sanity condition.

Unlabelled, with the four that the label eliminates marked:

| dilatation | shortest | characteristic polynomial | survives |
|---:|---:|---|:-:|
| 2.01536 | 4  | `x^6 - x^5 - 4x^3 - x + 1`              | yes |
| 2.47032 | 7  | `x^6 - x^5 - x^4 - 6x^3 - x^2 - x + 1`  | yes |
| 2.54205 | 6  | `x^6 - x^5 - 2x^4 - 4x^3 - 2x^2 - x + 1`| yes |
| 2.81739 | 8  | `x^6 - x^5 - 2x^4 - 8x^3 - 2x^2 - x + 1`| yes |
| 2.89005 | 9  | `x^6 - 2x^5 - x^4 - 4x^3 - x^2 - 2x + 1`| yes |
| 2.96557 | 10 | `x^6 - 5x^5 + 8x^4 - 8x^3 + 8x^2 - 5x + 1` | yes |
| 2.45317 | 7  | `x^6 - 3x^5 + 2x^4 - 2x^3 + 2x^2 - 3x + 1` | no |
| 2.75146 | 9  | `x^6 - 3x^5 + 2x^4 - 4x^3 + 2x^2 - 3x + 1` | no |
| 2.89005 | 10 | `x^6 - 3x^5 + x^4 - 2x^3 + x^2 - 3x + 1`   | no |
| 2.89005 | 10 | `x^6 - 4x^5 + 3x^4 + 3x^2 - 4x + 1`        | no |

The minimum, 2.01536 at length 4, is the same either way.  It is the
control, not the finding: what labelling changes is which classes close
up, not how low they go.

## Why the example used to take a hundred minutes

`check_norms()` derives its path-length bound from the dilatation window:
`max_path_length = maxnorm = floor(lambda^edges) + edges - 1`, which is
**734** at `lambda = 3` with six edges.  That is not an arbitrary stop.
`check_all_norms` compares the running matrix's total entry sum against
the same `maxnorm`, the sum starts at `edges` and grows by at least one
per fold, so a path longer than `maxnorm - edges` cannot stay under the
bound.  The length cap is the norm bound restated, and searching to it is
what makes the result complete for the window.

So the hundred minutes was not waste.  It was the price of completeness.

## The cap is measured, not proved

Counts by path-length cap, window 3, `badword_length(0)`:

| cap | unlabelled | labelled |
|---:|---:|---:|
| 6  | 2  | 2 |
| 8  | 5  | 4 |
| 10 | 10 | 6 |
| 12 | 10 | 6 |
| 14 | 10 | 6 |
| 16 | 10 | 6 |
| 18 | 10 | 6 |
| 20 | 10 | 6 |
| 22 | 10 | 6 |
| 24 | 10 | 6 |

Stable from 10 onwards.  But `path_length_exceeded()` -- the count of
branches the cap truncated -- **grows** rather than falling to zero:

| cap | branches cut, unlabelled | labelled | time |
|---:|---:|---:|---:|
| 16 | 104516   | 106322   | 0.47 s |
| 18 | 359032   | 362822   | 1.77 s |
| 20 | 1164826  | 1200110  | 5.83 s |
| 22 | 3482710  | 3780492  | 21.6 s |
| 24 | 9468958  | 10861282 | 70.0 s |

A zero there would have proved the cap sufficient, since the norm tests
would then have ended every path on their own.  It is nowhere near zero,
so there is no cheap cap that is provably enough, and the example's cap
of 12 rests on the table above rather than on an argument.  Say so when
quoting its numbers.

## The counts do not depend on pruning or on the gate test

At cap 12, every combination gives the same 10 and 6:

| badword_length | gates | unlabelled | labelled |
|---:|:-:|---:|---:|
| 0 | on  | 10 | 6 |
| 0 | off | 10 | 6 |
| 1 | on  | 10 | 6 |
| 1 | off | 10 | 6 |

## The old 8 and 2

The previous comment claimed 8 and 2.  Nothing tried reproduces it:
neither the gate test, nor bad-word pruning at lengths 0 and 1, nor any
cap from 10 to 24.  `badword_length(2)` -- the default in 2010, when
`check_norms()` still rebuilt the table -- was being tested uncapped when
this note was written; if it turns out to reproduce 8 and 2, that is the
explanation, since pruning removes exactly the non-minimal classes a
count of this kind is counting.

Either way it does not affect what the example is for.  Labelling
eliminating some classes is a *relative* statement, and both sides of the
comparison run the identical configuration, so a prune or a cap removes
from both.  What the episode shows is narrower: absolute counts move with
the configuration, so they should never be quoted without it.
