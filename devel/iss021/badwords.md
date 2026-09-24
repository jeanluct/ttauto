# Bad-word pruning: what it costs and what it buys

Measurements behind issue #21, taken 2026-09-24 on this machine with the
CMake Release build (`-O3 -ffast-math`).  Class counts and omitted counts
are deterministic; the timings are not, and are only meant to show orders
of magnitude.

## What changed

Before, the prune was gated on `max_badword_length > 0 && do_check_norms`,
so it could only fire in the norm-bounded mode, and the constructor
defaulted `max_badword_length` to 2.  Every caller in the repository
turned it off explicitly, except `examples/ttauto_labels.cpp`, which
never set it and so inherited the default.

Now `badword_length()` alone controls it, in either mode, and it defaults
to 0.

## It never moves the minimum

Across every configuration below, the lowest dilatation found is
identical with and without pruning.  That is the claim the paper makes
(`pubs/ttauto paper/ttauto.tex:1645-1651`: a repeated loop "can only
raise the dilatation, so the minimizer cannot involve such repeated
loops") and nothing had ever executed it.

## It does drop non-minimal classes

This is the part that matters, and it is why the prune is off by default.

Length-bounded search, no dilatation window, four punctures, stratum 2
(3 vertices):

| length | classes off | classes on | omitted | time off | time on |
|---:|---:|---:|---:|---:|---:|
| 8  | 21  | 15  | 29  | 0.005 s | 0.003 s |
| 10 | 62  | 39  | 93  | 0.028 s | 0.011 s |
| 12 | 181 | 98  | 293 | 0.141 s | 0.041 s |
| 14 | 526 | 267 | 917 | 0.637 s | 0.148 s |

Roughly half the classes disappear at the longer lengths.  Checked by
characteristic polynomial, the pruned set is always a proper subset of
the unpruned one -- 0 classes appear only with pruning.  The length-10
row is pinned in `testsuite/ttauto/test_badword_pruning.cpp`.

Note this must be compared by polynomial, not by dilatation: the same
class reached along a different path differs in the last few digits, so
a tight floating-point comparison spuriously reports a mismatch.

Five punctures behaves the same way (stratum 3, 3 vertices: 625 -> 316
classes at length 14; stratum 4, 9 vertices: 446 -> 334).

## A dilatation window hides the loss

Six punctures, stratum 2 (`traintrack(6,3)`, 138 vertices), norm-bounded:

| window | length | classes off | classes on | omitted |
|---:|---:|---:|---:|---:|
| 2.5 | 4 | 2  | 2  | 61    |
| 2.5 | 6 | 5  | 5  | 1577  |
| 5   | 6 | 11 | 11 | 1586  |
| 20  | 6 | 11 | 11 | 1586  |
| 2.5 | 8 | 5  | 5  | 49824 |

No class is lost at any of these.  The reason is that a path traversing
a loop twice has a higher dilatation than the collapsed one, so a window
tight enough to be interesting rejects it anyway.  The loss only shows
once the window is loose relative to the path length -- at length 8 with
no window at all, the same stratum gives 43 classes against 40.

So the old coupling to `check_norms` was hiding the problem rather than
preventing it: norm-bounded mode is *usually* safe in practice, but only
because the window does the work, and a user collecting everything below
a generous bound would still have lost classes silently.  That is what
the warning at the end of `ttauto::search` is for.

## Is it worth having

Yes, but the benefit varies by more than an order of magnitude, and it is
not where the earlier guess ("slow, or not worth it") put it.

- Small automaton, long paths: **4.3x faster** (n=4 stratum 2, length 14,
  0.637 s -> 0.148 s).  Most of the search is loops round a few vertices,
  and most of those loops repeat.
- Large automaton, short paths: **1.1x** (six punctures, 138 vertices,
  length 8, 8.68 s -> 7.71 s).  Paths are short relative to the graph, so
  few of them repeat anything.

The `badwords()` precomputation is negligible at these sizes (2 to 20
words, built in milliseconds); it is the table *search* per step that
costs, which is why `badword_length(1)` is as good as `badword_length(3)`
in every run above -- all the bad words found on these strata have
half-length 1.

Recommendation: leave it off by default, use it when hunting a minimum on
a small stratum with a long path bound, and never use it to enumerate.
