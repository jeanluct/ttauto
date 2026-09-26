# Issue #13: drawing the collapsed train track

State of the work on branch `iss013-embedded-plot`, paused 2026-09-26.
Read this before resuming.  `devel/plans/ttplot_tikz_plan.md` is the
plan for the older `iss013-plot-traintrack` branch and is superseded.

## Goal

Figures that closely approximate the hand-drawn (Inkscape) figures of the
ttauto paper, `ttauto.tex` in the parent repository, as TikZ.  The
convention is that of the paper's Figures 3, 4 and 7, not Figure 12
(the one radial figure) and not `doc/ttauto.tex`: punctures on an
invisible horizontal axis, the track properly embedded in the upper half
plane.  For now the multigons are collapsed to points and nothing is
labelled.

Two rules from the author that the construction follows:

- Arcs never cross; where two meet they are tangent.
- A multigon's prongs are equally spaced by angle, edges attach to prongs,
  and edges at the same prong are tangent there.

## Where things are

- `lib/traintracks/collapsed_layout.cpp`,
  `include/traintracks/collapsed_layout.hpp`: `make_collapsed_layout`,
  the whole drawing.  Each edge comes out as a chain of cubic Beziers.
- `examples/ttplot.cpp`: turns the layout into TikZ (`--snippet`,
  `--coding`, `--output`, `--scale`, `--labels`).
- `testsuite/traintracks/test_collapsed_layout.cpp`: asserts no proper
  crossing, no curve running past its own ends, and nothing below the
  axis, at all 3700 automaton vertices for n=3..7 (about 6 s).  A
  self-check proves it sees two arcs crossing after leaving the same
  point; an earlier checker excused all such pairs and so reported no
  crossings where there were hundreds.
- `latex/test_ttplot`: renders the Figure 12 track and the six tracks of
  the issue-2 bad path into `test_ttplot_all.pdf`.
- `find_crossings.cpp`, `gallery.py` here: list the tracks for a given n
  whose drawing crosses, and render them into a PDF, crossings circled.
- `n8_failures.txt` here: that list for n=8 at `6966969`, the baseline
  to compare against.

The two tracks of the paper's Figure 4, identified by the author, are
vertices 2 and 0 of `build_traintrack_list(4)[1]`:

    (a)  1111 1122 1112 1311 2311 1111 3311 1111
    (b)  1111 1311 2311 1111 3312 1111 3322 1111

In the paper, 4(b) has puncture 2 moved to the fourth position.  That is
a different cut of the boundary walk from the canonical one, most likely
the one obtained by carrying 4(a)'s punctures through the fold.  See
"Not started" below.

## How the drawing is built

1. `outer_embedding` gives the boundary walk, cut at puncture 1's loop,
   and so the punctures' order along the axis.  The tree of multigons is
   rooted at puncture 1.
2. Each edge cuts off a block of consecutive punctures.  The blocks are
   nested or disjoint, and each gets a box over its stretch of axis,
   narrowed a little at each depth so that nested boxes sit strictly
   inside one another and siblings do not touch.
3. An edge crosses exactly one box boundary, the top of its child's box,
   straight above the child and vertically.  It is two cubics joined
   smoothly there, one inside the child's box and one in the parent's, so
   curves in different boxes cannot meet: an arc diagram.
4. Built bottom up: a multigon sits half a unit (0.25) above its
   children's boxes, its box is closed above everything drawn in it, and
   its children's pieces are built nearest first.
5. Each prong star's rotation is the circular mean of where its edges go
   next: straight up for the edge to the parent, and the child box's
   entry point for the others.  Prong index runs anticlockwise, the sense
   of the boundary walk.
6. A piece's two control lengths are searched over scales 0.2 to 5 of a
   first guess, keeping the lowest piece that stays out of every
   sibling's box, above pieces already built, and above the axis.
7. Every control point is capped to the x-range of its piece's two ends,
   plus 0.1, so no curve runs past the point it is heading for.  This was
   the author's diagnosis of the last failures, and it took n=7 from 14
   crossing tracks to none.
8. Edges tangent at a shared prong keep their left-to-right order through
   their curvature there, which must increase along the slots, with a
   margin.  A final pass moves pairs apart, with each control length
   kept within 0.2 to 3 times its start.

## Measurements at `6966969`

| n | tracks | with crossings | tallest |
|---|---|---|---|
| 3..6 | 428 | 0 | 4.5 |
| 7 | 3272 | 0 | 6.6 |
| 8 | 31370 | 184, in 7 of 19 strata | 9.4 |

Along the way, n=3..6: 359 (true count once the checker was fixed) ->
185 (prong sense) -> 133 (gap rotation) -> 82 (box routing) -> 80
(rotation fitted to entry points) -> 65 (clearance measured, not
estimated) -> 46 (curvature order) -> 27 (with a margin, after building
bottom up) -> 0 (clearance checked on every sample, against every
sibling's box).

## Lessons that cost time

- **Check the checker.**  A test that excuses arcs sharing an endpoint
  tests almost nothing in a collapsed drawing.  Mutation-check anything
  that is meant to catch crossings.
- **Prong sense.**  Numbering prongs clockwise mirrors the prongs but not
  the slots within a prong, which forces a crossing.  It must match the
  walk.
- **Slot order at a punctured monogon** runs the same way as elsewhere;
  reversing it there made things far worse (133 -> 206).
- **Don't estimate a curve's height: sample it.**  Aiming the middle of a
  cubic at a peak let its ends dip into boxes.
- **Every lengthening loop needs a bound.**  Unbounded ones produced
  drawings 31000 and 1349 units tall.  The second came from the curvature
  pass: lengthening only takes a curvature towards zero, never past it.
- **Ties in "nearest first".**  Siblings at equal distance meant a piece
  was never checked against a box built after it; boxes are now checked
  against every sibling.
- **Don't read topology off a picture.**  Twice a crossing was blamed on
  the wrong pair of edges; locate it numerically.

## Next steps

1. **Cosmetics**, where the author wanted to go next.  The cap in step 7
   made the drawings angular: short controls leave corners and nearly
   straight runs (see both Figure-4 panels).  There are vertical stubs up
   to each box top, and flat-topped arcs.  Change one thing at a time,
   rerun `test_collapsed_layout`, and regenerate `latex/fig4.pdf` and
   the galleries.
2. **n=8**: 184 tracks still cross.  The first ones in the gallery look
   like the old n=7 failures: a long arc from a multigon near the left
   clips a neighbour as it comes down.  Cosmetic work may fix some.
3. **Promote the test** to n=8 once planar there; at n=8 it would take
   about two minutes, so probably as a slow test
   (`TTAUTO_ENABLE_SLOW_TESTS`).

Not started: labels; the full (uncollapsed) rendering of the paper's
Figures 3, 4 and 7, with teardrops and polygons; drawing a sequence of
tracks along a folding path with the punctures carried through each fold
(`transported_cut_dart`, `fold_block_swap`), as 4(b) seems to do.

## Regenerating the galleries

From the repository root, with the CMake build done:

    g++ -std=c++17 -O2 -Itestsuite -Iinclude -Iextern/jlt \
      -Iextern/jlt/extern/CSparse/Include devel/iss013/find_crossings.cpp \
      lib/libttauto.a build/libcsparse.a -lm -o /tmp/find_crossings
    /tmp/find_crossings 8 >| latex/n8_failures.txt     # about 2 minutes
    cd latex
    python3 ../devel/iss013/gallery.py 8 n8_failures.txt n8_failures.tex "31\,370"
    pdflatex n8_failures.tex

The Figure-4 pair: `examples/ttplot --snippet --coding "<coding>"
--output fig4a.tex` for each of the two codings above, then a small
document that `\input`s them inside `\scalebox{2.6}{...}`.

Galleries and examples for the author go in `latex/`, untracked;
overwrite old ones and clean up now and then.
