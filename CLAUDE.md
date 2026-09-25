# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`ttauto` is a C++17 library (plus example programs) that builds train-track
folding automata for homeomorphisms of punctured discs and searches them for
pseudo-Anosov (pA) candidates.  Authors: Jean-Luc Thiffeault and Erwan
Lanneau.  GPLv3.

Companion documents already in the repo; read them rather than rediscovering:

- `doc/ttauto.tex` (build with `latexmk -pdf` in `doc/`): the mathematics
  and implementation of the pseudo-Anosov test, from train tracks and folds
  to the Bestvina-Handel gate condition, with a worked n=4 example and a
  concept-to-file appendix.  Read this before touching the acceptance
  logic.
- `AGENTS.md`: style, ownership model, testing expectations, refactor rules.
- `doc/CODE_STRUCTURE.md`: class-by-class tour, glossary, worked "one fold
  through the stack" example, suggested reading orders.
- `doc/TESTING.md` and `testsuite/COVERAGE.md`: test matrix and header-to-test
  map.
- `devel/macro-status.md`: every non-guard feature macro and how to build
  with it.
- `extern/jlt/AGENTS.md`: overrides `AGENTS.md` for the bundled `jlt`
  submodule (header-only numerics; Catch2 tests under `extern/jlt/tests`).

## Build and test

CMake is the build system; the legacy SCons files were removed in 2026.

```bash
cmake -S . -B build
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Outputs are written in place, not into `build/`: library to `lib/`, example
binaries to `examples/`, test binaries to `tests/`.  Only `testsuite/**/*.cpp`
programs are registered with CTest, and those land in `build/testsuite/`.
Sources in `examples/`, `tests/`, and `testsuite/` are globbed, so a new
`.cpp` needs no CMake edit.

Target names are prefixed by folder:

```bash
cmake --build build --target example_ttauto_min_example   # examples/*.cpp
cmake --build build --target test_test_traintrack         # tests/*.cpp
cmake --build build --target testsuite_test_map_consistency  # testsuite/
ctest --test-dir build -R test_ttauto_search              # one CTest program
```

Tests are plain `main()` programs with no framework: they print and exit
nonzero on failure.  Fast manual matrix after touching `include/` or `lib/`:
`tests/test_traintrack`, `tests/test_folding_path`, `tests/test_permplus1`,
`tests/test_badwords`, `examples/ttauto_min_example`, `examples/ttauto_torus`.
`examples/ttauto_count` takes about 2 minutes; run it deliberately.
`examples/ttauto_labels` took about 100 minutes until it capped its path
length: `check_norms()` derives that bound from the dilatation window, and
searching to it is what makes the result complete for the window, so the
cost was real rather than waste.  See `devel/iss016/labelled_automaton.md`.

`examples/ttauto` is interactive; feed it a defaults-only stdin
(`examples/ttauto <<EOF` / `EOF`) or an answer file such as
`devel/iss002/ttauto.input`.  It and the search testsuite write Mathematica
files (`ttauto_n=*.m`, `pA_n=*.m`) into the current directory; delete those
before committing.

Compile flags are `-Wall -O3 -ffast-math`.  To build a non-default macro
variant, configure a separate build directory with
`-DCMAKE_CXX_FLAGS='-D<MACRO>'` (see `devel/macro-status.md`).

## Architecture

Two namespaces, two layers:

- `traintracks` (`include/traintracks/`, implemented in `lib/traintracks/`):
  the concrete combinatorial model.  `traintrack` owns `multigon`s (via
  `unique_ptr`), each `multigon` holds edge slots per prong, each `edge`
  holds non-owning back-links to its two endpoints.  The invariant is that
  edge endpoint records and multigon slots always agree; `check()` methods
  verify it and code fails fast (`std::cerr` + `std::exit(1)`) on
  violation.  `normalise()` puts a track in canonical form; fold indices and
  `coding()` are only meaningful after it.  `build.hpp` enumerates the
  starting tracks for a stratum.  `map.hpp` turns one fold into both a
  `mathmatrix_permplus1` transition matrix and a `jlt::freeauto<int>`
  train-track map; `map_labels.hpp` fixes the generator numbering (main
  edges `1..nmain`, infinitesimal edges after, sign = orientation).
- `ttauto` (`include/ttauto/`, header-only templates on `TrTr`):
  `ttfoldgraph` is the automaton, built by recursively applying every fold
  and deduplicating via codings; each branch stores its target vertex, its
  one-step matrix and its one-step map.  `folding_path` is a sequence of
  branch choices that composes those into a path matrix and path map.
  `ttauto` runs a DFS over folding paths with pruning (norm bounds,
  `badwords`, length caps) and groups accepted closed paths into `pAclass`
  objects keyed by characteristic polynomial.  `path_braid.hpp` reads the
  braid of a closed path.

Typical pipeline: `build_traintrack_list` -> pick a `traintrack` ->
`ttfoldgraph<traintrack>` -> optional `subgraphs(...)` -> `ttauto::search`
-> `pA_list()` -> `ttauto::folding_path_braid` for the braid of an
accepted path.  `examples/ttauto_min_example.cpp` is the shortest complete
instance.

Conventions worth knowing:

- A printed coding is one whitespace-separated block per walk step.  A
  block is four fields, `(prong, nprongs, edge, nedges)`, as in the paper;
  it gains the label as a third field when some multigon carries one, or
  when `print_coding(strm, dir, true)` forces it.  Fields run together as
  digits while each is a single digit, and are hyphen-separated otherwise
  (`1-1-12-1`), so the width is always carried per block and never
  inferred.  `traintracks::parse_coding` and `traintrack(const char*)`
  read both forms back; a bare list of integers is rejected.
- User-facing output (the `ttauto` program, Mathematica files, the devel
  notes) is 1-based; C++ vertex, branch and fold indices are 0-based.
- Per-class `static constexpr int debug = 0;` members gate verbose tracing.
- `jlt` types (`jlt::mathmatrix`, `jlt::vector`, `jlt::freeword`,
  `jlt::freeauto`) are used throughout; do not reimplement them.
- Mathematica post-processing lives in `mathematica/` (`TrainTracks.m`,
  saved `ttauto_output/`).  Notebooks go through the `dropoutput_nb` filter.

## Current work: issue #4 (branch `iss004-braid-from-path`)

Issue #4: read the braid off a closed folding path.  Done.  A closed path
defines a homeomorphism of the punctured disc; turning it into a word in
the braid generators needs the punctures to have positions, which is what
`traintracks::outer_embedding` (`embedding.hpp`) supplies: the multigons
and main edges form a tree, so the coding's cyclic orders are a complete
rotation system and the planar embedding is fixed up to reflection.  A
fold onto an unpunctured multigon moves nothing; one onto a punctured
monogon takes the moved edge round the puncture, and undoing that swaps
two adjacent blocks of punctures (`traintracks::fold_block_swap`).
`ttauto::folding_path_braid` (`path_braid.hpp`) accumulates those and
reports a `verified` flag: `braidword::growth`, the growth under the
Dynnikov action, must equal the Perron root of the path's own matrix.
Two traps, both of which caused wrong braids before being found: a branch
index is not a fold index, and each step must continue from the
automaton's own copy of the track, since at a cyclically symmetric track
that need not be the physical track you just folded.  See
`doc/ttauto.tex` sections "The braid of a closed path" and "Reading the
braid off a path", `examples/ttbraid`, and
`examples/ttauto_strata_braids.md`, which checks the braid of every
stratum minimiser for n=3..7 against the published table.

## Earlier work: issue #2 (branch `iss002-spurious-pAs`)

Issue #2: the search reported some pAs that are actually reducible.  An
irreducible matrix and a plausible dilatation are not sufficient; the
Bestvina-Handel gate condition must also hold.  Fixed 2026-09-19: the
search now requires a primitive matrix (not just irreducible) and runs
the gate test in `ttauto::record_pA` (default on,
`check_gates(false)` to disable; rejected classes in
`rejected_pA_list()`).  The machinery is `ttnumbering` (coding module),
`fold_map` and `gates` in `traintracks`.  Read
`devel/iss002/issue2_gates.tex` (or its PDF) for the mathematics and
`devel/iss002/ISSUE2_IMPLEMENTATION_PLAN.md` for what was done and what
remains (follow-ups).  `tests/ttauto_gate_census` lists every rejected
class across strata: for n=3..7 the gate test rejects five classes, all
"lower-puncture pA plus an idle puncture", and the primitivity
requirement removes two more with imprimitive matrices.  Test coverage
(86% of library lines, measured 2026-09-19) and how to re-measure it are
in `testsuite/COVERAGE.md` and `doc/TESTING.md`.  Older material in
`devel/iss002/`:

- `ISSUE2_BAD_PA_SUMMARY.md`: the evidence, current status, work log.
- `ISSUE2_GATES_PLAN.md`: phased plan (directed alphabet per vertex,
  derivative map `D`, gates as classes under `D^k(a) = D^k(b)`, gate graph
  connected via infinitesimal edges, regression, then promote helpers into
  the library).
- `analyze_n=6_5_bad_tt_train.md`: hand analysis of the trusted external
  `.train` output; the gate-connectivity failure is at its vertices 2 and 6.
- `n=6_5_bad_tt.train`, `n=6_5_bad_tt_data.m`, `n=6_5_bad_tt.nb`,
  `toby_hall-email_2009-04-30.pdf`, `notes_badbraid.pdf`: source data.
- `ttauto.input` / `ttauto.output`: answers that drive `examples/ttauto` to
  the bad case (6 punctures, stratum 5, subgraph 1 of 8, the 90-vertex one).

The canonical bad case is `n=6`, `trk=4`, `sgidx=0`, cycle
`{29,46,43,71,88,85,29}` (1-based; `{28,45,42,70,87,84,28}` in C++) with
branch sequence `{1,0,3,2,1,2}`.  Braid `1 2 3 -5 -4 -3` on 6 strings;
growth `2.01536` yet reducible.  That braid was worked out by hand in
2009; since issue #4 the code reads its own braid off the same path and
gets `1 2 1 2 3 4 5 5 -4 -3 -5 -5 -4 -3 -2 -1`, which Trains also calls
reducible at `2.01536` and which fixes one puncture, so the two agree in
everything conjugation preserves and differ only in which puncture is
idle.  It is pinned in `testsuite/ttauto/test_braid_extraction.cpp`.
Treat `.train` edge numbering as display-only and keep canonical internal
labels from `ttmap_labeler` for all logic.

Issue #3 (train track map) is closed and was the prerequisite; its word
conventions were corrected as part of issue #2 (the side letter is now the
target multigon's side, sides are permuted, orientations are canonical).
Testsuite programs are built with `-UNDEBUG`; use `CHECK()` from
`testsuite/check.hpp` in new tests.
`TODO.md` holds the author's loose task list.
