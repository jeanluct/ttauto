# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`ttauto` is a C++17 library (plus example programs) that builds train-track
folding automata for homeomorphisms of punctured discs and searches them for
pseudo-Anosov (pA) candidates.  Authors: Jean-Luc Thiffeault and Erwan
Lanneau.  GPLv3.

Companion documents already in the repo; read them rather than rediscovering:

- `AGENTS.md`: style, ownership model, testing expectations, refactor rules.
- `CODE_STRUCTURE.md`: class-by-class tour, glossary, worked "one fold
  through the stack" example, suggested reading orders.
- `TESTING.md` and `testsuite/COVERAGE.md`: test matrix and header-to-test
  map.
- `devel/macro-status.md`: every non-guard feature macro and how to build
  with it.
- `extern/jlt/AGENTS.md`: overrides `AGENTS.md` for the bundled `jlt`
  submodule (header-only numerics; Catch2 tests under `extern/jlt/tests`).

## Build and test

CMake is primary; SCons files are legacy and still work (`scons`).

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
`examples/ttauto_count` takes about 2 minutes and `examples/ttauto_labels`
about 100; run them deliberately.

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
  objects keyed by characteristic polynomial.

Typical pipeline: `build_traintrack_list` -> pick a `traintrack` ->
`ttfoldgraph<traintrack>` -> optional `subgraphs(...)` -> `ttauto::search`
-> `pA_list()`.  `examples/ttauto_min_example.cpp` is the shortest complete
instance.

Conventions worth knowing:

- Codings are printed as 5-digit units; the middle digit is the multiprong
  label controlled by `traintrack::label_multiprongs`.  The paper only uses
  4 digits.  Only `print_coding` honours the flag; string input does not.
- User-facing output (the `ttauto` program, Mathematica files, the devel
  notes) is 1-based; C++ vertex, branch and fold indices are 0-based.
- Per-class `static constexpr int debug = 0;` members gate verbose tracing.
- `jlt` types (`jlt::mathmatrix`, `jlt::vector`, `jlt::freeword`,
  `jlt::freeauto`) are used throughout; do not reimplement them.
- Mathematica post-processing lives in `mathematica/` (`TrainTracks.m`,
  saved `ttauto_output/`).  Notebooks go through the `dropoutput_nb` filter.

## Current work: issue #2 (branch `iss002-spurious-pAs`)

Issue #2: the search very rarely reports a pA that is actually reducible.
Matrix primitivity and a plausible dilatation are not sufficient; the
Bestvina-Handel gate condition must also hold.  Root cause established
2026-09-19: read `devel/iss002/issue2_gates.tex` (or its PDF) first, then
`devel/iss002/ISSUE2_IMPLEMENTATION_PLAN.md` for the approved fix plan.
The scratch diagnostic `devel/iss002/gates_check.cpp` reproduces the
result.  Older material in `devel/iss002/`:

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
growth `2.01536` yet reducible.  `tests/test_issue2_bad_path.cpp` is the
hardwired reproducer and current sandbox for the gate pipeline; it builds
with the normal CMake build and runs in seconds.  Its gate-connectivity
result is provisional: the `(multigon,prong)` vertex model does not yet
match the switch-level vertex model of the `.train` output, so do not use
its connectivity verdict as a classification criterion yet.  Treat `.train`
edge numbering as display-only and keep canonical internal labels from
`ttmap_labeler` for all logic.

Related: `jlt_todo.md` notes that the reproducer may be misusing
"infinitesimal" edges where "peripheral" edges are meant.  Issue #3 (train
track map) is closed and is the prerequisite this work builds on.
`TODO.md` holds the author's loose task list.
