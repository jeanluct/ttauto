# Issue #2: implementation plan for the gate test

Status: plan approved 2026-09-19; Step 0 done the same day (see its
"Outcome" below).  Supersedes `ISSUE2_GATES_PLAN.md`.  The mathematics and
the worked example are in `issue2_gates.tex` (built to `issue2_gates.pdf`)
in this directory; the scratch diagnostic that established the result is
`gates_check.cpp`.

## Summary of the problem

`ttauto` accepts a closed folding path as pseudo-Anosov when its main-edge
transition matrix is irreducible with dilatation above 1 (acceptance block in
`include/ttauto/ttauto.hpp`, `descend_graph()`, then `record_pA()` applies
the dilatation window).  That is necessary but not sufficient.
Bestvina-Handel also require that at every vertex the gates be connected
by infinitesimal edges: Prop. 3.3.2 of `Bestvina1995.pdf` (in this
directory) says that a disconnected gate graph gives a reducible map, and
their Section 3.4 (Theorems 3.1.4 and 4.1.4) says that connected gate
graphs give a pseudo-Anosov.
The bad path `{29,46,43,71,88,85,29}` (1-based) in the n=6, stratum-5,
90-vertex subgraph fails that test at exactly one vertex, the punctured
monogon with three edges that the map fixes.  Four control paths reported as
pA by the automaton pass it.

Why the library cannot express the test today:

- `traintrack::fold_infinitesimal_generator` returns the cusp's
  multigon/prong.  The side a fold actually traverses belongs to the
  *target* multigon, between `t_pr` and `t2_pr` in
  `traintrack::fold(multigon&,...)`.  Wrong in all six folds of the bad
  cycle.
- Infinitesimal labels are `mgv` positions (`multigon_prong_index`), which
  are not canonical (structurally equal multigons sort arbitrarily) and are
  not tracked through `normalise()`; every stored one-step map fixes them.
  Edge labels, by contrast, are canonical: they come from the DFS traversal
  (`recursive_get_weights`) started at monogon 0, chosen by
  `minimise_coding`.
- Word order and sign in `fold_traintrack_map` (`include/traintracks/map.hpp`)
  come from fold parity and are not tied to any edge orientation, so
  composed words are not edge paths.  Only the abelianisation was ever
  tested, and `check_fold_map_main_transition` is dead in release builds.
- `tests/test_issue2_bad_path.cpp` models vertices as `(multigon, prong)`
  and treats infinitesimal letters as directions; that marks every monogon
  disconnected regardless of dynamics.
- `ttfoldgraph::add_vertex` identifies a fold result with a stored vertex by
  `operator==` only (coding equality).  At a track with a nontrivial
  automorphism (e.g. vertex 88) the fold result's labelling may differ from
  the stored object's by that automorphism; the graph composes through it
  silently.  Harmless for matrices; a word-level map must use the same
  identification.
- Every `testsuite/**` target is compiled with `-DNDEBUG -O3`, so all
  current `assert`-based CTest tests pass vacuously.  This is why
  `test_map_consistency` never caught the wrong words.
- `master` is ahead of this branch and already contains
  `include/traintracks/coding.hpp` + `lib/traintracks/coding.cpp` (traversal
  internals extracted), normalisation-state guards, a fold-direction
  convention clean-up, and `examples/ttauto_scan_strata.sh/.md` with the
  slow regression `testsuite/ttauto/test_scan_strata_markdown.sh`.

## Design principles

- Namespace `traintracks` for everything about one train track and one
  map: canonical labels, the one-fold map, the gate test.  The test is a
  property of a train-track map, independent of the automaton, and
  `traintracks` is concrete (non-template), so the code lives in `lib/` and
  is unit-testable without a graph.  Namespace `ttauto` only gains the
  path-composition of derivative data and the call in the search.
- New code in new files; existing files get small local edits.  No
  ownership or pointer refactors (see `AGENTS.md`).
- Canonical labels for everything the map mentions, all derived from the
  one DFS that already defines edge labels and cusp numbering, started from
  the same monogon as `weights()`.  Then the graph's identification of a
  fold result with a stored vertex by edge labels is automatically an
  isomorphism on prongs and orientations, and no transport step is needed.
- Never compose words along a path.  Compose derivative data instead
  (Section 6 of the note).
- Every new test uses explicit checks that survive `-DNDEBUG`.

## Step 0: prerequisites (no new logic)

- Merge `master` into `iss002-spurious-pAs`.  Re-run `ctest` and
  `examples/ttauto_min_example`.
- Make the testsuite real.  In `CMakeLists.txt`, `ttauto_add_testsuite`
  adds `target_compile_options(... PRIVATE -UNDEBUG)`.  Add a tiny
  `CHECK(cond)` macro in a new `testsuite/check.hpp` that prints the failing
  expression and file:line and exits nonzero; new tests use it.  Expect some
  existing tests to start failing once their asserts are live; fix or record
  each one.
- Confirm on the merged tree that `gates_check.cpp` still reproduces the
  bad vertex (rebuild, run with no arguments).

Files: `CMakeLists.txt` (3 lines), new `testsuite/check.hpp` (~20 lines),
possibly small fixes in existing testsuite programs.

Outcome (2026-09-19):

- Merge of `master` was conflict-free (commit d78b03c).  Build, `ctest`
  (including the slow strata scan, which is enabled in the local build
  cache) and `examples/ttauto_min_example` all pass; `gates_check.cpp`
  still finds the bad vertex.
- `-UNDEBUG` added to `ttauto_add_testsuite`, plus `testsuite/check.hpp`
  with `CHECK`/`CHECK_MSG`, and `testsuite/` on the include path.
- With asserts live, `test_folding_path_and_badwords.cpp` no longer
  compiled: it called `dim1()`/`dim2()` on a `jlt::matrix`, which has
  `rows()`/`columns()`.  Fixed.  The stale binary had been passing.
- All five fast tests pass with live asserts; `nm` confirms
  `__assert_fail` is referenced in four binaries.  The fifth,
  `test_ttauto_search.cpp`, uses its own local check macro that returns 1,
  so it was never affected by `NDEBUG`.
- `test_map_consistency.cpp` passes with live asserts, including its
  hard-coded words `{1,5,2}` and the "one positive infinitesimal letter"
  check.  This confirms that the current map is self-consistently wrong,
  not that it is right; Step 2 rewrites these expectations.
- Master's coding module is `traintracks::detail::coding_engine`
  (`include/traintracks/coding.hpp`, `lib/traintracks/coding.cpp`) with
  static functions `coding`, `minimise_coding`, `cyclic_symmetry`,
  `print_coding`, `coding_from_monogon`, `recursive_coding`.  The weights
  and cusp traversals (`recursive_get_weights`, `recursive_find_cusp`)
  are still in `traintrack.cpp`.
- Master added `require_normalised(where)` guards to `fold`,
  `fold_cusp_location`, `weights` and `coding`.  Consequence for Step 2:
  the pre-`normalise()` labelling of the folded track must not go through
  those guarded public entry points; implement it as a `coding_engine`
  static (or a `traintrack` private) that walks the raw structure, and
  give it the start monogon explicitly.

## Step 1: canonical prong labels and edge orientation (`traintracks`)

Files: `include/traintracks/coding.hpp` (+~40 lines),
`lib/traintracks/coding.cpp` (+~150 lines),
`include/traintracks/traintrack.hpp` (declare accessor, ~10 lines),
`include/traintracks/map_labels.hpp` (comment: an infinitesimal index is
now a canonical prong label, ~5 lines).

- New struct `traintracks::ttlabels` in `coding.hpp`, built by a new static
  `detail::coding_engine::labels(const traintrack&, int mono)` (and the
  public wrapper `traintrack::labels()` using monogon 0).  For each prong label
  `q`: multigon index, prong index, `k`, `punctured`.  For each edge label:
  tail and head prong labels.  Helpers `side_from(q)`, `side_to(q)`,
  `directions_at(vertex)`.
- New recursion in `coding.cpp` mirroring `recursive_get_weights`: same
  start (monogon 0), same `cycle_edges` order, no descent into uncusped
  monogons.  On first entry into a multigon at prong `pin`, label all its
  prongs `pin, pin+1, ...` (mod `k`) consecutively.  A terminal uncusped
  monogon gets its label when its edge is traversed; monogon 0 gets label
  0.  Edge orientation = direction of first traversal (from the multigon
  that first reaches the edge; monogon 0's edge points outward).
- Side label `q` = side from prong `q` to the next prong in the
  `cycle_edges` direction (`p <- mod(p+1,k)`), the same direction that
  `t2_pr = mod(t_pr+dir,k)` uses in `fold`.
- The recursion takes its start monogon from the same place as `weights()`
  and asserts (fail-fast) that it equals what `fold_cusp_location` finds.
- Accessor `traintrack::labels() const`.  Also a variant that labels from an
  arbitrary start monogon, needed by Step 2.

## Step 2: correct one-fold map (`traintracks`)

Files: new `include/traintracks/fold_map.hpp` (~80 lines), new
`lib/traintracks/fold_map.cpp` (~200 lines), `include/traintracks/map.hpp`
(replace the body of `fold_traintrack_map`, ~40 lines changed),
`include/traintracks/traintrack.hpp` (declare one member; retire
`fold_infinitesimal_index/generator` or redefine them as "signed side of
the target"), `lib/traintracks/traintrack.cpp` (fix the docstrings; the
header's "global infinitesimal-loop convention" is false).

- `struct fold_map_data`, all in canonical labels: `moved`, `onto` (old
  edge labels), `edge_image[e]` (signed new label of each old edge), `side`
  (signed new side label traversed), `moved_first` (word order),
  `prong_image[q]` (new label of each old prong).
- `traintrack::fold_with_map(int f, fold_map_data&)`: compute pre-fold
  labels; locate `e0`, `e1`, target, `t_pr`, `t2_pr` via
  `fold_cusp_location` and `target_multigon`; perform the structural fold.
  Preferred correspondence: compute post-fold labels *before*
  `normalise()`, between `insert_edge` and `normalise()` in
  `fold(multigon&,...)`, starting from the minimising monogon that
  `minimise_coding` already computes.  Edge and multigon objects are then
  still in place and the pre/post correspondence is by pointer.  Fallback:
  after `normalise()`, match prongs by slot contents (edge-object identity
  survives `swap(multigon&,multigon&)`, ending indices do not), as
  `gates_check.cpp` does.
- `traintracks::fold_traintrack_map(tt0, f)` builds the
  `jlt::freeauto<int>` from `fold_map_data`: unmoved edges map to their
  signed new label; sides map through `prong_image`; the moved edge maps to
  `e0'.S.(+-e1')` or `(+-e1').S.e0'`.  Template signature unchanged, so
  `ttfoldgraph` needs no edit.  Do not pursue the "one fold with distinct
  weights" shortcut: weights collide because `fold` sets `w1 <- w0 + w1`.
- Rewrite the hard-coded word expectations in
  `testsuite/traintracks/test_map_consistency.cpp` and
  `tests/test_ttmap.cpp` in the same commit.
- Regression kept: the abelianisation equals `fold_transition_matrix`.

## Step 3: gates module (`traintracks`)

Files: new `include/traintracks/gates.hpp` (~120 lines), new
`lib/traintracks/gates.cpp` (~300 lines).

- `struct fold_derivative`: `D` on all signed letters (main and side) and
  the turns `T` of the one-step word, recorded at the *target* track's BH
  vertices.  A permutation fold has no turns; a fold onto an unpunctured
  target has one turn between two gates; a fold onto a punctured target has
  two turns at two prong-vertices.  Built from a one-step
  `jlt::freeauto<int>` (at most three letters per image) plus the target's
  `ttlabels`.
- `class gate_accumulator`: `push_back(const fold_derivative&)` maintains
  the composed `D` and the turn set `T <- T_i union D_i(T)`;
  `analyse(const ttlabels& start)` closes `T` under `D x D`, builds gates
  per BH vertex by `D^k` coincidence (iterate on the finite direction set to
  a fixed point), joins, and union-find connectivity.  Returns
  `gate_analysis` with `connected`, per-vertex records and `print()`.
- Word-based overload `analyse_gates(const traintrack&, const
  jlt::freeauto<int>&)` for tests and diagnostics only (what
  `gates_check.cpp` does); it must agree with the accumulator.
- Fail-fast (`std::cerr` + `exit(1)`, house style): realised turn inside a
  gate; turn whose ends are at different vertices; `D` not mapping a
  vertex's directions to a single vertex; side permutation not preserving
  `k` or puncturedness.
- Consequences to use as assertions (note, Corollaries 3.4 and 3.5): at a
  prong of a punctured multigon the only possible infinitesimal edges are
  incoming-peripheral to first main gate and last main gate to outgoing
  peripheral, so connectivity there means exactly one main gate at the
  prong and both peripheral turns realised; at an unpunctured k-gon at
  most one prong is refined, into two gates.  The general union-find stays
  the implementation; these are cross-checks on its output.
- Shape check (Bestvina-Handel Props. 3.3.3-3.3.4): when the gate graph is
  connected, every vertex must have at least two gates and its
  infinitesimal edges must form a single edge (two gates), a full polygon
  on adjacent gates, or a polygon with exactly one side missing.  Report
  any other shape as a diagnostic warning; on a connected candidate it
  indicates a bug in the map or the label bookkeeping, not a property of
  the braid.

## Step 4: integrate into the search (`ttauto`)

Files: `include/ttauto/folding_path.hpp` (~25 lines),
`include/ttauto/ttauto.hpp` (~40 lines).  `ttfoldgraph.hpp` untouched:
`fold_derivative` is derived from the stored `AMv` words at check time, cost
O(L (n + ninf)) per candidate.

- `folding_path::gates() const`: walk `fp`/`vp`, feed a
  `traintracks::gate_accumulator` from `ttg->traintrack_map(v,f)` and
  `ttg->traintrack(target).labels()`, analyse at the initial vertex.
- `ttauto::record_pA()`: after the `lambdamax`/`lambdamin` window checks
  and before class insertion, `if (do_check_gates && !p.gates().connected)
  { ++gate_rejected; return; }`.  Setter `check_gates(bool)`; the counter is
  printed with the other totals in `search()`.  `pAclass` unchanged.
  Badwords and `check_all_norms` are unaffected (they act before the
  closed-path test).
- Default: on (decision of 2026-09-19).  Rejections are counted and
  printed; `check_gates(false)` restores the old behaviour.

## Step 5: tests

- `testsuite/traintracks/test_labels.cpp`: label tables for n=3..5
  fixtures; counts equal `total_prongs()` and `edges()`; tails and heads are
  valid prong labels; labels of a copy equal the original's; labels are
  invariant under `normalise()` of an already-normalised track.
- `testsuite/traintracks/test_fold_map_paths.cpp`: for every vertex and
  fold of `n=3,trk=0` and `n=4,trk=1`, the one-step word is a continuous
  path in canonical labels and abelianises to `fold_transition_matrix`.
  Random paths (generator pattern from `tests/test_ttmap_from_path.cpp`,
  including vertices with nontrivial `cyclic_symmetry()`) compose to words
  that are continuous for three iterates and abelianise to
  `transition_matrix()`.
- `testsuite/traintracks/test_gates.cpp` (label `slow`, or keep it fast by
  building only the needed subgraph once): n=6, trk=4, first pruned
  subgraph.  The bad cycle is rejected at exactly the fixed 3-edge monogon;
  the four controls `{28,45,42,45,28}`, `{2,37,39,19,22,20,12,2}`,
  `{4,84,28,45,28,80,4}`, `{21,42,45,28,45,42,21}` (0-based; branch
  sequences selected by matrix as in `gates_check.cpp`) pass; the
  word-based and accumulator analyses agree.  Assert the puncture
  corollary on every vertex of every control (one main gate per punctured
  prong, both peripheral turns realised) and that the bad vertex fails it
  by a refined prong.  Print gate partitions so that at least one control
  exhibits a refined prong at an unpunctured multigon that is still
  connected; if none does, say so in the test output.
- `testsuite/ttauto/test_ttauto_gates.cpp`: `ttauto::search` with
  `max_pathlength(7)` on that subgraph.  The polynomial
  `x^7 - 2x^6 + x^5 - 4x^4 + 4x^3 - x^2 + 2x - 1` is present with gates off
  and absent with gates on.
- Over-rejection guard: run the gate test over every class reported by
  `examples/ttauto_min_example` (n=3..5 literature minima); none may be
  rejected.  Re-run `examples/ttauto_scan_strata.sh` and diff against
  `examples/ttauto_scan_strata.md`.  Any class that disappears must be
  explained (all its representatives fail the gate test), and the paper's
  tables (`pubs/ttauto paper`, commit bd47aa7) revisited.
- Resolve the fifth control `{1,18,65,18,56,57,57,1}`, where
  `gates_check.cpp` failed `tt == ttg.traintrack(next)` at 57 -> 58.  The
  scratch program matches branches by `target_vertex` while `add_vertex`
  enumerates non-identity folds; find the divergence before trusting the
  controls.
- Retire `tests/test_issue2_bad_path.cpp` and `devel/iss002/gates_check.cpp`
  once `test_gates.cpp` covers them, or keep `gates_check.cpp` as a
  documented diagnostic.  Update `.gitignore` for new `tests/` binaries.

## Step 6: documentation and close-out

- `CODE_STRUCTURE.md`: subsections for labels, fold_map, gates; the "one
  fold through the stack" walkthrough now records the target side and
  `record_pA` is gated.
- `testsuite/COVERAGE.md`, `TESTING.md`, `AGENTS.md` repository map,
  `CLAUDE.md` issue-2 section, `devel/macro-status.md` if a flag is added.
- `ISSUE2_BAD_PA_SUMMARY.md`: "resolved by" section.
- Close GitHub issue #2 pointing at the note and the tests.

## Order and commits

One commit per step, build green and `ctest` passing after each.  Steps 1
and 2 land together with the updated word expectations, since they change
the infinitesimal letters of every stored map.  Step 4 is the only
user-visible behaviour change and goes last, with its test.

## Risks

- Gates finer than prongs at unpunctured multigons make the BH condition
  stricter than "all sides realised".  The controls do not exercise this,
  so the `ttauto_min_example` and strata-scan sweeps are the real guard
  against over-rejection.
- Automorphic vertices: the canonical-label argument removes the transport
  step, but `test_fold_map_paths` must include such vertices to prove it.
- Published numbers may change (Step 5 diff).  That is the point of the
  exercise, but it needs a note in the paper.
