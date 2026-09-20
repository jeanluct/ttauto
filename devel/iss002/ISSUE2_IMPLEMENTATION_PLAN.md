# Issue #2: implementation plan for the gate test

Status: plan approved 2026-09-19; Steps 0 to 4 and the census of Step 5
done the same day (see the "Outcome" notes below).  Supersedes
`ISSUE2_GATES_PLAN.md`.  The mathematics and the worked example are in
`issue2_gates.tex` (built to `issue2_gates.pdf`) in this directory; the
scratch diagnostic that established the result is `gates_check.cpp`.

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

- Terminology: the new canonical index of prongs and edges is a
  *numbering* (`ttnumbering`, `traintrack::numbering()`, "prong number",
  "edge number").  "Label" is already taken in this code for the puncture
  label of a multigon (`multigon::label()`, `set_label`, `pure_braid`, the
  5th coding digit) and must not be reused for it.
- Namespace `traintracks` for everything about one train track and one
  map: canonical numbering, the one-fold map, the gate test.  The test is a
  property of a train-track map, independent of the automaton, and
  `traintracks` is concrete (non-template), so the code lives in `lib/` and
  is unit-testable without a graph.  Namespace `ttauto` only gains the
  path-composition of derivative data and the call in the search.
- New code in new files; existing files get small local edits.  No
  ownership or pointer refactors (see `AGENTS.md`).
- Canonical numbers for everything the map mentions, all derived from the
  one DFS that already defines edge order and cusp order, started from the
  same monogon as `weights()`.  Then the graph's identification of a fold
  result with a stored vertex by edge number is automatically an
  isomorphism on prongs and orientations, and no transport step is needed.
  This adds a fourth copy of that DFS (coding, weights, cusp location,
  numbering); consolidating them is a flagged follow-up (see the end of
  this file), not part of this work.
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

## Step 1: canonical prong numbering and edge orientation (`traintracks`)

Files: `include/traintracks/coding.hpp` (+~40 lines),
`lib/traintracks/coding.cpp` (+~150 lines),
`include/traintracks/traintrack.hpp` (declare accessor, ~10 lines),
`include/traintracks/map_labels.hpp` (comment: an infinitesimal index is
now a canonical prong number, ~5 lines).

- New struct `traintracks::ttnumbering` in `coding.hpp`, built by a new
  static `detail::coding_engine::numbering(const traintrack&, int mono)`
  (and the public wrapper `traintrack::numbering()` using monogon 0).  For
  each prong number `q`: multigon index, prong index, `k`, `punctured`.
  For each edge number: tail and head prong numbers.  Helpers
  `side_from(q)`, `side_to(q)`, `directions_at(vertex)`.
- New recursion in `coding.cpp` mirroring `recursive_get_weights`: same
  start (monogon 0), same `cycle_edges` order, no descent into uncusped
  monogons.  On first entry into a multigon at prong `pin`, number all its
  prongs `pin, pin+1, ...` (mod `k`) consecutively.  A terminal uncusped
  monogon gets its number when its edge is traversed; monogon 0 gets
  number 0.  Edge orientation = direction of first traversal (from the
  multigon that first reaches the edge; monogon 0's edge points outward).
- Side number `q` = side from prong `q` to the next prong in the
  `cycle_edges` direction (`p <- mod(p+1,k)`), the same direction that
  `t2_pr = mod(t_pr+dir,k)` uses in `fold`.
- The recursion takes its start monogon from the same place as `weights()`
  and asserts (fail-fast) that it equals what `fold_cusp_location` finds.
  It walks the raw structure and must not call the
  `require_normalised`-guarded entry points, because Step 2 runs it on a
  folded, not yet normalised track.
- Accessor `traintrack::numbering() const` (monogon 0); the
  `coding_engine` static takes the start monogon explicitly for Step 2.

## Step 2: correct one-fold map (`traintracks`)

Files: new `include/traintracks/fold_map.hpp` (~80 lines), new
`lib/traintracks/fold_map.cpp` (~200 lines), `include/traintracks/map.hpp`
(replace the body of `fold_traintrack_map`, ~40 lines changed),
`include/traintracks/traintrack.hpp` (declare one member; retire
`fold_infinitesimal_index/generator` or redefine them as "signed side of
the target"), `lib/traintracks/traintrack.cpp` (fix the docstrings; the
header's "global infinitesimal-loop convention" is false).

- `struct fold_map_data`, all in canonical numbers: `moved`, `onto` (old
  edge numbers), `edge_image[e]` (signed new number of each old edge),
  `side` (signed new side number traversed), `moved_first` (word order),
  `prong_image[q]` (new number of each old prong).
- `traintrack::fold_with_map(int f, fold_map_data&)`: compute the pre-fold
  numbering; locate `e0`, `e1`, target, `t_pr`, `t2_pr` via
  `fold_cusp_location` and `target_multigon`; perform the structural fold.
  Preferred correspondence: compute the post-fold numbering *before*
  `normalise()`, between `insert_edge` and `normalise()` in
  `fold(multigon&,...)`, starting from the minimising monogon that
  `minimise_coding` already computes.  Edge and multigon objects are then
  still in place and the pre/post correspondence is by pointer.  Fallback:
  after `normalise()`, match prongs by slot contents (edge-object identity
  survives `swap(multigon&,multigon&)`, ending indices do not), as
  `gates_check.cpp` does.
- `traintracks::fold_traintrack_map(tt0, f)` builds the
  `jlt::freeauto<int>` from `fold_map_data`: unmoved edges map to their
  signed new number; sides map through `prong_image`; the moved edge maps to
  `e0'.S.(+-e1')` or `(+-e1').S.e0'`.  Template signature unchanged, so
  `ttfoldgraph` needs no edit.  Do not pursue the "one fold with distinct
  weights" shortcut: weights collide because `fold` sets `w1 <- w0 + w1`.
- Rewrite the hard-coded word expectations in
  `testsuite/traintracks/test_map_consistency.cpp` and
  `tests/test_ttmap.cpp` in the same commit.
- Regression kept: the abelianisation equals `fold_transition_matrix`.

Outcome of Steps 1 and 2 (2026-09-19):

- `ttnumbering` and `coding_engine::numbering(tt, mono)` live in
  `coding.hpp/.cpp`, with `traintrack::numbering()` as the monogon-0
  wrapper; `ttnumbering` also carries `edge_ptr`, the identity of the edge
  objects, which `fold_with_map` uses to follow edges through the fold.
- `fold_map_data` and `traintrack::fold_with_map` are in
  `fold_map.hpp/.cpp`.  The fallback correspondence was implemented (match
  prongs after `normalise()` by the set of edge objects attached, plus
  multigon type and the multigon's edge set), not the pre-normalise
  variant: it needs no change to `fold(multigon&,...)` and the signature is
  unambiguous on every track tested.  `fold_traintrack_map` (free template
  and member) now builds the map from the record; a fold that cannot be
  performed yields the identity map, as before.
- `fold_infinitesimal_index/generator` were removed (no library callers
  remained); `map_labels.hpp` keeps the range arithmetic and says a side
  index is a canonical prong number.
- New tests: `test_numbering.cpp` (structure, agreement with the
  `weights()` order, copy and re-normalisation invariance, every one-fold
  neighbour of every n=3..5 stratum) and `test_fold_map_paths.cpp` (every
  branch of the n=3, 4 and four n=5 graphs: words continuous in the
  target's numbering, tails and heads mapped correctly, prong images a
  type-preserving bijection, exactly one side letter, matrix agreement;
  random paths and iterated closed paths likewise).
  `test_map_consistency.cpp` rewritten around the fold record with the new
  hand-checked n=3 words (`1 -> 1 -4 2`, path `[1,0]`:
  `2 -> 1 -4 2 -5 -2`); `tests/test_ttmap.cpp` rewritten as its printing
  companion.  All seven fast tests and the slow strata scan pass, so the
  graph's matrices are unchanged.
- Cross-check on the bad cycle: the library's composed map is continuous
  for three iterates, abelianises to the path matrix, keeps the fixed
  monogon's loop fixed and permutes the other five loops in a 5-cycle,
  matching `gates_check.cpp` up to the (now canonical) numbering and
  orientation conventions.
- Consequence for Step 5: `test_fold_map_paths.cpp` already covers the
  planned continuity and endpoint checks; the remaining test work is
  `test_gates.cpp` and `test_ttauto_gates.cpp`.

## Step 3: gates module (`traintracks`)

Files: new `include/traintracks/gates.hpp` (~120 lines), new
`lib/traintracks/gates.cpp` (~300 lines).

- `struct fold_derivative`: `D` on all signed letters (main and side) and
  the turns `T` of the one-step word, recorded at the *target* track's BH
  vertices.  A permutation fold has no turns; a fold onto an unpunctured
  target has one turn between two gates; a fold onto a punctured target has
  two turns at two prong-vertices.  Built from a one-step
  `jlt::freeauto<int>` (at most three letters per image) plus the target's
  `ttnumbering`.
- `class gate_accumulator`: `push_back(const fold_derivative&)` maintains
  the composed `D` and the turn set `T <- T_i union D_i(T)`;
  `analyse(const ttnumbering& start)` closes `T` under `D x D`, builds gates
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

Outcome of Step 3 (2026-09-19):

- `gates.hpp/.cpp` as planned: `fold_derivative` (built from any map plus
  the target numbering, so the word-based overload is the same code fed
  the composed map), `gate_accumulator`, `gate_analysis` with per-vertex
  reports, `analyse_gates`, `bh_vertex_of`.  Gates are computed by
  comparing `D^K` images with `K = (#directions)^2`; directions at a vertex
  come in cyclic order from `ttnumbering::directions_at_prong`, added to
  the numbering for this purpose (`prong_letters`).
- The shape check (Props. 3.3.3-3.3.4) is implemented with the cyclic
  order, not only counts, and reported per vertex as `shape_ok`; it is a
  diagnostic flag, not a fail-fast.
- `testsuite/traintracks/test_gates.cpp` (0.1 s, so it stays in the fast
  set): bad path rejected at exactly the fixed 3-edge monogon, by a refined
  prong, with the shape flag raised; the four controls (seven matching
  branch sequences in all) connected with allowed shapes; accumulator and
  word-based analyses agree; puncture corollary asserted; sweep of all 2630
  closed paths of length <= 5 on the n=5 first stratum: no fail-fast, and
  every connected primitive path has the allowed shape.  No control
  exhibits a refined prong at an unpunctured multigon; the sweeps in Step
  5 remain the guard for that case.

## Step 4: integrate into the search (`ttauto`)

Files: `include/ttauto/folding_path.hpp` (~25 lines),
`include/ttauto/ttauto.hpp` (~40 lines).  `ttfoldgraph.hpp` untouched:
`fold_derivative` is derived from the stored `AMv` words at check time, cost
O(L (n + ninf)) per candidate.

- `folding_path::gates() const`: walk `fp`/`vp`, feed a
  `traintracks::gate_accumulator` from `ttg->traintrack_map(v,f)` and
  `ttg->traintrack(target).numbering()`, analyse at the initial vertex.
- `ttauto::record_pA()`: after the `lambdamax`/`lambdamin` window checks
  and before class insertion, `if (do_check_gates && !p.gates().connected)
  { ++gate_rejected; return; }`.  Setter `check_gates(bool)`; the counter is
  printed with the other totals in `search()`.  `pAclass` unchanged.
  Badwords and `check_all_norms` are unaffected (they act before the
  closed-path test).
- Default: on (decision of 2026-09-19).  Rejections are counted and
  printed; `check_gates(false)` restores the old behaviour.

Outcome of Step 4 (2026-09-19):

- `folding_path::gates()` and the check in `record_pA` as planned;
  rejected candidates go to `rejected_pA_list()` (same `pAclass` keying
  as the accepted list) and are counted by `gate_rejected()`, which is
  cumulative over the whole search because the other statistics counters
  are reset for every initial vertex.  `check_gates(false)` restores the
  old behaviour.  `ttfoldgraph` untouched.
- `testsuite/ttauto/test_ttauto_gates.cpp`: on the n=6 stratum 3,3(2)
  main subautomaton with path length <= 7, 10 classes with the gate test
  off, 9 with it on; the class with dilatation 2.01536 (3 paths) is the one
  rejected and every other class is unchanged.
- Census program `tests/ttauto_gate_census.cpp` (in-place binary,
  ignored): mirrors `examples/ttauto_scan_strata.sh` (main subautomaton,
  same per-stratum path lengths) with the gate test on and lists every
  rejected class.  Results for n = 3..7 (about 8 s):

  | n | stratum | rejected lambda | paths | interpretation |
  |--:|--:|--:|--:|---|
  | 4 | 1 `(2)` | 2.61803, `(x-1)(x^2-3x+1)` | 3 | 3-braid sigma1 sigma2^-1 plus an idle puncture; a genuine class with the same dilatation remains accepted |
  | 6 | 3 `4(2)` | 1.61803, `-(x^2-1)(x^4-3x^2+1)` | 141 (3 classes) | matrix irreducible but not primitive (eigenvalues +-phi); disconnected at the unpunctured 4-gon |
  | 6 | 3 `4(2)` | 1.93185, `-(x^2-1)(x^4-4x^2+1)` | | same: polynomial in x^2, imprimitive, disconnected at the unpunctured 4-gon |
  | 6 | 3 `4(2)` | 2.15372 | | the n=5 stratum 4(1) minimum plus an idle puncture (3-edge monogon) |
  | 6 | 5 `3 3(2)` | 2.01536 | 3 | the known bad path (n=5 stratum 3 3(1) minimum plus an idle puncture) |
  | 6 | 5 `3 3(2)` | 2.54205 | 8 | two length-8 extensions of the bad cycle (search length raised to 8) |
  | 7 | 7 `3 4(2)` | 1.88320 | 6 | the n=6 stratum (4) minimum plus an idle puncture |

  Rejection rate (candidates = closed paths with a primitive matrix
  inside the dilatation window, counted by `gate_candidates()`), after
  the primitivity change: 26 of 11007 candidate paths over n=3..7
  (0.24%), 5 of 93 classes.  Per stratum the rate is 0 except n=4
  stratum 1 (3 of 47, 6.4%), n=6 stratum 3 (6 of 131, 4.6%), n=6 stratum
  5 (11 of 7401, 0.15%) and n=7 stratum 7 (6 of 2208, 0.27%).  Before the
  primitivity change the two imprimitive classes of n=6 stratum 3 reached
  the gate test too: 161 of 11142 paths (1.44%), 7 of 95 classes, 53% on
  that stratum.  The search prints the figure in its statistics block and
  exposes it through `gate_candidates()`, `gate_rejected()` and
  `gate_rejection_rate()`.

  Seven classes, 161 paths, none with an accepted representative.  Every
  rejection is at a prong of a multi-edge punctured monogon split into
  two main gates, except the two imprimitive ones, which are disconnected
  at an unpunctured 4-gon.  Cross-check with `braids.tex` (in this
  directory): with the gate test on, every n=6 row of the scan table
  agrees with the paper's per-stratum table, and every n=7 row except
  s12, where the paper's starred 2.21497 is superseded by 2.02598 at
  search length 10.  The imprimitive classes were never candidates in the
  paper: its enumeration requires a Perron root, and a polynomial in x^2
  has +-lambda.  On the n=5 first stratum, 0 of the 726 primitive
  closed paths of length <= 5 are disconnected.  So the bad path was NOT
  the only spurious pA, but the mechanism is always the same one.
- Consequence for the scan baseline `examples/ttauto_scan_strata.md`
  (three rows change; the old rows are kept struck through with a note,
  and the slow regression ignores such annotations):
  n=6 stratum 3: 1.61803 -> 1.88320; n=6 stratum 5: 2.01536 -> 2.08102
  (search length for this stratum raised from 6 to 8 in the scan script,
  user decision; at length 6 it would read 2.45317); n=7 stratum 7:
  1.88320 -> 2.47541.  Separately, the search length for n=7 stratum 12
  was raised from 8 to 10 (Lizi Guo found a lower-dilatation pA there):
  the row's minimum becomes 2.02598 at length 10, with nothing rejected by
  the gate test on that stratum.  The paper's table row for n=6
  stratum 3 (`ttauto.tex` ~1395, `1.61803`) is one of the entries the
  paper itself marks as below the systole (`\chkabs`); the gate test now
  explains it.
- Primitivity (decided and done 2026-09-19): `descend_graph` now requires
  `TM.is_primitive()` instead of `!TM.is_reducible()`, matching the
  paper's definition of a primitive closed path; the counter and the
  statistics line say "Primitive paths".  The two imprimitive classes of
  n=6 stratum 3 (1.61803, 1.93185; eigenvalues +-lambda) are no longer
  candidates at all, with or without the gate test.  While doing this,
  `jlt::mathmatrix::is_primitive` was found to do one squaring too few
  (it tested A^(2^(pmax-1)) against Wielandt's bound n^2-2n+2), returning
  false for every Wielandt matrix; fixed in `extern/jlt` with a Catch2
  test over n=3..9, and the submodule pointer updated.  Without that fix
  the primitivity rule would have dropped genuine pseudo-Anosovs with a
  large exponent of primitivity.
- The fifth control of the scratch program, `{1,18,65,18,56,57,57,1}`, is
  accepted by the library's gate test (class 2.79497 in
  `test_ttauto_gates`); the earlier failure was the scratch program's own
  branch bookkeeping.  Resolved.
- `tests/test_issue2_bad_path.cpp` removed (its vertex model was wrong and
  `test_gates.cpp` covers the case); `gates_check.cpp` kept as a documented
  diagnostic.

  Minimal example (2026-09-19): the n=4 stratum `(2)` rejection is the
  smallest instance.  Enumerating all closed paths of length <= 5 from
  every vertex, nothing fails for n=3 or on n=4 stratum `3(1)`, and
  nothing of length <= 3 fails for n <= 4; on n=4 `(2)` four paths of
  length 4 and ten of length 5 are primitive but disconnected.  The
  simplest is vertices `[0 2 0 3 0]` (0-based), branches `1 1 2 2`,
  matrix `[[0,0,1],[1,2,0],[0,1,2]]`, polynomial `(x-1)(x^2-3x+1)`,
  dilatation 2.61803: sigma1 sigma2^-1 on three strings plus an idle
  puncture, disconnected at the fixed 3-edge monogon (4 gates, 2
  components), the same pattern as the n=6 bad path.  Pinned in
  `test_gates` and written up in the note (section "A minimal example").

## Step 5: tests

- `testsuite/traintracks/test_numbering.cpp`: numbering tables for n=3..5
  fixtures; counts equal `total_prongs()` and `edges()`; tails and heads are
  valid prong numbers; the numbering of a copy equals the original's; the
  numbering is invariant under `normalise()` of an already-normalised
  track; edge order agrees with the order `weights()` uses.
- `testsuite/traintracks/test_fold_map_paths.cpp`: for every vertex and
  fold of `n=3,trk=0` and `n=4,trk=1`, the one-step word is a continuous
  path in canonical numbers and abelianises to `fold_transition_matrix`.
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

- `CODE_STRUCTURE.md`: subsections for numbering, fold_map, gates; the "one
  fold through the stack" walkthrough now records the target side and
  `record_pA` is gated.
- `testsuite/COVERAGE.md`, `TESTING.md`, `AGENTS.md` repository map,
  `CLAUDE.md` issue-2 section, `devel/macro-status.md` if a flag is added.
- `ISSUE2_BAD_PA_SUMMARY.md`: "resolved by" section.
- Close GitHub issue #2 pointing at the note and the tests.

## Actual scope (measured 2026-09-19 after commit da91acc)

`git diff --stat master...HEAD` (merge base 8add0f3), whole branch
including the coverage pass and the symmetry follow-up:

| Area | Files | Added | Removed | Net |
|---|---:|---:|---:|---:|
| Library (`include/`, `lib/`, CMake) | 15 | 1599 | 378 | +1221 |
| Tests (`testsuite/`, `tests/`) | 17 | 2430 | 433 | +1997 |
| Docs, `devel/`, `examples/`, `.gitignore` | 19 | 3268 | 154 | +3114 |
| Whole repository | 51 | 7155 | 947 | +6208 |

Library size (`include/**/*.hpp` + `lib/**/*.cpp`): 7048 -> 8265 lines,
+17%.  The estimate before Step 1 was 12-13% (~890 lines).  The overshoot
is `gates.cpp/.hpp` (495 lines, ~420 planned), `fold_map.cpp/.hpp` (445,
~280 planned) and `coding.cpp/.hpp` (+417, ~190 planned: follow-up A moved
the weights and cusp walks onto the numbering, and `same_labels` was
added).  The 378 removed lines are mostly the old `fold_traintrack_map`
body in `map.hpp` and the three recursive walks in `traintrack.cpp`.
Tests grew by more than the library, documentation by more than both
(about half of it the LaTeX note).

## Order and commits

One commit per step, build green and `ctest` passing after each.  Steps 1
and 2 land together with the updated word expectations, since they change
the infinitesimal letters of every stored map.  Step 4 is the only
user-visible behaviour change and goes last, with its test.

## Follow-ups

- DONE 2026-09-19: consolidate the four copies of the monogon-0 DFS.
  `ttnumbering` now also records the cusps in fold order (the walk checks
  the entry slot of each multigon first, then the other slots before
  descending, exactly as `recursive_find_cusp` did); `weights(mono)`,
  `weights(iter)`, `fold(f)` and `fold_cusp_location(f)` read the
  numbering, and `recursive_get_weights`, `recursive_set_weights` and
  `recursive_find_cusp` are gone (about 180 lines of
  `lib/traintracks/traintrack.cpp`).  `fold_with_map` reuses its pre-fold
  numbering for the cusp location.  Two walks remain: the coding (normal
  form, both directions) and the numbering (index).  The old walks live on
  as oracles in `test_numbering.cpp`, which compares edge order (from every
  uncusped monogon) and cusp order on every n=3..5 stratum and one-fold
  neighbour; the slow strata scan is byte-identical; census unchanged;
  `ttauto_gate_census 7 7` 7.4 s before and after (search-dominated).
- DONE 2026-09-19: derive the transition matrix from the fold record.
  `fold_map_data::transition_matrix()` builds the permutation-plus-one
  matrix from `edge_image`, `moved` and `onto`; `fold_transition_matrix`
  (free template and member) folds once and reads it; the unit-weight
  construction (n folds) moved to `testsuite/oracles.hpp` and is compared
  with the record for every fold index of every vertex of the n=3, 4, 5
  test graphs.  `check_fold_map_main_transition` removed.  Graph
  construction was already one fold per branch since Step 2, so no
  further speed-up was expected or measured.
- (superseded) Derive the transition matrix from the fold record.
  `fold_transition_matrix` folds `n` copies of the track with unit weights,
  and `fold_traintrack_map` calls it, so building a graph costs `n+1` folds
  per branch.  After Step 2 the fold record knows the permutation and the
  `+1` entry directly; one fold per branch suffices, and the old routine
  becomes a test oracle only.  An `n`-fold speed-up of `ttfoldgraph`
  construction, relevant to issue #14 (large graphs).

## Risks

Rewritten 2026-09-19 after the coverage-and-risks pass (see below); the
original four items are kept with their outcome.

- Gates finer than prongs at unpunctured multigons make the BH condition
  stricter than "all sides realised".  **Measured, never bites so far**: of
  the 113 stored representatives of every accepted class on every stratum
  for n=3..7 (census) and of the 726 primitive connected closed paths of
  length <= 5 in the n=5 first-stratum automaton (`test_gates`), none has
  a refined prong at an unpunctured multigon.  Every rejection so far is a
  refined prong at a puncture (Corollary "puncture" of the note).  The
  census prints the count, so a stratum where the stricter condition
  starts to bite will show up; those paths are the ones to confirm with
  Trains.
- Automorphic vertices.  **Verified**: `test_fold_map_paths` now asserts
  that the graphs it walks contain vertices with nontrivial cyclic symmetry
  (n=4 stratum 1: 2 of 4 vertices; n=6 stratum 3,3(2): 10 of 110) and
  reflection symmetric vertices, and every branch at them passes the
  continuity, endpoint and abelianisation checks.  The earlier comment
  that the n=5 strata contain cyclically symmetric vertices was wrong;
  they only contain reflection symmetric ones.  Those checks are
  automorphism-invariant, so `test_numbering` (2026-09-19, later the same
  day) pins the actual claim with `ttnumbering::same_labels`: over the
  n=4 s1, n=5 s1, n=6 s3 and s5 automata, every fold result carries the
  same labels as the track rebuilt from its coding and as the stored
  target vertex, and at every cyclically symmetric vertex exactly
  `cyclic_symmetry().order()` start monogons give the same labels.  The
  labels are a function of the coding; the identification by number is the
  identity on labels, not merely an isomorphism up to the automorphism.
  `ttnumbering::operator==` also compares `mgv` bookkeeping and is only
  for copies of one object.
- Redundancy (resolved): the numbering was briefly a fourth copy of the
  monogon-0 DFS; the weights and cusp walks now read it, and
  `test_numbering.cpp` keeps the old walks as oracles.
- Published numbers changed (four entries of the ttauto paper, one of the
  braids paper appendix); recorded in `pubs/braids/ERRATA.md` and pending
  for the ttauto paper.
- Norm-bounded search mode (new).  It is the mode the literature
  comparison (`ttauto_min_example`) relies on and it had no test.  **Now
  covered** by `testsuite/ttauto/test_min_dilatations.cpp`: n=3, 4, 5
  minima with polynomials, and the n=4 stratum-1 reducible candidate
  ((x-1)(x^2-3x+1), dilatation 2.61803) that only the gate test removes.
- Dependency bugs (new).  The `is_primitive` off-by-one in jlt showed that
  the acceptance chain trusts jlt for irreducibility, primitivity, the
  characteristic polynomial and the Perron root.  **Partly covered**:
  `test_ttauto_search` checks every stored representative's class
  polynomial against an independent principal-minors expansion and its
  dilatation against power iteration on the matrix; `is_reducible` was
  brute-force checked once (24022 matrices); `is_primitive` has Wielandt
  tests in jlt.  Not covered: jlt's own `polynomial` arithmetic beyond
  evaluation.
- Debug builds (new, found by the coverage build): `edge::nends` was a
  non-inline `static const int` and the testsuite did not link at `-O0`.
  Fixed (`static constexpr`).  The testsuite now builds and passes in
  Debug as well as Release; nothing in CI enforces that.

## Coverage and risks pass (done 2026-09-19)

Outcome of the pass proposed the same day.  One commit; library changes
limited to the `edge::nends` fix and a read-only `pAclass::paths()`
accessor.

- Measured line coverage before and after with gcov (table at the top of
  `testsuite/COVERAGE.md`; procedure in `TESTING.md`): 80% -> 86% of the
  2416 instrumented library lines; `traintrack.cpp` 64% -> 88%,
  `multigon.cpp` 60% -> 74%, `ttauto.hpp` 78% -> 89%.  What remains is
  mostly fail-fast exits, diagnostic printers, the broken string
  constructor and builders the automaton does not use.
- `test_traintrack_core`: every stratum for n=3..5 and every one-fold
  neighbour (29 tracks); idempotent normalisation on coding and numbering,
  coding round trip, mirror image from the reversed coding, labels,
  `pure_braid`, printers, `print_coding` parsed back.  The string
  constructor is documented as broken instead of tested.
- `test_mathmatrix_permplus1`: 88 one-fold matrices of the n=4 and n=5
  automata against the unit-weight oracle, dense round trip, `order()` by
  dense powers, products both sides, printers.
- `test_folding_path_and_badwords`: badwords content (each entry is the
  square of a pattern-idempotent closed path; the length-1 layer is
  complete; the table is deterministic, 524 bad words) and the path
  algebra (`subpath`, `operator*`, `ending_equals`, `cycle_path`).
- `test_ttauto_search`: exact eight classes on n=5 stratum 1 at length 8
  with polynomials and lengths, 3555 candidates, 0 rejected, oracles on
  all 26 representatives.
- New `test_min_dilatations` (norm-bounded mode), see Risks.
- `test_fold_map_paths` walks n=4 stratum 1 and n=6 stratum 5 as well and
  asserts automorphic vertices are present (12 cyclic, 24 reflection).
- `test_gates` and the census count refined prongs at unpunctured
  multigons among accepted paths: 0 of 726 and 0 of 113.
- Oracles added to `testsuite/oracles.hpp`: `charpoly_by_minors` (Bareiss
  principal minors) and `power_iteration_radius`.
- Census timing with the extra statistic: 8.7 s for n=3..7 (was 7.4 s).
