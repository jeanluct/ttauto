# Issue #12 Graph-Isotopy Canonicalization Plan

## Restart Snapshot (2026-05-01)

This section is the authoritative quick-restart state for future work.

### Current behavior in tree

1. `traintrack::operator==` is now compile-time selectable:
   - graph-iso mode (default): orientation-preserving structural matcher
   - legacy mode: canonical-coding equality
2. Build flag:
   - CMake option: `TTAUTO_USE_GRAPH_ISO_EQUALITY` (default `ON`)
   - public compile define: `TRAINTRACKS_USE_GRAPH_ISO_EQUALITY`
3. Equality-mode switch commit history (latest relevant):
   - `3a5c93e` add equality mode switch and mode-aware tests
   - `a327805` add explanatory comments/docs in touched source/tests
4. Canonical fold-map regression was added:
   - `219367c` `testsuite/traintracks/test_canonical_fold_map_consistency.cpp`
5. Issue tracker updates posted:
   - `https://github.com/jeanluct/ttauto/issues/12#issuecomment-4360080783`
   - `https://github.com/jeanluct/ttauto/issues/12#issuecomment-4360174783`

### Validation matrix at snapshot

Both modes pass full testsuite deterministically:

- graph-iso mode:
  - configure/build: `cmake -S . -B build && cmake --build build -j`
  - test: `ctest --test-dir build --output-on-failure`
  - result: `9/9` pass
- legacy mode:
  - configure/build: `cmake -S . -B build-legacy -DTTAUTO_USE_GRAPH_ISO_EQUALITY=OFF && cmake --build build-legacy -j`
  - test: `ctest --test-dir build-legacy --output-on-failure`
  - result: `9/9` pass

### Known runtime/perf note

`testsuite_test_canonical_fold_map_consistency` currently takes ~30 seconds in
each mode. The test is deterministic and stable, but uses exhaustive matrix
relabel-canonicalization for robustness.

### Next technical target from this snapshot

Optimize canonical fold-map cross-representation comparison so the new test keeps
the same semantic guarantee with lower runtime:

1. Prefer exploiting `mathmatrix_permplus1` structure directly (single
   permutation + optional plus-one location) instead of full bi-permutation
   search over dense matrices.
2. Build a canonical fingerprint from sparse decode tuples
   `(row_perm, plus1_row, plus1_col)` under admissible relabel actions.
3. Keep existing test logic and acceptance criteria unchanged; only replace the
   fingerprint implementation.

Status update (2026-05-01): completed in `1e8676a`.

## Release Gate Checklist (Graph-Iso Default)

This is the short sign-off matrix to rerun before closing issue #12 work on
this branch. The goal is to validate default production behavior, not to expand
legacy-mode coverage.

### 0) Clean configure/build

```bash
cmake -S . -B build
cmake --build build -j
```

Expected:

- configure/generate succeeds
- build completes with no link errors
- `testsuite/test_canonical_fold_map_consistency` builds

### 1) Deterministic testsuite gate (primary)

```bash
ctest --test-dir build --output-on-failure
```

Expected:

- all tests pass (`9/9` at current snapshot)
- includes:
  - `testsuite_test_issue12_coding_uniqueness`
  - `testsuite_test_graph_iso_canonical_witness`
  - `testsuite_test_canonical_fold_map_consistency`

### 2) Fast repo smoke binaries (recommended by AGENTS)

```bash
./tests/test_permplus1
./tests/test_traintrack
./tests/test_folding_path
./tests/test_badwords
./examples/ttauto_torus
```

Expected:

- all commands exit `0`
- no crash/assert in release build

Note: `./examples/ttauto_min_example` has historically been noisy/long in some
sessions; use it manually when needed, but keep gate time bounded.

### 3) Issue #12 behavior sanity spot-check

Run:

```bash
./build/testsuite/test_issue12_coding_uniqueness
```

Expected:

- exit `0`
- confirms structural equality for the 29/71 pair and graph-iso-level
  deduplication behavior assumptions used by the regression.

### 4) (Optional) legacy guard while flag still exists

Only while `TTAUTO_USE_GRAPH_ISO_EQUALITY` remains in tree:

```bash
cmake -S . -B build-legacy -DTTAUTO_USE_GRAPH_ISO_EQUALITY=OFF
cmake --build build-legacy -j
ctest --test-dir build-legacy --output-on-failure
```

Expected:

- all tests pass in compatibility mode

When legacy mode is removed from this branch, remove this step entirely.

## Pre-close TODOs

Before declaring production-final on this branch:

1. Add at least one negative-control canonical-fold test where non-isotopic
   tracks are required to disagree in canonical fold-generator sequence.
2. Remove legacy equality mode if/when approved (switch, macro, and mode-guard
   branches in tests).
3. Re-run the release gate checklist after legacy removal to ensure no hidden
   dependency on compatibility paths.

## Goal

Implement a structural, orientation-preserving isotopy oracle/canonicalizer for
`traintrack` that does not rely on existing coding-numbering conventions, then
use it to resolve the 29/71 equivalence case from issue #12.

## Scope for Phase 1 (historical)

1. Preserve orientation (no reflection equivalence in this phase).
2. Ignore puncture flags and multigon labels in matching/canonicalization.
3. Treat tracks as topological graphs with prong incidence preserved by gadget
   expansion.
4. Keep existing legacy `coding()` untouched.

Status: completed and then tightened. Current semantics are stricter than this
historical initial scope (puncture and label preservation are enforced).

## Core Representation

### Why

A single-vertex-per-multigon model loses which train-track edges attach to
which prong. We need prong-distinguishing structure while still reducing to a
plain graph-style canonicalization problem.

### Blow-up model

For each multigon with `k` prongs:

- `k >= 3`: create `k` boundary nodes connected in a peripheral `k`-cycle.
- `k == 2`: use a fixed asymmetric bigon gadget with two attachment ports.
- `k == 1`: use a fixed asymmetric monogon gadget with one attachment port.

Each train-track branch edge is then represented by an edge connecting the two
attachment-port nodes corresponding to its two incident prongs.

This preserves:

- prong incidence,
- edge multiplicity,
- topological adjacency independent of numbering.

## Orientation-Preserving Constraint

We must forbid local orientation reversal around expanded multigon gadgets.

Implementation strategy:

1. Store per-gadget cyclic successor metadata (clockwise order).
2. During matching/canonicalization, only permit automorphism candidates that
   preserve this successor relation.

This keeps reflection handling separate (already represented at `ttfoldgraph`
symmetry layer).

## Implementation Steps

### Step 1 - New graph-iso module (internal first)

Add new files:

- `include/traintracks/graph_iso.hpp`
- `lib/traintracks/graph_iso.cpp`

Define internal data structures:

- node type / edge type enums,
- blown-up graph container,
- metadata for gadget ports and cyclic order constraints.

Add builder:

- `build_iso_graph(const traintrack&)`.

### Step 2 - Isotopy matcher + canonical key

Implement orientation-preserving backtracking matcher on blown-up graphs:

- `bool is_isotopic_oriented(const traintrack&, const traintrack&)`.

Implement canonical key/signature generator from the same structure:

- `traintrack::intVec canonical_key_oriented(const traintrack&)` (or string key
  if easier during prototype, then map to `intVec`).

### Step 3 - Integration path (safe rollout)

Initial rollout (recommended):

1. Add a temporary experimental API call path for issue #12 tests only.
2. Compare results vs current `canonical_coding()`.

If results are stable and correct:

1. Switch `traintrack::canonical_coding()` wrapper to the new key backend.
2. Keep current coding module available for regression comparison until issue is
   closed.

## Test Plan

Add/adjust tests under `testsuite/traintracks/`:

1. Issue #12 pair:
   - 29/71 should be equal under the new oriented-isotopy canonical key if this
     is the intended relation.
2. Negative controls:
   - clearly non-isomorphic tracks must remain non-equal.
3. Numbering invariance:
   - prong/multigon renumbering variants of the same track must match.
4. Orientation guard:
   - reflected-only equivalents should not match in orientation-preserving mode.

Run at minimum:

- `testsuite_test_issue12_coding_uniqueness`
- `testsuite_test_traintrack_core`
- `testsuite_test_map_consistency`
- `testsuite_test_folding_path_and_badwords`
- `testsuite_test_ttauto_search`

## Risks and Mitigations

1. Combinatorial blow-up for naive automorphism search.
   - Mitigate with early partitioning by local invariants (degree, gadget role,
     incident edge-type counts, cycle position class).
2. Ambiguity in monogon/bigon special gadgets.
   - Keep gadgets asymmetric and fixed by construction; add targeted tests.
3. Semantic drift while issue is unresolved.
   - Roll out behind experimental path first; switch default only after test
     agreement.

## Success Criteria

1. 29/71 behavior is explained and deterministic under explicit
   orientation-preserving graph-isotopy semantics.
2. Canonical key is numbering-independent and stable across rebuilds.
3. Existing core deterministic tests remain green.
4. Reflection remains a separate notion unless explicitly enabled.

Status at snapshot: met for current branch state.
