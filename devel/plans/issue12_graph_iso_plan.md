# Issue #12 Graph-Isotopy Canonicalization Plan

## Goal

Implement a structural, orientation-preserving isotopy oracle/canonicalizer for
`traintrack` that does not rely on existing coding-numbering conventions, then
use it to resolve the 29/71 equivalence case from issue #12.

## Scope for Phase 1

1. Preserve orientation (no reflection equivalence in this phase).
2. Ignore puncture flags and multigon labels in matching/canonicalization.
3. Treat tracks as topological graphs with prong incidence preserved by gadget
   expansion.
4. Keep existing legacy `coding()` untouched.

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
