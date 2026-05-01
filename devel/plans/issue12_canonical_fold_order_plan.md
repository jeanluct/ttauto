# Issue #12 Canonical Fold/Cusp Ordering Plan

## Restart Snapshot (2026-05-01)

This file now doubles as the fold-ordering restart checklist.

### What is implemented

1. Canonical witness API in graph-iso layer:
   - `canonical_witness canonical_witness_oriented(const traintrack&)`
2. Canonical fold/cusp APIs in `traintrack`:
   - `fold_cusp_location_canonical(...)`
   - `fold_infinitesimal_index_canonical(...)`
   - `fold_infinitesimal_generator_canonical(...)`
3. Equality mode switch does not disable canonical fold APIs.
   - They are available in both operator== modes.
4. Regression coverage includes:
   - `test_graph_iso_canonical_witness`
   - `test_issue12_coding_uniqueness`
   - `test_canonical_fold_map_consistency` (new)

### Commits to anchor current state

- `234fd63` add canonical witness and canonical fold-index APIs
- `d8dbb99` prune witness search by signature classes
- `3a5c93e` compile-time equality mode switch (graph-iso vs legacy coding)
- `a327805` expanded in-source docs for switch/tests
- `219367c` canonical fold map consistency regression

### Known caveat

The new fold-map consistency regression currently uses exhaustive matrix
relabel-canonicalization and takes ~30s. This is acceptable for now but is the
top optimization candidate for next pass.

### Exact restart commands

Default (graph-iso mode):

1. `cmake -S . -B build`
2. `cmake --build build -j`
3. `ctest --test-dir build --output-on-failure`

Legacy mode:

1. `cmake -S . -B build-legacy -DTTAUTO_USE_GRAPH_ISO_EQUALITY=OFF`
2. `cmake --build build-legacy -j`
3. `ctest --test-dir build-legacy --output-on-failure`

Expected result at snapshot: `9/9` tests pass in both modes.

## Goal

Restore a canonical fold/cusp ordering under the new oriented graph-isotopy
equality backend, with semantics aligned to strict label/puncture-aware
equivalence classes.

## Agreed Semantics

1. Orientation-preserving only.
2. Preserve multigon labels.
3. Preserve puncture flags.
4. Canonical fold/cusp order must be canonical within the same equality class
   used by `operator==`.

## Strategy

### 1) Canonical witness in `graph_iso`

Add a deterministic witness structure:

- canonical multigon order,
- inverse rank map,
- per-multigon prong shift.

Add API:

- `canonical_witness canonical_witness_oriented(const traintrack&)`.

The witness is selected by lexicographic minimization of a structural
certificate built from relabeled endpoint multiplicities and multigon
attributes (prongs, puncture, label).

### 2) Canonical cusp keys from witness

For each original cusp location `(m,p,e)` build canonical key:

- `m_can = witness.multigon_rank[m]`
- `p_can = (p + witness.prong_shift[m]) mod prongs(m)`
- `e` unchanged

Sort by `(m_can, p_can, e)` to define canonical cusp order.

### 3) Canonical fold APIs in `traintrack`

Add public APIs (non-breaking, alongside existing ones):

- `fold_cusp_location_canonical(...)`
- `fold_infinitesimal_index_canonical(...)`
- `fold_infinitesimal_generator_canonical(...)`

Map canonical fold index back to original fold index using the canonical cusp
order permutation, preserving clockwise/counterclockwise parity semantics.

## Tie-break Rule

If structural certificates tie, choose the lexicographically smallest original
multigon order vector. This avoids search-order dependence.

## Validation

1. Existing deterministic tests remain green.
2. Equivalent tracks (e.g. issue pair) produce identical canonical cusp/fold
   ordering.
3. Relabeled tracks are no longer equivalent under strict mode and should not
   share canonical fold order.

## Next Action Queue (ordered)

1. Replace exhaustive canonicalization in
   `test_canonical_fold_map_consistency.cpp` with a faster structural
   fingerprint based on `mathmatrix_permplus1` sparse data.
2. Add one additional negative-control pair where canonical fold-generator
   sequences differ (to ensure test would fail on wrong witness/permutation).
3. If runtime is reduced significantly, include this test in any future
   expanded CI matrix by default.
