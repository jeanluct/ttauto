# Issue #12 Canonical Fold/Cusp Ordering Plan

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
