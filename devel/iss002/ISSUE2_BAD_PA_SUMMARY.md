# Issue #2 - "Bad pA" Repro Summary

## Problem

Issue #2 is that ttauto can very rarely report a pseudo-Anosov candidate that
is actually reducible. The known canonical case is in `devel/iss002/`:

- `n=6_5_bad_tt_data.m`
- `n=6_5_bad_tt.nb`
- `n=6_5_bad_tt.train`
- `toby_hall-email_2009-04-30.pdf`

The bad cycle in the 90-vertex graph is:

- 1-based: `{29, 46, 43, 71, 88, 85, 29}`
- 0-based (C++): `{28, 45, 42, 70, 87, 84, 28}`

## Unified evidence and interpretation

The issue tracker and Toby Hall's email tell the same story from two angles:

- The tracker (`#2`) records the symptom: occasional spurious pA detections.
- The trusted external run (`n=6_5_bad_tt.train`) for this braid reports
  growth `2.01536`, entropy `0.700796`, and isotopy class **Reducible**.
- Toby's explanation: irreducible/primitive transition-matrix behavior is not,
  by itself, enough to certify pA. For an efficient graph map, the gate
  connectivity condition via infinitesimal edges must also hold (BH criterion,
  cited in email as Prop. 3.3.2).
- In this specific example, gate connectivity fails (noted at vertices 2 and
  6), so reducibility is expected even though the dilatation data looks pA-like.
- Toby also notes deleting the 4th string gives a 5-braid `1 2 -4 -3` that is
  genuinely pA with the same growth `2.01536`, reinforcing that growth alone is
  not a discriminator here.

So the real gap is not matrix assembly; it is missing gate/reducibility checks
on top of the matrix criterion.

## Current reproducer status

`tests/test_issue2_bad_path.cpp` now provides a direct reproducer that:

1. builds the `n=6` automata,
2. locates the 90-vertex subgraph containing the bad cycle,
3. follows the bad cycle and composes the train-track map,
4. prints map and matrix, and
5. verifies `matrix(path) == matrix(from_map)`.

The reproduced transition matrix has spectral radius `2.015357...`, matching
the trusted growth value to displayed precision.

Additional diagnostics now in the reproducer:

- Hardwired bad realization on the known graph:
  - `n=6`, `trk=4`, `sgidx=0`
  - cycle `{29,46,43,71,88,85,29}`
  - branch sequence `{1,0,3,2,1,2}` (0-based)
- Primitive-matrix confirmation is explicit in output.
- Derivative map printout includes the key example:
  - `D(1) = -4`, `D(-1) = -2`

## Where we stand now

We implemented an initial gate pipeline in `tests/test_issue2_bad_path.cpp`:

1. per-vertex directed-generator inventory,
2. gate partitioning by eventual equality under `D`,
3. infinitesimal-edge-based connectivity check.

This currently reports broad non-connectivity and is likely **too strict / not
the same vertex model** as the trusted `.train` output. In particular:

- `.train` vertex 1 has edge directions `1 -1 12` and gate-links
  `-1<->12`, `1<->12`, which is connected.
- our current `(multigon,prong)`-based model does not yet reproduce that local
  picture.

So the next correction is conceptual/modeling, not numerical:

- distinguish clearly among **main edges**, **peripheral edges**, and
  **infinitesimal joins**;
- align our local-vertex/gate objects with the `.train` interpretation before
  treating connectivity output as authoritative.

## What this implies for implementation

To close issue #2 robustly, pseudo-Anosov classification should require both:

1. existing matrix criterion, and
2. gate/infinitesimal-edge connectivity criterion.

This aligns with the tracker notes (post-issue-3) that gate detection is the
next substantive step.

## Next steps

- Lock the local vertex model so gate sets and infinitesimal joins match the
  `.train` per-vertex format (especially vertices 1, 2, and 6).
- Re-run gate connectivity under that corrected model and confirm the expected
  failures are localized (not global).
- Only then add regression assertions for issue #2 classification.

## Work log (latest session)

Completed in code (`tests/test_issue2_bad_path.cpp`):

- Added derivative map printout for signed main generators.
- Added primitive-matrix assertion on the bad path matrix.
- Added per-vertex directed-generator inventory including infinitesimal labels.
- Added gate partition computation using eventual equality under `D^k`.
- Added a first gate-connectivity implementation based on infinitesimal joins.

Key observed outputs on the bad case:

- Matrix remains primitive (`yes`) while isotopy class is known reducible.
- Derivative data reproduces expected examples such as
  `D(1) = -4`, `D(-1) = -2`.
- The first connectivity implementation flags many/all vertices as
  non-connected, which conflicts with `.train` where only a subset fail.

Interpretation of current mismatch:

- The implemented local vertex model (`(multigon,prong)`) likely does not yet
  match the switch-level vertex model used in `.train` gate output.
- Also, the `.train` wording indicates infinitesimal edges act as connectors
  **between gate elements**, not as standalone gate elements to be treated the
  same way as main/peripheral directions.

Current status summary:

- Reproducer is strong and stable for the bad path and derivative diagnostics.
- Gate pipeline is partially in place but not yet aligned with `.train`
  semantics, so connectivity results are provisional and should not yet be used
  as classification criteria.

## Known-good run snapshot

- Branch: `iss002-spurious-pAs`
- Commit anchor: `ab300f2`
- Command:
  - `scons tests/test_issue2_bad_path && ./tests/test_issue2_bad_path`
- Key output checks:
  - `Primitive transition matrix: yes`
  - `D(1) = -4, D(-1) = -2`
  - `Directed generators by vertex (multigon,prong):`
  - `vertex (6,1): -13 -1 2 5 13`
  - `Gates by vertex (multigon,prong) on all generators (main+infinitesimal):`
  - `vertex (6,1): {-13}, {-1 2}, {5}, {13}`
