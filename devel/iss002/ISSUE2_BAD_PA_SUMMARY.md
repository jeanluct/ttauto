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

## What this implies for implementation

To close issue #2 robustly, pseudo-Anosov classification should require both:

1. existing matrix criterion, and
2. gate/infinitesimal-edge connectivity criterion.

This aligns with the tracker notes (post-issue-3) that gate detection is the
next substantive step.

## Next steps

- Make the reproducer branch-disambiguation match the exact bad example path
  data from `n=6_5_bad_tt_data.m` when multiple parallel branches share the
  same source/target vertices.
- Add gate extraction on this case and confirm the expected non-connectivity.
- Add a regression check so this path cannot be classified as pA in search
  output.
