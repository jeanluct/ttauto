# Issue #3 Plan - Train-Track Map Construction

This plan is for issue #3 ("find the train track map") on branch
`iss003-tt-map`.

## Goal

Produce a correct train-track map by accumulating per-fold maps during
`traintrack::fold` workflows and composing them along folding paths.

## Key decisions

1. Source of truth: fold-time composition (not matrix inference).
2. Orientation convention: signed generators (`+e`, `-e`).
3. Delivery strategy:
   - Phase A: main edges first (MVP, correctness gate)
   - Phase B: add infinitesimal edges with stable labeling
4. Validation strategy: both automated tests and manual matrix/map checks.

## What stood out during repo review

- `include/traintracks_util.hpp` currently contains placeholder logic:
  - `int infinitesimal = -ngen; // placeholder`
  - This confirms issue #3 is partially scaffolded and not finalized.
- `tests/test_freeword.cpp` exists and appears directly relevant to issue #3,
  but it is not currently part of default SCons test targets.
- Existing comments already identify ambiguity from matrix-based reconstruction
  (especially orientation and infinitesimal-edge handling).
- Notes/paper artifacts exist in:
  - `devel/iss003/ttmaps_from_automata.pdf`
  - `ttmaps_from_automata.pdf`
  These should guide expected composition behavior and notation.

## Implementation phases

### Phase 0a - freeword suitability review (design checkpoint)

1. Review `freeword`/`free_auto` API and semantics before map work:
   - composition order,
   - inverse/sign conventions,
   - generator indexing and bounds expectations,
   - container inheritance tradeoffs.
2. Decide whether current design is suitable as-is for issue #3 MVP, or if a
   minimal interface cleanup is required first.
3. Record that decision in this file (or issue comment) before implementation.

Deliverable: explicit suitability decision and constraints for `freeword`.

### Phase 0b - freeword strengthening (pre-map hardening)

1. Make `include/freeword.hpp` self-contained (explicit includes, no transitive
   include assumptions).
2. Fix obvious correctness hazards (e.g., constructor/shadowing/state setup).
3. Add/extend focused tests in `tests/test_freeword.cpp`.
4. Add `tests/test_freeword.cpp` to `tests/SConscript` so it runs with normal
   test workflows.

Deliverable: stable and test-backed `freeword` foundation.

### Phase 0 - Baseline and wiring

1. Add `tests/test_freeword.cpp` to `tests/SConscript`.
2. Build and run it as-is; capture current behavior and any failures.
3. Keep this baseline output as a reference while refactoring.

Deliverable: reproducible baseline for map-related tests.

### Phase 1 - Semantics lock-in

1. Document map semantics in code comments near fold/map methods:
   - composition order,
   - generator numbering,
   - orientation sign convention,
   - what is included in MVP (main edges only).
2. Ensure map composition location is explicit (prefer folding path composition
   using per-step maps from fold operations).

Deliverable: one unambiguous in-code contract.

### Phase 2 - Main-edge map correctness (MVP)

1. Make each fold return/apply the correct small map on main-edge generators.
2. Compose these small maps along a path and compare with transition-matrix
   action on main edges.
3. Add tests:
   - identity path -> identity map,
   - short known fold sequence -> expected images,
   - sign/orientation behavior for inverses.

Deliverable: main-edge mapping validated by tests and matrix agreement.

### Phase 3 - Infinitesimal-edge support

1. Replace placeholder infinitesimal generator logic.
2. Introduce deterministic infinitesimal labels based on traintrack-local
   structure (e.g., multigon index + prong index), consistent under folding.
3. Extend per-fold small maps to include infinitesimal-edge images.
4. Add tests for infinitesimal consistency through composed folds.

Deliverable: composed map includes both main and infinitesimal edges.

### Phase 4 - Integration cleanup

1. Remove or de-prioritize matrix-guessing paths for map recovery where
   fold-time composition now provides canonical maps.
2. Keep matrix machinery for numeric checks, but not as primary map source.
3. Improve diagnostics for map mismatches in test/debug paths.

Deliverable: one primary path for map construction, one for verification.

### Phase 5 - Validation and closeout

1. Run recommended matrix from `TESTING.md`.
2. Run map tests including `test_freeword` and any new focused tests.
3. Run one manual sanity scenario (print matrix + map for known path).
4. Post issue update with:
   - what is complete,
   - what is intentionally deferred.

Deliverable: issue #3 either closed or split into clearly scoped follow-ups.

## Constraints and guardrails

- Preserve current behavior unless the change is directly part of map logic.
- Avoid broad ownership rewrites; keep changes incremental and test-first.
- Keep commits small and scoped by phase.
- Prefer deterministic labeling and explicit composition over inferred heuristics.

## Risks to watch

- Composition order mistakes (left-vs-right action confusion).
- Orientation sign mismatches that still satisfy matrix checks accidentally.
- Infinitesimal-label instability across folds and normalisation.

## Immediate next step

Start Phase 0: wire `test_freeword` into default tests and establish baseline
behavior before touching map construction logic.

## Progress notes

- Added `ttmap_labeler` for unified main/infinitesimal generator numbering.
- Added matrix-from-map consistency checks in `test_freeword` and aligned
  transpose conventions.
- Replaced the old infinitesimal placeholder selector with a geometry-based
  index lookup:
  - `traintrack::fold_infinitesimal_index(f)`
  - resolves cusp index using existing fold-order traversal
  - maps cusp prong to global peripheral generator index
- Added test coverage that `fold_infinitesimal_index(f)` stays in range for all
  fold indices in the two `test_freeword` scenarios.
- Added per-step map/matrix consistency assertions across reachable fold maps
  in the two `test_freeword` scenarios.
- Important caveat: this is a first geometry-based implementation, not the
  final semantic proof for infinitesimal orientation/sign under all folds.
