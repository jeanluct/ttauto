# Issue #2 Detailed Plan - Gates Including Infinitesimal Edges

## Objective

Close issue #2 (spurious pA detections) by implementing a gate-based
reducibility check that augments the existing transition-matrix criterion.

The bad reference case is the hardwired path in `tests/test_issue2_bad_path.cpp`
for `n=6`, `trk=4`, `sgidx=0`, cycle `{29,46,43,71,88,85,29}`.

## Core definition to implement

At each vertex (identified as a multigon-prong pair `(m,p)`), consider all
directed edges that start at that vertex. A gate is an equivalence class under

- `a ~ b` iff there exists `k > 0` such that `D^k(a) = D^k(b)`.

Here `D` is the derivative map (first directed symbol in image path).

This includes comparisons between:

- main directions,
- infinitesimal directions,
- inverses of incoming directions (since those also start at the same vertex).

## Clarifying what "gates joined by infinitesimal edges" means

This phrase is easy to misread, so here is the intended graph-theoretic meaning
used for implementation:

1. Fix one vertex `v = (multigon,prong)`.
2. Compute gates at `v` as equivalence classes of directed generators starting
   at `v`, using `a ~ b` iff `D^k(a)=D^k(b)` for some `k>0`.
3. Build a **gate graph** whose nodes are those gates.
4. Add an edge between two gate nodes when an infinitesimal edge has endpoint
   directions lying in those two gates.

So:

- We do **not** require every pair of gates to be directly linked.
- We do require the gate graph at the vertex to be connected (path-connected).
- Singleton gate sets are allowed; a vertex with one gate is trivially
  connected.

Concrete pattern:

- Gates at a vertex: `{A}`, `{B}`, `{C}`.
- If infinitesimal links give only `{A}-{B}`, then `{C}` is isolated and the
  connectivity condition fails.

## Why this plan

Issue #2 evidence (tracker + Toby Hall notes + `.train` output) says matrix
primitivity is not sufficient; gate/infinitesimal connectivity is required.
So we need a canonical, traintrack-internal gate computation independent of
external ad hoc edge numbering.

## Canonical labeling strategy (internal)

Do **not** use `.train` numbering as source of truth.

Use canonical symbols derived from traintrack structure:

1. Main generators: existing signed labels `+/-1..+/-nmain`.
2. Infinitesimal generators: existing signed labels from `ttmap_labeler`
   (`+/- (nmain+1 .. nmain+ninf)`) with fixed orientation convention.
3. Vertex identifier: `(multigon index, prong index)`.

This gives deterministic labels across runs and branches.

## Data model to add (Phase A)

In the issue-2 reproducer first (then promote to library helpers if stable):

1. Build `start_vertex(g)` for any signed generator `g` (main and
   infinitesimal), returning `(m,p)`.
2. Build `outgoing_generators(vertex)` for each vertex, including:
   - outgoing main directions,
   - outgoing infinitesimal direction,
   - inverse of incoming infinitesimal direction if represented separately.
3. Print per-vertex directed generator sets for sanity.

Deliverable: deterministic per-vertex directed alphabet.

## Derivative map extension (Phase B)

Current reproducer computes `D` for main generators only.

Extend to all directed generators:

1. Main part:
   - `D(g)` = first symbol of `AM.get_action(g)` (already implemented).
2. Infinitesimal part:
   - derive first-symbol action from train-track map convention used in
     `fold_traintrack_map` composition,
   - ensure orientation/sign convention is consistent with issue #3 choices.
3. Add consistency checks:
   - every `D(g)` is valid nonzero signed generator in domain,
   - iteration stays in domain.

Deliverable: total derivative map on the full directed alphabet.

## Gate computation per vertex (Phase C)

For each vertex `(m,p)`:

1. Let `S(m,p)` be directed generators starting there.
2. Build equivalence classes on `S(m,p)` with criterion
   `a ~ b` iff `D^k(a)=D^k(b)` for some `k>0`.
3. Use bounded iteration on finite state space (safe bound from domain size).
4. Output full partition including singleton gates.

Deliverable: machine-computed gate partition at each vertex.

## Connectivity check (Phase D)

Implement the BH-style local connectivity test from gate data:

1. Build local infinitesimal graph at each vertex:
   - nodes = gates at that vertex,
   - edges = infinitesimal adjacency relations from traintrack geometry.
2. Check required connectivity condition per vertex.
3. Aggregate to a global pass/fail reducibility indicator for this criterion.

Operational implementation detail in current reproducer:

- Build one global graph with nodes `(vertex, gate_id)`.
- Add one undirected edge per infinitesimal generator connecting the node for
  `(start vertex of +i, gate containing +i)` and the node for
  `(start vertex of -i, gate containing -i)`.
- For each fixed vertex, test whether all its gate nodes lie in the same
  connected component of this global graph.

Deliverable: explicit reasoned failure on bad case (expected).

## Validation against known bad case (Phase E)

For `n=6_5` bad path:

1. Confirm matrix remains primitive.
2. Confirm gate output structurally matches expected non-connectivity pattern.
3. Compare qualitatively with `n=6_5_bad_tt.train` gate report.
4. Record any differences in labeling only (not structure).

Deliverable: reproducible false-positive explanation in code output.

## Regression and integration (Phase F)

1. Add assertions in `test_issue2_bad_path.cpp` that this case fails the gate
   connectivity criterion.
2. If stable, extract reusable helpers (derivative/gates/connectivity) from the
   test into `include/traintracks/` / `include/ttauto/`.
3. Add one positive control case where matrix+gates both indicate pA.

Deliverable: issue-2 regression guard in CI/test workflow.

## Suggested implementation order

1. Phase A (vertex-directed alphabet) in reproducer.
2. Phase B main+infinitesimal derivative map.
3. Phase C gates per vertex with full partition printout.
4. Phase E validation on bad case before coding connectivity logic.
5. Phase D connectivity check.
6. Phase F regression assertions + optional helper extraction.

## Risks and mitigations

- **Sign/orientation mismatch**
  - Mitigation: print full `D` table and verify simple hand-checked relations.
- **Ambiguity in infinitesimal conventions**
  - Mitigation: keep one explicit canonical convention in comments and tests.
- **Confusing label translation vs `.train` output**
  - Mitigation: treat `.train` numbering as display-only; keep canonical
    internal IDs for all logic.
- **False agreement from too-short iterate bound**
  - Mitigation: bound by full directed-domain size + 1.

## Immediate next coding step

Add infinitesimal directed generators to the per-vertex outgoing sets in
`tests/test_issue2_bad_path.cpp`, then print those sets before gate
computation.

## Where to look next for theoretical understanding

Need to understand better the theory behind gates.  See B&H Prop. 3.1.3, and
Toby Hall's book pages 10--, with an emphasis on Lemma 3 (page 18) and
Corolary 4, which are key.  Also look at the `trains` code directly.

## Current implementation snapshot

What has been implemented already in the bad-path reproducer:

- Hardwired bad case reconstruction (`n=6`, `trk=4`, `sgidx=0`, fixed cycle and
  branch sequence).
- Transition-matrix checks (`matrix(path)==matrix(from_map)`) and primitive
  assertion.
- Derivative map printout for signed main generators.
- Directed-generator inventory by local vertex including infinitesimal labels.
- Gate partitioning by eventual derivative coincidence.
- Initial gate-connectivity check wired to infinitesimal joins.

What this currently gets right:

- Good diagnostic visibility: map, matrix, derivative, and local gate data all
  print in one place.
- Confirms matrix criterion alone is satisfied on the bad example.

What is still unresolved:

- Connectivity output does not match trusted `.train` locality pattern (where
  only a few vertices fail), indicating a model mismatch.
- Need to reconcile our local-vertex definition with `.train` switch vertices.
- Need to finalize whether infinitesimals are represented as gate elements,
  connectors, or both in the exact BH/trains convention we want to replicate.

## Revised immediate tasks

1. Build an explicit correspondence table between current local objects and the
   `.train` "Edges at vertex are ..." lines for the bad case.
2. Re-express gate computation on that matched vertex set.
3. Rebuild the infinitesimal join graph using that same matched locality.
4. Re-check that failures localize to the expected problematic vertices.
