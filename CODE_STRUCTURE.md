# ttauto Code Structure (Main Repository)

This document describes the structure of the main `ttauto` C++ codebase and explains what the major classes and functions are for.

Scope: this covers code in `include/`, `lib/`, `examples/`, `tests/`, and `testsuite/`. It does not document internals of `extern/jlt`.

## What This Project is About

For the mathematics behind the pipeline (folds as free-group substitutions,
the transition matrix, and the Bestvina-Handel gate condition that decides
whether a closed path is pseudo-Anosov) see `doc/ttauto.tex`.

At a high level, the project explores **train-track automata** associated with mappings on punctured discs and searches for pseudo-Anosov candidates.

In concrete terms, it does three main things:

1. Builds train-track objects representing a combinatorial/topological state.
2. Builds a directed graph whose edges correspond to valid fold operations.
3. Traverses that graph to find closed paths with matrix/map properties consistent with pseudo-Anosov behavior.

## High-Level Layout

- `include/traintracks/`: core train-track data model (traintrack/multigon/edge), map conversions, helper math structures.
- `lib/traintracks/`: implementations of core transformation operations (construction, fold, relinking, normalization).
- `include/ttauto/`: automaton graph, path representation, search algorithm, result grouping.
- `examples/`: runnable programs showing interactive and scripted usage.
- `tests/`: executable tests and consistency checks.
- `testsuite/`: deterministic CTest targets (with optional slow integration checks).

## Terms Used in This Codebase

- **Train track**: combinatorial object built from pronged singularities and branches (edges), used to encode dynamical invariance.
- **Multigon**: a single pronged singularity piece in a train track; has a number of prongs and edge slots at each prong.
- **Prong**: local vertex on a multigon boundary where one or more edge slots attach.
- **Cusp**: local foldable location between adjacent edge slots on a prong.
- **Fold**: local transformation that folds one edge onto another at a cusp and changes the train track.
- **Transition matrix**: matrix counting how branches map under a fold or sequence of folds.
- **Train-track map**: free-group automorphism representation of branch images under folds.
- **Folding automaton**: directed graph of train tracks connected by valid folds.
- **Folding path**: a sequence of folds in a train track graph.
- **Closed path**: folding path whose final vertex equals its initial vertex.

## Typical End-to-End Execution Flow

Most programs follow this pipeline:

1. Build initial train tracks (strata representatives) via `build_traintrack_list(...)`.
2. Select one initial `traintrack`.
3. Construct `ttfoldgraph<traintrack>` from it.
4. Optionally decompose with `subgraphs(...)`.
5. Configure and run `ttauto<traintrack>::search(...)`.
6. Inspect results through `pA_list()` and print/export utilities.

## Worked Example: One Fold Through the Stack

This section gives a concrete mental model of how one fold is represented across geometry, graph structure, and algebraic data.

### Step A: Start at a Train-Track Vertex

- You have one normalized `traintrack` object `tt` (for example, `ttg.traintrack(v)` inside the automaton).
- At this point, cusp indices are meaningful because normalization fixes the traversal convention used by cusp-ordering helpers.

### Step B: Choose a Fold Index `f`

- A fold index identifies one directed fold at one cusp.
- Internally, `fold_cusp_location(f, mmc, pc, ec)` resolves `f` into:
  - `mmc`: the multigon containing the cusp,
  - `pc`: prong index on that multigon,
  - `ec`: edge-slot index of the cusp start.
- `fold(f)` then applies the local mutation and re-normalizes the track.

### Step C: Build Algebraic Descriptions of the Same Fold

For the same `tt` and `f`, the code builds two algebraic objects:

- `fold_transition_matrix(tt, f)`:
  - gives the main-edge transition matrix for that one fold, read off the fold record (`fold_map_data::transition_matrix`: the edge permutation plus the extra entry for the moved edge running over `onto`),
  - represented compactly as `mathmatrix_permplus1`; the original unit-weight construction (n folds) survives only as `oracles::unit_weight_transition_matrix` in `testsuite/oracles.hpp`.
- `fold_traintrack_map(tt, f)` (via `traintrack::fold_with_map`, `include/traintracks/fold_map.hpp`):
  - gives the train-track map on generators (main edges + sides), in the canonical numbering (`ttnumbering`, `include/traintracks/coding.hpp`) of the track before and after the fold,
  - the moved edge's image is the three-letter path `moved copy . side . onto edge` (or reversed), where the side is the side of the *target* multigon traversed between the two target prongs; every other edge maps to its signed new number and every side to its new prong number.

Consistency rule: main-edge counts extracted from the map must agree with the transition matrix (`transition_matrix_from_map`; checked against the unit-weight oracle in `testsuite/traintracks/test_map_consistency.cpp` and `test_fold_map_paths.cpp`).

### Step D: Store as One Automaton Branch

In `ttfoldgraph`:

- source vertex = original track state,
- branch label = fold index `f`,
- target vertex = folded/normalized track,
- branch data includes:
  - target vertex id,
  - one-step transition matrix,
  - one-step train-track map.

So each graph edge is not just connectivity; it carries all one-step algebraic data needed for path composition.

### Step E: Compose Along a Folding Path

`folding_path` stores a sequence of branch choices. From that sequence it computes:

- path transition matrix (product of one-step matrices),
- path train-track map (composition of one-step maps),
- induced vertex path (start/end vertices, closure test).

This is why `folding_path` is the core DFS state in `ttauto`: it is both combinatorial (which edges were taken) and algebraic (what map/matrix they compose to).

### Step F: Search Accept/Reject in `ttauto`

During DFS in `ttauto`:

- pruning checks reject many partial paths early (norm bounds, badwords, depth limits, etc.),
- closed paths with a primitive matrix (irreducible and aperiodic) and dilatation above 1 become candidates,
- `record_pA` applies the dilatation window and then the Bestvina-Handel gate test (`folding_path::gates()`, `traintracks/gates.hpp`); candidates with a disconnected gate graph are reducible and go to `rejected_pA_list()`,
- accepted candidates are grouped into `pAclass` objects keyed by polynomial/dilatation.

In short: one fold becomes one graph edge with matrix+map payload; many edges compose into one candidate dynamical class.

## Common Reading Paths (How to Onboard Quickly)

If you are new to this code, these reading orders are effective.

### Path 1: "I want the big picture first"

1. `examples/ttauto_min_example.cpp`
2. `include/ttauto/ttfoldgraph.hpp`
3. `include/ttauto/folding_path.hpp`
4. `include/ttauto/ttauto.hpp`
5. `tests/test_ttmap.cpp`

### Path 2: "I need to modify fold mechanics"

1. `include/traintracks/traintrack.hpp`
2. `lib/traintracks/traintrack.cpp` (`fold`, cusp resolution helpers)
3. `include/traintracks/map.hpp`
4. `tests/test_traintrack.cpp`
5. `tests/test_ttmap.cpp`

### Path 3: "I need to change search behavior"

1. `include/ttauto/ttauto.hpp`
2. `include/ttauto/badwords.hpp`
3. `include/ttauto/folding_path.hpp`
4. `tests/test_badwords.cpp`
5. `examples/ttauto.cpp`

## Train Track Classes (`traintracks` namespace)

The foundation is three interacting classes:

- `traintracks::traintrack` in `include/traintracks/traintrack.hpp`.
- `traintracks::multigon` in `include/traintracks/multigon.hpp`.
- `traintracks::edge` in `include/traintracks/edge.hpp`.

### Data Relationships and Invariants

- `traintrack` owns the collection of multigons (`mgv`).
- each `multigon` stores the incident edges attached at each `(prong, edge-slot)` position.
- each `edge` stores endpoint metadata: which multigons it touches, at which prongs/slots.

Purpose of this design: fold/swap/relabel operations need to update many cross-links quickly while preserving consistency. The invariant is that every edge endpoint record and every multigon edge slot agree with each other.

### `traintrack`: Main Transformable Topological Object

`traintrack` represents one complete train track and is the core object used by graph construction and search.

Key public responsibilities:

- Construct a track:
  - `traintrack(N)`, `traintrack(N,K...)`, `traintrack(N,Kv)`: build standard strata shapes.
  - `traintrack(const intVec& code)`: rebuild from coding vector.
- Query structure and combinatorics:
  - `edges()`, `multigons()`, `monogons()`, `punctures()`, `cusps()`, `total_prongs()`.
- Canonicalize and compare:
  - `normalise()`: put track into canonical form (important before comparisons and fold indexing assumptions).
  - `coding(int dir=1)`: serialize normalized structure to comparable coding.
  - `operator==`: isotopy-style equality using normalized coding/multigon compatibility.
  - `cyclic_symmetry()`: detect rotational symmetry and return corresponding permutation structure.
- Perform folds:
  - `fold(int f)`: apply fold by global fold index in current cusp ordering.
  - `fold_cusp_location(...)`: map fold index to concrete `(multigon, prong, cusp-edge-slot)`.
  - `fold_with_map(f, fm)`: fold and fill a `fold_map_data` record (moved/onto edges, target prongs, side letter, edge and prong images, numberings before and after).
  - `numbering()`: canonical prong/edge numbering with orientations (`ttnumbering`), rooted at monogon 0 like `weights()` and the coding.
- Integrate with matrix/map layer:
  - `fold_transition_matrix(int f)`.
  - `fold_traintrack_map(int f)`.
- Manage weights and diagnostics:
  - `weights()` getter/setter traversal in coding order.
  - `check()` for structural consistency.
  - `print()`, `print_singularity_data()`, `print_coding()` for diagnostics/export.

Important internal helpers (private):

- `recursive_build`: reconstruct topology from coding blocks.
- `weights()`, `fold(f)` and `fold_cusp_location(f)` no longer walk the
  track themselves: they read the canonical numbering
  (`coding_engine::numbering`, which records edge order, prong numbers and
  the cusp order used by fold indices in one depth-first walk).  Only two
  walks remain: the coding (normal form, both directions) and the
  numbering (index).

Coding implementation note:

- Canonical coding logic now lives in `include/traintracks/coding.hpp` and
  `lib/traintracks/coding.cpp` (`traintracks::detail::coding_engine`).
- The same module provides `ttnumbering` and `coding_engine::numbering`, the
  canonical numbering of prongs and edges (with orientations) produced by
  the same depth-first walk as the coding and `weights()`.  Side `q` runs
  from prong `q` to the next prong of its multigon.  "Number" is used, not
  "label", because `label` is the puncture label of a multigon.
- `include/traintracks/fold_map.hpp` / `lib/traintracks/fold_map.cpp`:
  `fold_map_data` and `traintrack::fold_with_map`, the record of one fold in
  canonical numbering from which the train-track map is built.  Prongs are
  followed through `normalise()` by the set of edge objects attached to
  them (edge objects persist; multigon objects and edge ending indices do
  not).
- `traintrack::{coding, normalise, cyclic_symmetry, print_coding}` delegate to
  that coding module.

### `multigon`: Local Piece and Edge-Slot Manager

`multigon` is a local component representing one polygonal piece and the list of edge slots around each prong.

Key responsibilities:

- Attach and modify local edge slots:
  - `attach_edge(...)`: place an edge at a given slot.
  - `insert_edge(...)`: insert into a slot and shift later slots.
  - private `erase_edge_pointer(...)`: remove slot and renumber metadata.
  - private `point_to_edge(...)`: assign pointer in slot with checks.
- Traverse local cyclic order:
  - `cycle_edges(...)`: walk around a multigon in clockwise/anticlockwise order.
  - `cycle_prongs(...)`: rotate prong indexing.
  - `edge_sequence(...)`: summarize edges-per-prong from an offset.
- Local canonicalization/comparison:
  - `normalise()`: rotate prongs to maximal edge-sequence form.
  - `operator==`, `operator<`: compare local structure for sorting/equality workflows.
- Invariant checking:
  - `check()` validates edge slot and endpoint metadata consistency.

Special note on swap behavior:

- `swap(multigon&, multigon&)` in `lib/traintracks/multigon.cpp` is critical because it swaps contents and then repairs all edge endpoint metadata so that both swapped multigons remain consistent.

### `edge`: Endpoint Metadata and Weight Carrier

`edge` is intentionally low-level. It stores:

- a weight (`wt`), and
- two endpoint records (which multigon, which prong, which edge-slot index at that endpoint).

Key methods and purpose:

- `target_multigon(...)`: from one endpoint, find the opposite endpoint and metadata.
- `detach_from_multigon(...)`: detach one endpoint and request slot cleanup from the owning multigon side.
- `renumber_ending(...)`: update stored slot index when multigon slots shift.
- `relink_ending(...)`: move endpoint attachment to another multigon/slot.
- `check()`: verify endpoint-level consistency assumptions.

## Map and Matrix Utilities (`traintracks`)

### `include/traintracks/map.hpp`

This file is the bridge between geometric folds and algebraic representations.

Core functions:

- `fold_transition_matrix(const TrTr&, int f)`: computes one-fold main-edge transition matrix.
- `fold_traintrack_map(const TrTr&, int f)`: one-fold train-track map including side generators, built from `fold_with_map` on a copy (see `fold_map.hpp`).
- `transition_matrix_from_map(const TrTr&, const jlt::freeauto<int>&)`: projects map back to main-edge transition matrix.

Why this matters: the automaton stores both matrix and map data per branch; these functions keep conventions synchronized.

### `include/traintracks/gates.hpp` and `lib/traintracks/gates.cpp`

The Bestvina-Handel gate test for a train-track map (issue #2; the
mathematics is in `devel/iss002/issue2_gates.tex`).

- `fold_derivative(AM, target_numbering)`: the derivative `D` (first letter
  of each image) and the turns taken by the image words of the main edges,
  read in the target track's numbering.  Works for one-step and composed
  maps.
- `gate_accumulator`: `push_back` one `fold_derivative` per branch of a
  path; composes `D` and the realised turns (`T <- T_i U D_i(T)`) without
  ever forming the composed words.  `analyse(N)` on a closed path closes the
  turns under `D x D`, partitions the directions at each Bestvina-Handel
  vertex (an unpunctured multigon, or one prong of a punctured multigon)
  into gates by eventual coincidence under `D`, joins gates by realised
  turns, and reports connectivity and the allowed infinitesimal-edge shape
  per vertex (`gate_analysis`, `gate_vertex_report`).
- `analyse_gates(N, AM)`: word-based version of the same test, for tests and
  diagnostics.
- Fail-fast on impossible data: a realised turn inside a gate, a turn
  joining two vertices, a derivative that does not map directions to
  directions.

### `include/traintracks/map_labels.hpp`

`ttmap_labeler` defines edge (free group generator) indexing conventions used everywhere in map code.

- main edges are `1..nmain`.
- infinitesimal edges are `nmain+1..nmain+ninf`.
- sign encodes orientation (negative means inverse orientation).

Purpose: avoid ad-hoc index handling and keep map/matrix conversions consistent.

### `include/traintracks/mathmatrix_permplus1.hpp`

`mathmatrix_permplus1` is a compact matrix representation for matrices that are either:

- pure permutation matrices, or
- permutation matrices with one extra `+1` entry.

Key API and purpose:

- constructor from dense matrix: validate shape constraints and compress structure.
- `full()`: re-expand to dense matrix for generic code paths/tests.
- `is_perm()`, `is_identity()`: quick structural predicates.
- `row_perm()`, `column_perm()`, `plus1_row()`, `plus1_col()`: access compressed structure.
- `order()`: order of the permutation part.
- `operator*` overloads: efficient multiplication with dense matrices.

### `include/traintracks/build.hpp` and `lib/traintracks/build.cpp`

These files generate initial train-track representatives (strata seeds).

Public functions:

- `build_traintrack_list(int N, int N2 = 0)`: enumerate tracks for fixed puncture count (and optional punctured-bigon count).
- `build_traintrack_list_sweep_bigons(int N)`: aggregate over admissible bigon counts.

Purpose: provide standard start states for automaton construction.

## Automaton and Search Layer (`ttauto` namespace)

This layer is template-based and generally instantiated as `TrTr = traintracks::traintrack`.

### `ttfoldgraph<TrTr>` (`include/ttauto/ttfoldgraph.hpp`)

`ttfoldgraph` is the directed graph of train tracks under folds.

Stored per vertex/branch:

- vertex train-track object,
- outgoing target vertices,
- outgoing transition matrices,
- outgoing train-track maps,
- outgoing fold count.

Key behavior:

- `add_vertex(...)`: recursively builds the graph by trying all fold slots and discarding identity transitions.
- `find_symmetries()` / `sort_by_symmetries()`: computes reflection/cyclic symmetry relations and reorders vertices accordingly.
- branch accessors: `foldings(v)`, `target_vertex(v, br)`, `transition_matrix(v, br)`, `traintrack_map(v, br)`.
- decomposition utility: `subgraphs(...)` extracts invariant subgraphs via sparse matrix decomposition.

Purpose: provide the search space for pseudo-Anosov candidate detection.

### `folding_path<TrTr>` (`include/ttauto/folding_path.hpp`)

`folding_path` represents a sequence of fold choices and its induced vertex sequence in one fixed `ttfoldgraph`.

Key methods:

- edit path: `push_back`, `pop_back`, `clear`, `cycle_path`.
- query shape: `length`, `closed`, `initial_vertex`, `final_vertex`, `number_of_foldings`.
- derive algebra:
  - `transition_matrix()` for full path,
  - `traintrack_map()` for full path.
- substructure checks: `subpath(int)`, `ending_equals(...)`.
- equality/hash: closed-path equality supports cyclic reindexing of the same loop.

Purpose: this is the central traversal state object used by the search engine.

### `ttauto<TrTr>` (`include/ttauto/ttauto.hpp`)

`ttauto` performs depth-first search (DFS) on folding paths and records classes of accepted closed paths.

Public configuration knobs:

- `check_norms(...)`: enable matrix-bound pruning mode.
- `min_dilatation(...)`, `max_dilatation(...)`: acceptance window for candidate dilatation.
- `max_pathlength(...)`: hard path-length bound.
- `badword_length(...)`: configure repeated-pattern pruning.
- `max_paths_to_save(...)`, `max_paths_to_print(...)`, `print_path_every(...)`.
- `output_file(...)`: optional Mathematica-form output destination.

Public execution/results:

- `search(int tt00 = 0)`: run search starting from vertex ordering offset.
- `pA_list()`: retrieve grouped results.
- `print_pA_list()`, `print_pA_list_MathematicaForm()`: reporting helpers.

Important internal flow:

- `find_pAs()`: initialize per-start-vertex DFS state and counters.
- `descend_graph()`: one DFS step; applies pruning, checks closure, and handles backtracking.
- `check_all_norms()`: matrix-based lower-bound pruning checks.
- `record_pA()`: apply the dilatation window and the gate test, then insert/update the result class keyed by characteristic polynomial (`add_current_path`).
- `check_gates(bool)`: enable/disable the gate test (default on); `rejected_pA_list()`, `gate_candidates()`, `gate_rejected()` and `gate_rejection_rate()` expose the rejections (cumulative over the search); the statistics block prints "Gate test = R rejected of C candidates (P%)".

Purpose: this is the main "engine" of the repository.

### `pAclass<TrTr>` (`include/ttauto/pAclass.hpp`)

`pAclass` groups candidate results that share one characteristic polynomial/dilatation and stores representative paths.

Key methods:

- `add_path(...)`: add a closed path and transition matrix, honoring path-cap policy.
- `number_of_paths()`, `shortest()`, `longest()`: summary stats.
- `print_paths(...)`: human-readable path list.
- `print_pA_MathematicaForm(...)`: machine-friendly export form.

Purpose: deduplicate and summarize search output.

### `path` (`include/ttauto/path.hpp`)

Small bounded integer-sequence type used by `folding_path` for fold and vertex sequences.

Purpose:

- enforce range-aware comparison semantics,
- support enumeration (`operator++`) of fixed-length sequences.

### `badwords` (`include/ttauto/badwords.hpp`)

Builds pattern-based "bad word" filters to prune DFS branches that correspond to repeated/redundant patterns.

Key API:

- `badwords(const ttfoldgraph<TrTr>&, int maxplen)`.
- helper `pattern_equal(...)` for zero/nonzero pattern comparisons.

Purpose: reduce search cost without changing the core fold graph.

## Entry Points in `examples/`

- `examples/ttauto.cpp`: interactive CLI driver (stratum selection, optional subgraph splitting, search, optional file export).
- `examples/ttauto_min_example.cpp`: compact scripted driver for low-dilatation checks on small puncture counts.
- `examples/ttauto_torus.cpp`, `examples/ttauto_count.cpp`, `examples/ttauto_labels.cpp`: additional scenarios/sweeps.

## Test Files Worth Reading First

- `tests/test_traintrack.cpp`: core transformation invariants (attach/insert/fold/normalize/copy/check).
- `tests/test_folding_path.cpp`: path closure/equality/cyclic behavior/hash behavior.
- `tests/test_ttmap.cpp`: consistency between fold maps, transition matrices, and composition conventions.
- `tests/test_permplus1.cpp`: correctness of `mathmatrix_permplus1` representation and multiplication.
- `tests/test_badwords.cpp`: construction and reporting of badword filters.

## Practical "Where to Change What"

- Topological modifications and transformations, coding, cusp/fold mechanics: `include/traintracks/traintrack.hpp`, `lib/traintracks/traintrack.cpp`, `include/traintracks/multigon.hpp`, `lib/traintracks/multigon.cpp`, `include/traintracks/edge.hpp`.
- Fold map and matrix conventions: `include/traintracks/map.hpp`, `include/traintracks/fold_map.hpp`, `include/traintracks/map_labels.hpp`, and `tests/test_ttmap.cpp`.
- Gate test (pseudo-Anosov versus reducible): `include/traintracks/gates.hpp`, `lib/traintracks/gates.cpp`, `testsuite/traintracks/test_gates.cpp`.
- Automaton construction/symmetry/decomposition: `include/ttauto/ttfoldgraph.hpp`.
- DFS pruning and candidate acceptance logic: `include/ttauto/ttauto.hpp`, `include/ttauto/badwords.hpp`.
- Result grouping/serialization: `include/ttauto/pAclass.hpp`.
