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

- `include/traintracks/`: core train-track data model (traintrack/multigon/edge), map conversions, the planar embedding, braid words, helper math structures.
- `lib/traintracks/`: implementations of core transformation operations (construction, fold, relinking, normalization).
- `include/ttauto/`: automaton graph, path representation, search algorithm, result grouping, braid of a closed path.
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
- **Proper embedding**: the track drawn in the disc with the punctures on the real axis and no main edge crossing the segments below them; it is what gives the punctures an order.
- **Boundary walk**: the walk around the complementary region that touches the boundary of the disc; it meets each main edge twice and each side once.
- **Braid word**: a word in the generators `sigma_i`, held by `traintracks::braidword`, signed for inverses.
- **Full twist**: `delta^n` where `delta = sigma_1 ... sigma_{n-1}`; central, and the ambiguity in the braid of a folding path.
- **Dynnikov coordinates**: a coordinate system on measured loops whose growth under a braid gives its dilatation; used here only to check a braid against the path it came from.

## Typical End-to-End Execution Flow

Most programs follow this pipeline:

1. Build initial train tracks (strata representatives) via `build_traintrack_list(...)`.
2. Select one initial `traintrack`.
3. Construct `ttfoldgraph<traintrack>` from it.
4. Optionally decompose with `subgraphs(...)`.
5. Configure and run `ttauto<traintrack>::search(...)`.
6. Inspect results through `pA_list()` and print/export utilities.
7. Read the braid of an accepted closed path with `ttauto::folding_path_braid`.

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

### `include/traintracks/embedding.hpp` and `lib/traintracks/embedding.cpp`

The proper embedding of a track, as far as it is combinatorial.

- `outer_embedding(num, cut_dart = 0)`: walk the boundary region and
  return a `tt_embedding` holding the walk, `puncture_order`,
  `position_of`, the exterior cusps and `cut_dart`.  Pass `0` for the
  canonical cut, the loop of the root monogon.
- `all_directions_at_prong(num, q)`: like `ttnumbering::directions_at_prong`
  but including the sides of unpunctured multigons, which face walking
  needs and gates do not.
- `transported_cut_dart(before, fm, cut_dart)`: carry a cut across a fold.

Purpose: give the punctures an order along the real axis.

### `include/traintracks/collapsed_layout.hpp` and `lib/traintracks/collapsed_layout.cpp`

- `make_collapsed_layout(num, emb)`: the collapsed drawing, every
  multigon a point and every main edge one cubic arc, with the punctures
  on the axis in `emb`'s order and each multigon's prongs equally spaced.
- `cubic_point(c, t)`: a point of one of those arcs.

Purpose: what `examples/ttplot` draws, kept in the library so that its
planarity can be tested.

### `include/traintracks/braid.hpp` and `lib/traintracks/braid.cpp`

- `braidword`: a word in the braid generators, with `inverse()`,
  `reduce()`, `permutation()`, `exponent_sum()`, `block_swap()`,
  `delta()`, and `growth()`, the last from the action on Dynnikov
  coordinates and independent of everything else here.
- `fold_block_swap(before_num, fm, before)`: what one fold does to the
  punctures, read off the track before it, as a `fold_swap`.
- `fold_braid(...)`, `rotation_braid(n,k)`: the same as a braid word.

Purpose: express a fold, and a whole path, as a braid.

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

- `build_graph(...)`: builds the graph from one train track by trying all fold slots and discarding identity transitions, with `new_vertex(...)` creating each vertex.  It keeps an explicit worklist rather than recursing once per vertex, which used to exhaust the stack on large strata (issue #14), and indexes vertices by their coding rather than scanning the tracks stored so far.  Both the vertex numbering and the graph are exactly what the recursion produced.
- `find_symmetries()` / `sort_by_symmetries()`: computes reflection/cyclic symmetry relations and reorders vertices accordingly.  Reflection symmetry also goes through a coding index; comparing every pair of vertices cost more than building the graph.
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

There are two search modes.  Length-bounded is the default: set
`max_pathlength(...)` and leave `check_norms` off.  Norm-bounded needs
`max_dilatation(...)` and `check_norms()`, and is the mode that reproduces
the published minimum dilatations.

- `check_norms(...)`: enable the norm-bounded mode.  It recomputes the
  path-length bound from the dilatation window, so it discards any value
  set by `max_pathlength(...)`; so does a later `max_dilatation(...)`.
- `min_dilatation(...)`, `max_dilatation(...)`: acceptance window for candidate dilatation.
- `max_pathlength(...)`: path-length bound, 0 for none.  Set it after
  `check_norms()` and `max_dilatation()` if you want it to survive.
- `badword_length(...)`: prune paths that traverse a closed loop twice,
  in either mode; 0 (the default) disables it.  Valid only when hunting a
  minimum: repeating such a loop can only raise the dilatation, so the
  minimiser never contains one, but the repeated path is usually a
  pseudo-Anosov class of its own that the prune silently drops.  Measured
  in `devel/iss021/badwords.md`.
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

### `folding_path_braid` (`include/ttauto/path_braid.hpp`)

- `folding_path_braids(p)`: every braid a closed path can stand for; more
  than one only at a cyclically symmetric initial vertex.
- `folding_path_braid(p, &verified)`: the braid, with `verified` saying
  whether its growth matched the Perron root of the path's matrix.
- `detail::fold_index_of_branch(ttg, v, branch)`: recover the fold index a
  branch stands for, by matrix and target coding, since `ttfoldgraph` does
  not record it.

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
- `examples/ttbraid.cpp`: the braid of a closed folding path (issue #4).

## Test Files Worth Reading First

- `tests/test_traintrack.cpp`: core transformation invariants (attach/insert/fold/normalize/copy/check).
- `tests/test_folding_path.cpp`: path closure/equality/cyclic behavior/hash behavior.
- `tests/test_ttmap.cpp`: consistency between fold maps, transition matrices, and composition conventions.
- `tests/test_permplus1.cpp`: correctness of `mathmatrix_permplus1` representation and multiplication.
- `tests/test_badwords.cpp`: construction and reporting of badword filters.
- `testsuite/traintracks/test_embedding.cpp`: boundary-walk identities at every automaton vertex, and the puncture order of the two hand-drawn tracks of `doc/ttauto.tex`.
- `testsuite/ttauto/test_braid_extraction.cpp`: braids of closed paths checked against the Perron root of their own path, plus three pinned words.

## Practical "Where to Change What"

- Topological modifications and transformations, coding, cusp/fold mechanics: `include/traintracks/traintrack.hpp`, `lib/traintracks/traintrack.cpp`, `include/traintracks/multigon.hpp`, `lib/traintracks/multigon.cpp`, `include/traintracks/edge.hpp`.
- Fold map and matrix conventions: `include/traintracks/map.hpp`, `include/traintracks/fold_map.hpp`, `include/traintracks/map_labels.hpp`, and `tests/test_ttmap.cpp`.
- Gate test (pseudo-Anosov versus reducible): `include/traintracks/gates.hpp`, `lib/traintracks/gates.cpp`, `testsuite/traintracks/test_gates.cpp`.
- Automaton construction/symmetry/decomposition: `include/ttauto/ttfoldgraph.hpp`.
- DFS pruning and candidate acceptance logic: `include/ttauto/ttauto.hpp`, `include/ttauto/badwords.hpp`.
- Result grouping/serialization: `include/ttauto/pAclass.hpp`.
- Proper embedding, puncture positions, braid extraction: `include/traintracks/embedding.hpp`, `include/traintracks/braid.hpp`, `include/ttauto/path_braid.hpp`, `testsuite/traintracks/test_embedding.cpp`, `testsuite/ttauto/test_braid_extraction.cpp`.

## Braids from Folding Paths (issue #4)

A closed folding path defines a homeomorphism of the punctured disc, and
`ttauto::folding_path_braid` reads it off as a word in the braid
generators.  What makes that possible is that the multigons and main edges
form a tree, so the cyclic order of prongs recorded by the coding is a
complete rotation system and the track has one planar embedding up to
reflection.  `traintracks::outer_embedding` walks the boundary region of
that embedding and puts the punctures in the order they occupy along the
real axis.  A fold that lands on an unpunctured multigon keeps the track
properly embedded and moves nothing; a fold that lands on a punctured
monogon takes the moved edge once around the puncture, and undoing that
exchanges the punctures of the moved subtree with the target puncture,
which are next to each other in the walk.

Conventions that had to be calibrated rather than derived: the sense of the
boundary walk, fixed by the two tracks drawn by hand in `doc/ttauto.tex`
Fig. 2 and asserted in `testsuite/traintracks/test_embedding.cpp`; and the
handedness of each swap, which follows the fold direction.

Every braid is checked against the path it came from.  `braidword::growth`
iterates the braid on Dynnikov coordinates, which shares nothing with the
rest of the library, and for a pseudo-Anosov path that growth must equal
the Perron root of the path's transition matrix.  `folding_path_braid`
reports the outcome through its `verified` flag.

Two things the extraction has to get right that are easy to miss.  The
branch indices of the automaton are not fold indices: `build_graph`
numbers branches by their rank among the folds with a non-identity
transition matrix and throws the fold index away, and a subgraph
renumbers them again, so `fold_index_of_branch` identifies the fold by its
transition matrix and the coding of the track it produces.  And each step
must continue from the automaton's own copy of the track rather than from
the one just folded: the two have the same coding, but at a cyclically
symmetric vertex they need not be the same physical track, and it is the
automaton's copy the next branch index refers to.  Getting that wrong
gives a braid that is right up to a root of the full twist, which is
enough to change the dilatation.

The braid is only ever defined up to the full twist, since nothing here
records a framing at the boundary of the disc, so the final rotation takes
whichever way round is shorter and the words come out with exponent sum
zero.

Verified this way: 494 paths at three punctures to length 8, 1284 at four
to length 6, 1062 at five to length 5 and 2070 at six to length 5, every
vertex of every stratum, with no failures.
