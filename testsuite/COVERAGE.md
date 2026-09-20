# Testsuite Coverage Map

This file maps public headers to deterministic CTest programs in `testsuite/`.

## Measured line coverage (2026-09-19)

gcov on a Debug build of the ten fast testsuite programs (see "Measuring
line coverage" in `TESTING.md`; header lines are merged over all
translation units).  Before/after the coverage pass of 2026-09-19.

| File | Lines | Before % | After % | Left uncovered |
|---|---:|---:|---:|---|
| `lib/traintracks/traintrack.cpp` | 182 | 64 | 88 | string constructor (exits "Broken?"), error exits |
| `lib/traintracks/coding.cpp` | 249 | 79 | 83 | `ttnumbering::print`, error exits |
| `lib/traintracks/multigon.cpp` | 152 | 60 | 74 | `insert_edge(int,int)`, `print_details`, error exits |
| `lib/traintracks/build.cpp` | 213 | 69 | 69 | 3- and 4-multigon constructors, `build_traintrack_list_sweep_bigons` |
| `lib/traintracks/fold_map.cpp` | 183 | 89 | 89 | `fold_map_data::print`, fail-fast paths |
| `lib/traintracks/gates.cpp` | 186 | 92 | 92 | `bh_vertex_of`, sized accumulator constructor, fail-fast paths |
| `include/traintracks/traintrack.hpp` | 70 | 85 | 96 | one error branch |
| `include/traintracks/mathmatrix_permplus1.hpp` | 136 | 69 | 76 | the `exit(1)` branches of the dense constructor |
| `include/traintracks/map.hpp` | 26 | 73 | 73 | error branches of `transition_matrix_from_map` |
| `include/traintracks/map_labels.hpp` | 16 | 75 | 75 | |
| `include/traintracks/coding.hpp` | 16 | 100 | 100 | |
| `include/traintracks/edge.hpp` | 53 | 91 | 91 | |
| `include/traintracks/multigon.hpp` | 41 | 95 | 95 | |
| `include/traintracks/util.hpp` | 17 | 100 | 100 | |
| `include/ttauto/ttauto.hpp` | 329 | 78 | 89 | `debug` branches, badword pruning inside the search, symmetric-norm variant |
| `include/ttauto/folding_path.hpp` | 189 | 76 | 80 | error exits, `find_vertices` on a bad path |
| `include/ttauto/pAclass.hpp` | 91 | 97 | 90 | print cap branches |
| `include/ttauto/ttfoldgraph.hpp` | 183 | 95 | 95 | |
| `include/ttauto/badwords.hpp` | 61 | 100 | 100 | |
| `include/ttauto/path.hpp` | 21 | 90 | 90 | |
| **total** | 2416 | 80 | 86 | |

Most of what remains uncovered is deliberate: fail-fast `exit(1)` paths,
diagnostic printers, the broken string constructor, and builders the
automaton never uses.  The `pAclass` figure dropped because `paths()`
added lines that the census, not the testsuite, exercises.

## traintracks headers

- `include/traintracks/mathmatrix_permplus1.hpp`
  - Primary: `testsuite/traintracks/test_mathmatrix_permplus1.cpp`
  - Coverage focus:
    - every one-fold transition matrix of the n=4 and n=5 (first stratum)
      automata: dense round trip, agreement with the unit-weight oracle,
      +1 entry present exactly for legal folds, `order()` verified by dense
      powers, left/right products against dense matrices, printers
    - hand-built permutation, permutation+1, identity and cyclic cases

- `include/traintracks/traintrack.hpp`
  - Primary: `testsuite/traintracks/test_traintrack_core.cpp`
  - Secondary: `testsuite/traintracks/test_map_consistency.cpp`
  - Coverage focus:
    - every stratum for n=3..5 and every one-fold neighbour: `check()`,
      idempotent `normalise()` (coding and numbering unchanged), coding
      round trip, mirror image from the reversed coding, symmetry
      accessors, prong/puncture/monogon/cusp counts
    - `print_coding()` output parsed back to the coding vector (the
      intended behaviour of the string constructor, which is not exercised
      because it exits with "Broken?")
    - `set_label()`, `pure_braid()`, `print()`, `print_singularity_data()`,
      `printMathematicaForm()`
    - weights set through the iterator are read back in the same order

- `include/traintracks/build.hpp`
  - Primary: `testsuite/traintracks/test_traintrack_core.cpp`
  - Secondary: `testsuite/traintracks/test_map_consistency.cpp`
  - Coverage focus:
    - representative train-track list generation for deterministic fixtures

- `include/traintracks/coding.hpp` (`ttnumbering`)
  - Primary: `testsuite/traintracks/test_numbering.cpp`
  - Coverage focus:
    - prong and edge numbering is a bijection consistent with the track
    - edge order agrees with `weights()`
    - invariance under copy and re-normalisation
    - labels are a function of the coding (`same_labels`): every fold
      result of the n=4 s1, n=5 s1, n=6 s3 and s5 automata matches the
      track rebuilt from its coding and the stored target vertex; at
      cyclically symmetric vertices exactly `cyclic_symmetry().order()`
      start monogons give the same labels

- `include/traintracks/fold_map.hpp`, `include/traintracks/map.hpp`
  - Primary: `testsuite/traintracks/test_map_consistency.cpp`,
    `testsuite/traintracks/test_fold_map_paths.cpp`
  - Coverage focus:
    - one-step fold map/matrix consistency; side letter on the target multigon
    - hand-checked one-step and composed words for the n=3 example
    - every one-step and composed word is a continuous edge path in the
      canonical numbering, maps tails to tails and heads to heads, and
      abelianises to the transition matrix (all branches of the n=3, 4, 5
      graphs and of the 110-vertex n=6 stratum 3,3(2) graph, random paths,
      iterated closed paths)
    - the visited graphs contain vertices with nontrivial cyclic symmetry
      (n=4 stratum 1, n=6 stratum 5) and reflection symmetric vertices, so
      the automorphic-vertex identification is exercised, and the test
      fails if they disappear

- `include/traintracks/map_labels.hpp`
  - Primary: `testsuite/traintracks/test_map_consistency.cpp`
  - Coverage focus:
    - main/infinitesimal label indexing used in fold map checks

- `include/traintracks/gates.hpp`
  - Primary: `testsuite/traintracks/test_gates.cpp`
  - Coverage focus:
    - Bestvina-Handel gate test on the n=6 stratum 3,3(2) subautomaton: the
      known spurious pseudo-Anosov is rejected at exactly the fixed 3-edge
      monogon, four genuine pseudo-Anosovs pass
    - the minimal example on n=4 stratum (2): the length-4 path with
      polynomial (x-1)(x^2-3x+1) is rejected at its fixed 3-edge monogon
      with the same gate pattern; exhaustive counts show nothing fails at
      length <= 3 for n <= 4 and nothing at all for n=3 or n=4 stratum 3(1)
    - accumulator-based and word-based analyses agree
    - puncture corollary and the Props. 3.3.3-3.3.4 shape pattern on every
      connected primitive closed path of a small automaton; the sweep also
      reports how many accepted paths have a gate partition at an
      unpunctured multigon finer than its prongs (zero as of 2026-09-19)

- `include/traintracks/edge.hpp`
- `include/traintracks/multigon.hpp`
  - Covered indirectly via `traintrack` mutation checks:
    - `testsuite/traintracks/test_traintrack_core.cpp`
    - `testsuite/traintracks/test_map_consistency.cpp`

## ttauto headers

- `include/ttauto/ttfoldgraph.hpp`
  - Primary: `testsuite/traintracks/test_map_consistency.cpp`
  - Secondary: `testsuite/ttauto/test_folding_path_and_badwords.cpp`
  - Coverage focus:
    - automaton construction and fold-branch accessors

- `include/ttauto/folding_path.hpp`
- `include/ttauto/path.hpp`
  - Primary: `testsuite/ttauto/test_folding_path_and_badwords.cpp`
  - Secondary: `testsuite/traintracks/test_map_consistency.cpp`,
    `testsuite/traintracks/test_fold_map_paths.cpp`
  - Coverage focus:
    - `subpath()` from both ends, `operator*`/`operator*=`,
      `ending_equals()`, `cycle_path()`, `initial_vertex()` on an empty
      path, `clear()`, matrix of a composite = product of the pieces
    - closed-path cyclic equality and hash behaviour; open paths compare
      directly

- `include/ttauto/badwords.hpp`
  - Primary: `testsuite/ttauto/test_folding_path_and_badwords.cpp`
  - Coverage focus:
    - content of `badwords(ttg, L)`: every entry at (v, l) is the square
      h*h of a closed path h of length l+1 from v whose transition matrix
      has the same zero pattern as its square; the length-1 layer is
      complete; recomputation gives the same table

- `include/ttauto/ttauto.hpp`
  - Primary: `testsuite/ttauto/test_ttauto_search.cpp`
  - Coverage focus:
    - n=5 first stratum, length 8: the exact list of eight classes with
      dilatations, characteristic polynomials, shortest and longest
      representatives; 3555 gate candidates, none rejected
    - every stored representative: matrix primitive, class polynomial equals
      an independent characteristic polynomial (principal minors, Bareiss),
      dilatation equals the power-iteration spectral radius, gates connected
    - output helpers

- `include/ttauto/ttauto.hpp` (norm-bounded mode)
  - Primary: `testsuite/ttauto/test_min_dilatations.cpp`
  - Coverage focus:
    - `check_norms()` + `max_dilatation()` as in `examples/ttauto_min_example`
      for n=3, 4, 5: the literature minima (2.61803; 2.61803, 2.29663;
      1.72208, 1.72208, 2.15372, 2.01536) with their polynomials and
      shortest lengths
    - the n=4 first stratum also produces a reducible candidate with
      polynomial (x-1)(x^2-3x+1) and the same dilatation, which the gate
      test moves to `rejected_pA_list()`; with `check_gates(false)` it is
      reported as a pA
    - `check_all_norms()` and `find_maxnorm()` are exercised only here

- `include/ttauto/ttauto.hpp` (gate test in the search)
  - Primary: `testsuite/ttauto/test_ttauto_gates.cpp`
  - Coverage focus:
    - `check_gates` on/off on the n=6 stratum 3,3(2) subautomaton: the
      2.01536 class moves to `rejected_pA_list()`, all others unchanged

- `include/ttauto/pAclass.hpp`
  - Primary: `testsuite/ttauto/test_ttauto_search.cpp`
  - Coverage focus:
    - per-class invariants (`number_of_paths`, `shortest`, `longest`, dilatation)
    - class-level data returned through `ttauto::pA_list()`

## Oracles

`testsuite/oracles.hpp` holds slow, independent implementations used to
cross-check the library and its dependencies:

- `unit_weight_transition_matrix(tt, f)`: the original n-fold way of
  reading a one-fold transition matrix (against the fold record).
- `charpoly_by_minors(M)`: det(xI - M) from sums of principal minors with a
  fraction-free Bareiss determinant (against `jlt::mathmatrix::charpoly()`,
  which is the trace-based recursion and returns det(M - xI)).
- `power_iteration_radius(M)`: spectral radius of a primitive matrix
  (against `ttauto::findroot()` on the characteristic polynomial).

## Notes

- CTest executes only tests defined from `testsuite/`.
- Testsuite targets are compiled with `-UNDEBUG`; use `CHECK()` from
  `testsuite/check.hpp` in new tests.
- The testsuite also builds and passes in a Debug (`-O0`) configuration;
  a coverage build exposed a link error there (`edge::nends` was a
  non-inline `static const`, fixed 2026-09-19).
- Slow scan-strata markdown regression is provided by
  `testsuite/ttauto/test_scan_strata_markdown.sh` and is enabled with
  `-DTTAUTO_ENABLE_SLOW_TESTS=ON`.
- Existing `tests/` programs remain available for exploratory or feature-branch work.
