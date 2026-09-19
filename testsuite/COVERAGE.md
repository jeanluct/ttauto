# Testsuite Coverage Map

This file maps public headers to deterministic CTest programs in `testsuite/`.

## traintracks headers

- `include/traintracks/mathmatrix_permplus1.hpp`
  - Primary: `testsuite/traintracks/test_mathmatrix_permplus1.cpp`
  - Coverage focus:
    - permutation/permutation+1 construction
    - dense round-trip via `full()`
    - left/right multiplication agreement with dense matrices

- `include/traintracks/traintrack.hpp`
  - Primary: `testsuite/traintracks/test_traintrack_core.cpp`
  - Secondary: `testsuite/traintracks/test_map_consistency.cpp`
  - Coverage focus:
    - construction, normalization, coding round-trip
    - fold mutation invariants (`check()` before/after)
    - weight traversal consistency
    - fold cusp and infinitesimal index resolution in-range

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

- `include/traintracks/fold_map.hpp`, `include/traintracks/map.hpp`
  - Primary: `testsuite/traintracks/test_map_consistency.cpp`,
    `testsuite/traintracks/test_fold_map_paths.cpp`
  - Coverage focus:
    - one-step fold map/matrix consistency; side letter on the target multigon
    - hand-checked one-step and composed words for the n=3 example
    - every one-step and composed word is a continuous edge path in the
      canonical numbering, maps tails to tails and heads to heads, and
      abelianises to the transition matrix (all branches of the n=3, 4, 5
      graphs, random paths, iterated closed paths)

- `include/traintracks/map_labels.hpp`
  - Primary: `testsuite/traintracks/test_map_consistency.cpp`
  - Coverage focus:
    - main/infinitesimal label indexing used in fold map checks

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
  - Secondary: `testsuite/traintracks/test_map_consistency.cpp`
  - Coverage focus:
    - path composition and map/matrix accumulation
    - closed-path cyclic equality and hash behavior

- `include/ttauto/badwords.hpp`
  - Primary: `testsuite/ttauto/test_folding_path_and_badwords.cpp`
  - Coverage focus:
    - deterministic badword matrix dimensions and construction on fixed fixture

- `include/ttauto/ttauto.hpp`
  - Primary: `testsuite/ttauto/test_ttauto_search.cpp`
  - Coverage focus:
    - deterministic bounded search smoke on a small fixed automaton
    - structural invariants of discovered pseudo-Anosov classes
    - output helper calls remain valid on non-empty results

- `include/ttauto/pAclass.hpp`
  - Primary: `testsuite/ttauto/test_ttauto_search.cpp`
  - Coverage focus:
    - per-class invariants (`number_of_paths`, `shortest`, `longest`, dilatation)
    - class-level data returned through `ttauto::pA_list()`

## Notes

- CTest executes only tests defined from `testsuite/`.
- Testsuite targets are compiled with `-UNDEBUG`; use `CHECK()` from
  `testsuite/check.hpp` in new tests.
- Slow scan-strata markdown regression is provided by
  `testsuite/ttauto/test_scan_strata_markdown.sh` and is enabled with
  `-DTTAUTO_ENABLE_SLOW_TESTS=ON`.
- Existing `tests/` programs remain available for exploratory or feature-branch work.
