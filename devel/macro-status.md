# Macro Status (Non-Guard Macros)

This file tracks build/runtime feature macros used in the main `ttauto` codebase
(excluding include guard macros).

## Core / traintracks-layer macros

- `TRAINTRACKS_NO_SHARED_PTR`
  - Default: **not defined**
  - Enabled when: compile with `-DTRAINTRACKS_NO_SHARED_PTR`
  - Effect: switches core pointer aliases from shared pointers to raw pointers in
    compatibility paths.
  - Status: **supported compatibility mode** (kept intentionally for portability).

`TRAINTRACKS_OLD_HASH` has been removed.

## Automaton / ttauto-layer macros

- `TTAUTO_HASH_BADWORDS`
  - Default: **not defined**
  - Note: currently forced off by `#undef TTAUTO_HASH_BADWORDS` in
    `include/ttauto/badwords.hpp`.
  - Effect (when enabled): use hashed container path for badword storage.

- `TTAUTO_USE_FORTRAN`
  - Default: **not defined**
  - Effect: enables optional Fortran-backed norm checks in `ttauto` search logic.

- `TTAUTO_CHECK_SYMMETRIC_NORM`
  - Default: **not defined**
  - Effect: enables extra symmetric-norm consistency checks in `ttauto` pruning.

## Example-local macro

- `TTAUTO_RELEASE`
  - Scope: local to `examples/ttauto.cpp`
  - Effect: toggles release warning/banner text block.

## How to test non-default macro builds

Yes - this is easy to test with CMake by configuring a dedicated build directory
with `CMAKE_CXX_FLAGS` defines.

Examples:

```bash
# raw-pointer compatibility path
cmake -S . -B build-macro-nosptr -DCMAKE_CXX_FLAGS='-DTRAINTRACKS_NO_SHARED_PTR'
cmake --build build-macro-nosptr --target test_test_traintrack -j

# optional symmetric norm checks in automaton layer
cmake -S . -B build-macro-symnorm -DCMAKE_CXX_FLAGS='-DTTAUTO_CHECK_SYMMETRIC_NORM'
cmake --build build-macro-symnorm --target test_test_badwords -j

# optional Fortran-enabled path (compile-time check)
cmake -S . -B build-macro-fortran -DCMAKE_CXX_FLAGS='-DTTAUTO_USE_FORTRAN'
cmake --build build-macro-fortran --target test_test_badwords -j
```

Notes:

- Prefer dedicated build directories per macro mode to avoid stale-object reuse
  across incompatible macro settings.
- Build one or a few representative targets per macro mode for quick checks,
  then scale to full `cmake --build <dir> -j` if needed.

## Fixed preprocessor toggles (`#if 0` / `#if 1`) and removal rationale

These are not named macros, but they are compile-time toggles worth tracking.

- `include/ttauto/ttauto.hpp` — warning block near end of `search()` (`#if 0`)
  - Current behavior: warning about badword skipping is compiled out.
  - Removal guidance: **probably safe to delete**.
  - Rationale: this is user-facing logging only; deleting it does not change
    search/pruning semantics, only whether a warning can be printed.

- `include/ttauto/ttauto.hpp` — `pseudoAnosov` stats line (`#if 0`)
  - Current behavior: line is disabled and already labeled redundant.
  - Removal guidance: **safe to delete**.
  - Rationale: dead output-only code with explicit in-code note that it is not
    meaningful at that stage.
  - Additional comment by J-LT: this is a tricky issue given the fact that we
    now know that there are cases where the matrix is irreducible, but the map
    is not pseudo-Anosov.  Need further checks (issue #2).

- `include/ttauto/ttfoldgraph.hpp` — symmetry-sort strategy selector (`#if 1/#else`)
  - Current behavior: uses contiguous self/left/right block order; alternate
    pairwise ordering retained in `#else`.
  - Removal guidance: **keep unless intentionally retiring alternate ordering**.
  - Rationale: both branches encode valid policies; this is design/analysis
    behavior, not dead debug scaffolding.

- `include/traintracks/traintrack.hpp` + `lib/traintrack/build.cpp` —
  disabled `traintrack(const intVec& Kv)` constructor/implementation (`#if 0`)
  - Current behavior: constructor is unavailable.
  - Removal guidance: **higher-risk, keep for now**.
  - Rationale: comments document historical ambiguity with coding-based
    construction. Deleting loses that provenance and easy reactivation point.

- `include/traintracks/traintrack.hpp` — convenience fold overload (`#if 0`)
  - Current behavior: one-argument `fold_transition_matrix` convenience wrapper
    is disabled; two-output API is used instead.
  - Removal guidance: **higher-risk, keep for now**.
  - Rationale: small API-shape choice; currently disabled but potentially useful
    for call-site simplification experiments.

- `include/traintracks/edge.hpp` — `attach_to_multigon` declaration/implementation (`#if 0`)
  - Current behavior: API is disabled with explanatory shared_ptr caveats.
  - Removal guidance: **higher-risk, keep for now**.
  - Rationale: block documents an ownership pitfall and why naive reintroduction
    would be unsafe (`shared_from_this` concerns).

- `lib/traintrack/build.cpp` — strata debug print block (`#if 0`)
  - Current behavior: optional debug dump is disabled.
  - Removal guidance: **safe to delete**.
  - Rationale: pure diagnostics with no behavioral effect.

- `tests/test_badwords.cpp` — fixture selector (`#if 1/#else`)
  - Current behavior: default fixture (`n=6`, `traintrack(n,3)`, `maxplen=5`)
    is active; Ham-Song-style alternative kept under `#else`.
  - Removal guidance: **safe to simplify** (delete `#else` or split into two
    explicit tests).
  - Rationale: test-configuration convenience only; no production-path impact.
