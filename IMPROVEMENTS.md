# IMPROVEMENTS.md

This file lists non-blocking improvement ideas observed while working on the
repository. None of these are required for correctness right now; they are
candidate follow-ups.

## High-value, low-risk follow-ups

- Continue incremental `typedef` -> `using` modernization in headers.
  - Keep it mechanical and scoped per file.
  - Avoid mixing with behavior changes.

- Clean up stale `#if 0` blocks where code is clearly obsolete.
  - `include/traintracks/edge.hpp`
  - `include/traintracks/traintrack.hpp`
  - `include/ttauto/ttauto.hpp`
  - Remove only when comments indicate no intended future use.

- Tighten message consistency in fail-fast diagnostics.
  - Example: comparison errors in `path` currently mention `operator=`.
  - Keep existing fail-fast style (`std::cerr` + `std::exit(1)`) but improve
    message accuracy.

- Expand `TESTING.md` with rough expected runtime for long-running examples.
  - Keep recommendations practical for local dev and CI.

## Ownership and mutation hardening (optional)

- Continue narrowing direct edge metadata manipulation.
  - Current progress moved relink/renumber paths into small edge helpers.
  - Remaining work should preserve semantics and stay incremental.

- Add a sanitizer profile for optional debugging runs.
  - ASan/UBSan builds can improve confidence in mutation-heavy paths.
  - Keep off by default; document one command for local use.

- Consider one additional stress test focused on repeated swap/fold sequences
  with bounded steps and invariant checks.

## Build and toolchain

- Keep SCons C++ baseline at C++17 unless a compatibility target requires less.

- If cross-platform support is active, verify static CSparse linking behavior on
  non-Linux toolchains and document any differences.

## Nice-to-have refactors (defer unless needed)

- Evaluate whether `path` should eventually use composition instead of
  inheritance from `jlt::vector<int>`.
  - This is a larger API design change and should be done only with tests in
    place.

- If future ownership redesign is desired (e.g., ID-based backlinks), open a
  dedicated issue and proceed test-first in small commits from stable `master`.
