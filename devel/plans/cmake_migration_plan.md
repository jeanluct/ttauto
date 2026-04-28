# CMake Migration Plan

## Goal

Migrate the main `ttauto` build from SCons to CMake with a workflow optimized for in-place development and branch-local feature testing.

## Agreed Constraints

1. Backward compatibility with the existing SCons flow is not a priority.
2. Default CMake behavior should place built binaries in source-adjacent folders (`examples/`, `tests/`, `lib/`) for in-place usage.
3. No formal `ctest` suite is required initially.
4. Test sources should be auto-discovered from `tests/*.cpp` so new branch-local test programs compile without editing `CMakeLists.txt`.
5. Keep SCons temporarily; remove it near the end once CMake workflow is stable.

## High-Level Approach

Use a staged migration:

- Introduce CMake first.
- Preserve current developer behavior (manual executable test runs, feature probes as small programs).
- Defer cleanup/removal of SCons until final stabilization.

## Phase 1: Introduce Top-Level CMake Build

Create a root `CMakeLists.txt` that:

- sets `CMAKE_CXX_STANDARD 17` (required),
- applies baseline compile flags matching project defaults (`-Wall -O3 -ffast-math` where appropriate),
- defines include paths:
  - `include/`
  - `extern/jlt`
  - `extern/jlt/extern/CSparse/Include`
- builds the core `ttauto` library from `lib/traintracks/*.cpp`,
- links CSparse statically.

## Phase 2: Build Programs Automatically

### Examples

- auto-discover `examples/*.cpp`,
- build one executable per source file,
- link each against `ttauto`.

### Tests

- auto-discover `tests/*.cpp` (not only `test_*.cpp`),
- build one executable per source file,
- link each against `ttauto`.

This keeps branch workflow simple: dropping a new `tests/*.cpp` file is enough for automatic compilation.

## Phase 3: In-Place Output Layout

Set output directories by target category:

- example executables -> `examples/`
- test executables -> `tests/`
- static libraries -> `lib/`

Rationale: preserves current "run from tree" workflow and minimizes friction during feature work.

## Phase 4: CSparse Integration Details

Prefer CMake target-based linking over hard-coded archive paths.

Options:

1. Add CSparse as a subdirectory and link to its static target.
2. If that proves brittle, build CSparse separately in a sub-build and import the static artifact.

Primary preference is (1), with explicit static linkage behavior in top-level options.

## Phase 5: Developer Documentation Update

Update `README.md` and `TESTING.md` to document CMake-first usage:

- configure/build commands,
- where binaries are produced,
- manual test execution flow using built `tests/*` executables,
- note that `tests/*.cpp` are auto-built by default.

Keep a short SCons note during transition.

## Phase 6: Stabilization and Final SCons Removal

After the CMake workflow is validated in regular use:

- remove SCons files and SCons-based docs,
- keep only CMake build instructions,
- verify fresh clone onboarding with CMake only.

## Non-Goals for Initial Migration

- no formal `ctest` orchestration yet,
- no conversion of current executable tests into a unit-test framework,
- no effort to preserve old command-line SCons options unless they are still actively needed.

## Risks and Mitigations

- **Risk:** `tests/*.cpp` glob may pick up non-standalone helper files.
  - **Mitigation:** keep only standalone programs in `tests/`; if helper `.cpp` appears later, move it under a non-globbed subdirectory.
- **Risk:** platform-specific linker differences (especially static CSparse).
  - **Mitigation:** test Linux first, then adjust target properties/options for other toolchains as needed.

## Success Criteria

Migration is considered successful when:

1. A fresh configure+build with CMake produces library, examples, and tests.
2. New files added under `tests/*.cpp` compile automatically with no CMake edits.
3. Built executables appear in `examples/` and `tests/` by default.
4. Core manual test flow from `TESTING.md` runs using CMake-built binaries.
5. SCons can be removed without losing required developer workflows.
