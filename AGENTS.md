# AGENTS.md - Practical Guide for Coding Agents

## Scope and Priority

- This file applies to the repository root (`ttauto-code`).
- `extern/jlt/AGENTS.md` overrides this file for paths under `extern/jlt/`.
- Preserve behavior by default; treat broad rewrites as separate, explicit tasks.

## Repository Map

- `include/` - public headers for the `ttauto` core.
- `lib/` - `ttauto` implementation (`multigon`, `traintrack`, builders, etc.).
- `examples/` - executable programs used for smoke/integration checks.
- `tests/` - executable tests (plain `main()`, no unit-test framework).
- `extern/jlt/` - bundled dependency with its own build/tests and AGENTS file.
- `extern/jlt/extern/CSparse/` - bundled CSparse source/build directory.

## Build Commands (SCons)

- Build all:
  - `scons`
- Clean:
  - `scons -c`
- Build one target quickly:
  - `scons tests/test_traintrack`
  - `scons examples/ttauto_min_example`
- Optional modes:
  - `scons profile=1`
  - `scons static=1`
  - `scons win32=1`

## Test Commands (Main Repo)

- Recommended fast validation set:
  - `./tests/test_permplus1`
  - `./tests/test_traintrack`
  - `./tests/test_folding_path`
  - `./tests/test_badwords`
  - `./examples/ttauto_min_example`
  - `./examples/ttauto_torus`
- Interactive example with defaults:
  - `./examples/ttauto <<'EOF'`
  - `EOF`
- Long-running programs (use intentionally):
  - `./examples/ttauto_count` (~2 minutes)
  - `./examples/ttauto_labels` (~100 minutes)

See `TESTING.md` for the project-maintained recommended matrix and timing notes.

## External jlt / CSparse Notes

- For `extern/jlt` edits, follow `extern/jlt/AGENTS.md` and run its CMake/CTest flow.
- If CSparse archive is missing, build it with CMake:
  - `cmake -S extern/jlt/extern/CSparse -B extern/jlt/extern/CSparse/build`
  - `cmake --build extern/jlt/extern/CSparse/build`
- Current SCons links CSparse statically from:
  - `extern/jlt/extern/CSparse/build/libcsparse.a`
- Runtime note: examples/tests should not require `libcsparse.so` at runtime.

## Build/Toolchain Facts

- Default compiler flags include `-Wall -O3 -ffast-math`.
- Current default language level in `SConstruct` is `-std=c++17` (for modern GCC branch).
- Keep changes portable with existing compatibility macros:
  - `TTAUTO_NO_SHARED_PTR`
  - `TTAUTO_OLD_HASH`

## Code Style (Observed Conventions)

### Headers and includes

- Keep the license banner in touched source/header files.
- Prefer include order:
  1) standard headers,
  2) external/project-wide (`<jlt/...>`),
  3) local (`"..."`).
- Avoid adding unused includes.

### Naming and layout

- Follow existing file-local style; do not normalize everything globally.
- Existing style is mixed but stable:
  - types/classes: lowercase/descriptive (`traintrack`, `multigon`)
  - methods/functions: snake_case (`print_coding`, `cycle_prongs`)
  - macros: `ALL_CAPS`
- Keep brace/indent style of the file; many files use compact 2-space style.

### Error handling and invariants

- Core code often uses fail-fast checks (`std::cerr` + `std::exit(1)`) for internal invariants.
- Return `false` for operational failure paths where existing APIs already do so.
- Avoid introducing exceptions in legacy paths unless explicitly requested.

## Ownership Model (Important)

Current stabilized model in `ttauto`:

- `traintrack` owns multigons (`mgv`, currently `std::unique_ptr<multigon>` in modern mode).
- `multigon` stores incident edges via `edgep` (shared/raw depending on macro mode).
- `edge` stores non-owning back-links to multigons.

Guidance:

- Do not introduce broad ownership rewrites in routine tasks.
- Prefer small helper refactors that keep semantics unchanged.
- Add or keep local invariant checks around relink/swap paths when touching them.

## Refactor Strategy That Worked Here

When changing pointer/relink logic in this repo, this incremental pattern was effective:

1. Document ownership intent in headers.
2. Add postcondition checks in mutation/relink paths.
3. Replace direct low-level field writes with tiny edge helpers.
4. Add bounded mutation stress in `tests/test_traintrack`.
5. Rebuild and run the full recommended test matrix.

Avoid: large one-shot pointer-model conversions across `edge`/`multigon`/`traintrack`.

## Testing Expectations for Changes

- Header/lib logic change (`include/`, `lib/`): run at least one relevant test binary.
- Mutation or graph-topology change: run all four core tests plus one example smoke test.
- Build-system change: run full `scons` and at least one example executable.
- `extern/jlt` change: run its local CMake/CTest targets.

## Git and Workflow Expectations

- Keep commits small, scoped, and descriptive.
- Do not revert unrelated working-tree changes.
- Avoid destructive git commands unless explicitly requested.
- If generated files appear during testing, clean them before committing.
