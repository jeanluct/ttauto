# Macro Status (Non-Guard Macros)

This file tracks build/runtime feature macros used in the main `ttauto` codebase
(excluding include guard macros).

## Core / traintracks-layer macros

- `TRAINTRACKS_NO_SHARED_PTR`
  - Default: **not defined**
  - Enabled when: `scons win32=1` (set in `SConstruct`)
  - Effect: switches core pointer aliases from shared pointers to raw pointers in
    compatibility paths.
  - Status: **supported compatibility mode** (kept intentionally for portability).

`TRAINTRACKS_OLD_HASH` has been removed.

## Automaton / ttauto-layer macros

- `TTAUTO_HASH_BADWORDS`
  - Default: **not defined**
  - Note: currently forced off by `#undef TTAUTO_HASH_BADWORDS` in
    `include/badwords.hpp`.
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

Yes - this is easy to test with `scons` by passing `CPPFLAGS`.

Examples:

```bash
# raw-pointer compatibility path
scons -c tests/test_traintrack
scons CPPFLAGS='-DTRAINTRACKS_NO_SHARED_PTR' tests/test_traintrack

# optional symmetric norm checks in automaton layer
scons -c tests/test_badwords
scons CPPFLAGS='-DTTAUTO_CHECK_SYMMETRIC_NORM' tests/test_badwords

# optional Fortran-enabled path (compile-time check)
scons -c tests/test_badwords
scons CPPFLAGS='-DTTAUTO_USE_FORTRAN' tests/test_badwords
```

Notes:

- Prefer `CPPFLAGS` for preprocessor defines; overriding `CXXFLAGS` can accidentally
  replace required default flags.
- Build one or a few representative targets per macro mode for quick checks,
  then scale to full `scons` if needed.
