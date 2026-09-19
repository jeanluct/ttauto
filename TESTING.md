# How to test the ttauto code

## Build first (CMake)

From the repository root:

```bash
cmake -S . -B build
cmake --build build -j
```

This compiles:

- all `examples/*.cpp` into `examples/`
- all `tests/*.cpp` into `tests/`
- all `testsuite/**/*.cpp` into `build/testsuite/`

`ctest` runs deterministic fast tests by default.

Testsuite programs are compiled with `-UNDEBUG` (see `ttauto_add_testsuite`
in `CMakeLists.txt`), so their `assert()` calls stay live in the default
Release build.  Before this was added every assert-based test passed
vacuously.  New tests should prefer `CHECK(cond)` / `CHECK_MSG(cond, msg)`
from `testsuite/check.hpp`, which print the failing expression and its
location and exit nonzero.

```bash
ctest --test-dir build --output-on-failure
```

To include slow integration tests, configure with:

```bash
cmake -S . -B build -DTTAUTO_ENABLE_SLOW_TESTS=ON
ctest --test-dir build --output-on-failure
```

To run only slow tests:

```bash
ctest --test-dir build --output-on-failure -L slow
```

The `tests/` folder remains available for direct-run feature and exploratory
programs.

For header/component coverage mapping, see `testsuite/COVERAGE.md`.

## Programs that can be used as tests but take a long time to run

`examples/ttauto_count`   (2 minutes)
`examples/ttauto_labels`   (100 minutes)

## Measuring line coverage

Use a separate Debug build with gcov instrumentation; never the in-place
`build/` directory:

```bash
cmake -S . -B build-cov -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_C_FLAGS=--coverage -DCMAKE_CXX_FLAGS=--coverage \
      -DCMAKE_EXE_LINKER_FLAGS=--coverage
cmake --build build-cov -j --target $(cmake --build build-cov --target help \
      | grep -o 'testsuite_[a-z_0-9]*' | sort -u)
ctest --test-dir build-cov -LE slow
gcov --json-format -o build-cov/CMakeFiles/ttauto.dir/lib/traintracks \
      build-cov/CMakeFiles/ttauto.dir/lib/traintracks/*.gcda
```

Build only the `testsuite_*` targets: the examples and `tests/` programs
are placed in-place in `examples/` and `tests/` and would be overwritten
with instrumented binaries.  The library archive `lib/libttauto.a` is also
in-place and *is* overwritten, so after a coverage run do

```bash
rm -rf build-cov lib/libttauto.a && cmake --build build -j
```

to restore the Release archive.  Header-only code (`include/ttauto/**`)
is covered through the testsuite objects
(`build-cov/CMakeFiles/testsuite_*.dir/`); merge the per-object JSON by
line to get a per-file figure.  The per-file table as of 2026-09-19 is at
the top of `testsuite/COVERAGE.md`.

## Gate-test census

`tests/ttauto_gate_census [nmin [nmax]]` (default 3 6; n=7 takes about
9 s) repeats the strata scan with the Bestvina-Handel gate test on and lists
every class rejected by it, per stratum, with dilatation, characteristic
polynomial and representative paths.  It also counts, among the stored
representatives of the accepted classes, those whose gate partition at an
unpunctured multigon is finer than its prongs (zero on every stratum up to
n=7 as of 2026-09-19).  Use it after any change to the fold map or the gate
test; the expected rejections are recorded in
`devel/iss002/ISSUE2_IMPLEMENTATION_PLAN.md` (Step 4 outcome).

## Programs that should be used for testing

`tests/test_permplus1`
`tests/test_traintrack`
`tests/test_folding_path`
`tests/test_badwords`
`examples/ttauto_min_example`
`examples/ttauto_torus`

The program `examples/ttauto` has interactive input but the defaults can be
used:
```
examples/ttauto << EOF
EOF
```
This produces a file `pA_n=5_1_1_inv.m`.
