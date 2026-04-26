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

No CTest harness is required in the current workflow; test programs are run directly.

## Programs that can be used as tests but take a long time to run

`examples/ttauto_count`   (2 minutes)
`examples/ttauto_labels`   (100 minutes)

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
