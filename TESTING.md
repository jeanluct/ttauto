# How to test the ttauto code

## Programs that can be used as tests but take a long time to run

`examples/ttauto_count`   (2 minutes)
`examples/ttauto_label`   (? minutes)

## Programs that should be used for testing

`tests/test_permplus1`
`tests/test_traintrack`
`tests/test_folding_path`
`tests/test_badwords`
`examples/ttauto_min_examples`
`examples/ttauto_torus`

The program `examples/ttauto` has interactive input but the defaults can be
used:
```
examples/ttauto << EOF
EOF
```
This procuces a file `pA_n=5_1_1_inv.m`.
