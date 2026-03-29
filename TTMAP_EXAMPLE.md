3 monogons (= 3 peripheral edges) labeled 1,2,3
2 main edges labeled a, b

edge a goes from 1 to 2
edge b goes from 2 to 3

folding steps:

1) edge a folded onto b clockwise at the only cusp on monogon 2

ttmap for this fold is

a -> a -2 b
b -> -b
1 2 3 -> 1 3 2

1) edge b folded onto a counterclockwise at the only cusp on monogon 2

ttmap for this fold is

a -> -a
b -> a -2 b 
1 2 3 -> 2 1 3

The composition of these two operations is

a -> -b 2 -a
b -> a -2 b -3 -b

---

## Code-side fold index mapping (from `tests/test_freeword`)

For the `n=3` example (`trk=0`), one-step maps at initial vertex are:

- `f=0` (counterclockwise)
  - `1 -> 1 5 2`
  - `2 -> -2`
- `f=1` (clockwise)
  - `1 -> -1`
  - `2 -> 1 5 2`

Here `1,2` are main generators and `5` is the selected infinitesimal one.
This gives a concrete bridge between hand labels (`a,b,1,2,3`) and fold indices.

TODO:
- Confirm whether your step (1) should match code `f=1` or `f=0` after label
  relabeling (`a,b` vs generator `1,2`).
- Confirm the intended second step index and orientation on this same track,
  then lock expected composition in test assertions.

Status from `tests/test_freeword`:
- Step (1) is now locked: `f=1` gives
  - `1 -> -1`
  - `2 -> 1 5 2`
- Candidate step (2) as `f=0` after step (1) currently gives
  - `1 -> 1 5 2`
  - `2 -> -2`
This matches your current guess and is now locked in assertions.

With the current orientation convention, the two-step composition `[1,0]` is:

- `1 -> -2 5 -1`
- `2 -> 1 5 2 5 -2`
