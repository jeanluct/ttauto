# Issue #14: large train-track automata crash during construction

Status: carried out 2026-09-20 on branch
`iss014-large-ttgraph-segfaults`, in four commits.  The previous branch
is preserved as `iss014-large-ttgraph-segfaults-dubious`; its
coding-traversal guard is not a fix and was not carried over.

**Result.**  The nine-puncture stratum-8 automaton, the reproducer for
this issue, now builds on the default 8 MB stack: 71253 vertices in
165 s, 1.4 GB peak.  Construction of all twelve seven-puncture strata
fell from 6.81 s to 0.73 s.  Every automaton for n = 3 to 7 is
byte-identical to before, codings, branch targets, matrices and maps
alike, so the vertex numbering quoted in the tests, the issue-2 note
and the papers is unchanged.

**The plan below missed the larger half.**  It identified the recursion
and the vertex lookup, and both needed fixing, but once they were done
the reproducer stopped allocating and sat at 1.4 GB inside
`find_symmetries()`, which compared every vertex's reversed coding
against every earlier vertex, recomputing codings for each pair.  At
1012 vertices that already cost five times as much as building the
graph.  It is the same quadratic pattern and took the same fix, in a
fourth commit.

On whether this is worth doing now: both changes matter almost entirely
beyond seven punctures, and essentially not at all below that, since
every automaton the papers rest on builds fine today.  The scoping
discussion is on the issue.

## The problem

Building the automaton for nine punctures, stratum 8
(`1. 1. 1. 1. 1. 1. 1. 1. 1. 3 3 (5)`) segfaults.  Two runs on master
`90f69a6`, concurrent, same deterministic input, differing only in the
stack limit:

| run | stack limit | outcome |
|---|---|---|
| A | 8 MB (default) | SIGSEGV after 11 min 41 s, 237.7 MB max RSS |
| B | `ulimit -s unlimited` | no crash; still building at a 15 min cap, 314.6 MB max RSS |

At 11 min 36 s run B held 237.1 MB, the same point where run A died five
seconds later.  The cause is stack exhaustion in
`ttfoldgraph::add_vertex`, which calls itself once per newly discovered
vertex, so the C stack depth is bounded only by the vertex count.  Each
live frame carries a full train-track copy (`trtr0`) and the fold's
free-group map (`AM`), so the depth ceiling is low.

Backtraces land in `recursive_coding` only because that is the innermost
code running when the stack runs out; every frame ends in
`fold` -> `normalise` -> `minimise_coding` -> `coding_from_monogon`.

Raising the stack is a stopgap: it moves the threshold rather than
removing it, it is not under the program's control (`ulimit` belongs to
the caller's shell), it fails silently as a bare SIGSEGV, and it leaves
the second problem untouched.  That second problem is the vertex lookup:
`std::find` over `trtrv` compares with `operator==`, which recomputes
codings on both sides, so construction is quadratic in the vertex count.
This is where the eleven minutes go.

## What the code does today

In `include/ttauto/ttfoldgraph.hpp`:

- `add_vertex` is lines 128-186, 59 lines, private, with exactly two call
  sites: the first constructor (line 108) and its own recursion (line
  178).  The second constructor takes prebuilt vectors and never calls
  it.
- `delete_vertex` and `delete_vertices` have no callers anywhere in the
  repository; both are marked "no longer in use".
- `swap_vertices` is called only from `sort_by_symmetries`, which runs
  inside `find_symmetries()` after the build has finished.

Two consequences shape the fix.  A coding-to-index map is needed only
while the graph is being built, so it can be a local of the build
routine and never has to survive deletion or renumbering.  And an
explicit stack frame needs only `(vertex index, next fold index)`,
because the track itself is already stored in `trtrv[idx]`.

## Step 1: replace the recursion with an explicit worklist

One commit.  `include/ttauto/ttfoldgraph.hpp` only.

- Replace `add_vertex` with `void build_graph(const TrTr& start)` and a
  small private helper `int new_vertex(const TrTr&)` holding the five
  `push_back` calls that currently create a vertex.
- The worklist is a `std::vector` of frames, each `{int idx; int f;}`,
  used as a stack.  Seed it with `new_vertex(start)`.
- Loop: look at the top frame; if its `f` has reached `nfoldsmax`, pop
  and continue; otherwise take `f` and increment it, copy `trtrv[idx]`,
  fold the copy, compute the matrix, and if the matrix is not the
  identity record the branch, look up or create the target, push the
  target index onto `tv[idx]`, and if the target is new push its frame.
- Update the call site at line 108.  Keep the `debug` prints.

### Why the vertex numbering is preserved

The current numbering is pre-order depth-first.  The worklist reproduces
it exactly because the child's frame is pushed immediately after its
branch is recorded, so LIFO order processes that child and its whole
subtree before the parent's next fold, and because each vertex receives
exactly one frame (a second encounter finds it in the index and pushes
no frame).

The one reordering is safe: today `tv[idx].push_back(tidx)` happens after
the recursive call returns, whereas the worklist pushes it before the
child is explored.  These are equivalent, because the exploration of a
child can only append to `tv` of other vertices; a second visit to `idx`
returns early without touching `tv[idx]`.  Branch order within a vertex
stays fold order in both forms.

### Hazards to respect

- `trtrv` grows during the build, so never hold a reference or pointer
  into it across a `push_back`.  Take the copy of the source track and
  compute `AM` and the matrix before creating the target vertex, exactly
  as the recursion does.
- The same applies to the frame: `stack.back()` is invalidated by
  pushing a new frame, so read `idx` and `f` by value first.
- `TM` is a member used as scratch.  It is consumed before the child is
  created, as it is today, but making it a local of the loop would
  remove the hazard entirely and is worth doing while the function is
  being rewritten.

## Step 2: replace the linear lookup with a coding-to-index map

Second commit, same file.

- Add a hash functor for the coding vector (`TrTr::intVec`, which
  derives from `std::vector<int>`), combining the elements in the usual
  way.  The existing `folding_path<TrTr>::hash` is the style reference,
  though it is a plain sum and a proper combine is preferable here.
- Hold `std::unordered_map<intVec,int,coding_hash>` as a local of
  `build_graph`, inserting each vertex as it is created and querying it
  in place of `std::find`.
- The start track must be normalised for its coding to be canonical;
  folded tracks are normalised by the fold itself.

### Why keying on the coding is faithful

`operator==` is `same_multigons()` and coding equality.  Keying on the
coding alone is equivalent provided that distinct vertices always have
distinct codings, which the tests below assert directly.  Note that
issue #12 concerned the opposite failure, isotopic tracks receiving
*different* codings, which would create duplicate vertices; that is
closed and is unaffected by this change, since the map preserves the
current equality semantics exactly.

### What it buys

Each lookup currently computes up to `V` codings on the stored side plus
one per comparison on the query side; the map computes one coding per
candidate track.  Construction goes from quadratic to near-linear in the
vertex count.

## Tests

New file `testsuite/ttauto/test_ttfoldgraph_build.cpp`, with
`CHECK`/`CHECK_MSG` from `testsuite/check.hpp`:

1. **The recursion is gone.**  Build a mid-sized automaton (an n=7
   stratum) on a POSIX thread created with a reduced stack via
   `pthread_attr_setstacksize`, and check it completes with the expected
   vertex count.  The stack size and stratum are chosen by measurement
   during implementation: shrink the stack until the pre-fix code fails,
   then set the test comfortably below that, with the post-fix code
   passing by a wide margin.  The iterative build needs only a few
   kilobytes of C stack, since the remaining recursion (the coding walk)
   has depth of order the number of multigons.
2. **Codings are unique per vertex.**  For each automaton of n=3..6,
   insert every vertex coding into a `std::set` and check the set size
   equals `vertices()`.  This is what makes the map faithful, and it
   would catch a regression of issue #12 as well.
3. **The graph is unchanged.**  Exact expected vertex counts and total
   branch counts per stratum for n=3..6, read off the current build.

`CMakeLists.txt` gains `find_package(Threads)` and links
`Threads::Threads` into that one test target.

The test in item 1 must be verified to fail against the pre-fix code
before it is trusted.  A test that cannot fail is worse than none, as
the `-DNDEBUG` episode in this repository showed.

## Verification

- `ctest --test-dir build --output-on-failure`, including the slow
  strata-scan regression, all green.  The existing suite already pins
  the vertex numbering this change must preserve: `test_gates`
  hard-codes the n=6 cycle 28/45/42/70/87/84 and the n=4 path
  0/2/0/3/0, and `test_ttauto_gates`, `test_min_dilatations`,
  `test_ttauto_search` and the scan regression all depend on exact
  results.
- After step 1, dump every automaton for n=3..7 (vertex codings, branch
  targets, transition matrices) with the old and the new code and diff
  them; they must be byte-identical.  A scratch program, not committed.
- Time `tests/ttauto_gate_census 7 7` before and after; it is about
  8.7 s today.
- Acceptance: build n=9 stratum 8 with the default 8 MB stack and have
  it complete.  **Done**: 71253 vertices, 165 s, 1.4 GB peak.

## What was measured

| quantity | before | after |
|---|---|---|
| n=6 construction, all 7 strata (376 vertices) | 0.222 s | 0.073 s |
| n=7 construction, all 12 strata (3272 vertices) | 6.81 s | 0.73 s |
| smallest stack building n=7 stratum 6 | > 384 KB | < 32 KB |
| n=9 stratum 8 | SIGSEGV at 11 min | 71253 vertices in 165 s |

One correction to the earlier estimate: the vertex lookup was about
three quarters of construction time at n=7, not the factor of thousands
guessed from the asymptotics.  `operator==` tests `same_multigons()`
first, which is cheap and usually rejects, so the scan rarely reached
the coding comparison it was assumed to be dominated by.

## Risks

- **A silent renumbering.**  Mitigated by the diff of full graph dumps
  after step 1, before the lookup semantics change, and by the existing
  hard-coded indices in the suite.
- **Coding collisions merging two vertices.**  Mitigated by test 2.  If
  it ever fired, the map would have to key on the coding and fall back
  to `operator==` on collision.
- **Memory rather than stack becomes the limit.**  We still do not know
  how many vertices stratum 8 has; run B had not finished at fifteen
  minutes.  The crash may turn into a very large allocation.  The
  acceptance run will tell us, and the answer belongs in this file.
  Note that step 2 makes this slightly worse: the map holds a coding per
  vertex, on the order of a few hundred bytes at nine punctures.  That
  is the right trade unless memory is what stops us, which is precisely
  the unknown here.
- **The `pthread` dependency** in one test target, on a project that is
  otherwise plain C++17 plus CMake.  The repository has cross-MinGW
  history, visible in the Windows entries of `.gitignore`, so guard the
  reduced-stack test for the case where `find_package(Threads)` comes up
  empty rather than letting it become the first thing that breaks such a
  build.

## Out of scope

- Caching the minimal coding computed inside `normalise()` so that
  `coding()` need not recompute it.  A further constant-factor win,
  worth its own issue.
- Removing the dead `delete_vertex` and `delete_vertices`.
- The diagnostics gating from the `-dubious` branch; if progress
  reporting is wanted for long builds, add it after the rewrite.
