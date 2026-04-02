Ok, I need a refresher on what's being output by ttauto.cpp.  Two m-files can get created, one with the pA's, the other ones with the graph.  What changes when we have the map_only flag?

Great question. Here’s the quick refresher.

ttauto.cpp can emit two Mathematica files:
- ttauto_...m (optional; when you answer “Write subgraphs to file = yes”)
  - contains the automaton graph/subgraph data.
- pA_...m (always when you run search with max path length > 0)
  - contains pseudo-Anosov candidate output from the search.

What changes with map_only vs matrix_only:
For ttauto_...m (graph file):
- First section (vertex codings + traintrack graph drawing) is unchanged.
- Second section (branch payload) changes:
  - matrix_only: each edge stores a permutation(+1) matrix payload.
  - map_only: each edge stores a free-group automorphism map payload, e.g.
    - {1->2,{1 -> -2, 2 -> 2 7 1, ...}}
    - i.e. source->target, {generator -> word, ...}

For pA_...m (search output):
- Intended mathematical output (candidates, polynomials/dilatations/paths, etc.) should be the same.
- Internally, with map_only, matrix data may be recovered from maps when needed instead of being pre-stored.
So practically:
- map_only mainly changes the payload representation in the optional ttauto_...m graph file.
- The pA result file should remain equivalent in content, barring ordering/noise differences
