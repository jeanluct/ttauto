{
(*
Train track graph for n=4 monogons and 3 edges.

Initial train track:

multigon 0 is a monogon connected to 1 edge
multigon 1 is a monogon connected to 1 edge
multigon 2 is a monogon connected to 1 edge
multigon 3 is a monogon connected to 3 edges

Subgraph has 3 vertices.
*)
{
 {
  {"1111 1113 1123 1111 1133 1111",{1->4,2->4,3->4}},
  {"1111 1122 1112 1112 1122 1111",{1->3,2->4,3->4}},
  {"1111 1112 1122 1122 1112 1111",{1->4,2->3,3->4}}
 },
 {
  {1->2,{{2,1,3},{2,2}}},
  {1->3,{{1,3,2},{2,2}}},
  {2->2,{{3,2,1},{3,2}}},
  {2->1,{{2,3,1},{1,1}}},
  {2->2,{{1,2,3},{3,2}}},
  {2->1,{{2,1,3},{1,3}}},
  {3->1,{{1,3,2},{3,1}}},
  {3->3,{{3,2,1},{3,2}}},
  {3->1,{{3,1,2},{3,3}}},
  {3->3,{{1,2,3},{3,2}}}
 }
}
,
(*
Train track graph for n=4 monogons and 3 edges.

Initial train track:

multigon 0 is a monogon connected to 1 edge
multigon 1 is a monogon connected to 1 edge
multigon 2 is a monogon connected to 2 edges
multigon 3 is a monogon connected to 2 edges

Subgraph has 1 vertex.
*)
{
 {
  {"1111 1112 1122 1112 1122 1111",{1->3,2->4,3->4}}
 },
 {
  {1->1,{{1,2,3},{1,2}}},
  {1->1,{{1,2,3},{3,2}}}
 }
}
}
