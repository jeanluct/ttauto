{
(*
Train track graph for n=5 monogons and 5 edges.

Initial train track:

multigon 0 is a monogon connected to 1 edge
multigon 1 is a monogon connected to 1 edge
multigon 2 is a monogon connected to 1 edge
multigon 3 is a monogon connected to 1 edge
multigon 4 is a monogon connected to 1 edge
multigon 5 is a tetragon with edge sequence 2 1 1 1 

Subgraph has 3 vertices.
*)
{
 {
  {"1111 1411 2411 1111 3411 1111 4412 1111 4422 1111",{1->7,2->6,3->6,4->8,5->9,6->7,6->9,7->8,8->9}},
  {"1111 1112 1122 1411 2411 1111 3411 1111 4411 1111",{1->5,2->7,3->8,4->9,5->6,6->7,6->9,7->8,8->9}},
  {"1111 1122 1112 1411 2411 1111 3411 1111 4411 1111",{1->5,2->7,3->8,4->9,5->6,6->7,6->9,7->8,8->9}}
 },
 {
  {1->2,{{4,5,1,2,3},{2,4}}},
  {1->3,{{5,4,1,2,3},{2,5}}},
  {2->1,{{4,5,2,1,3},{3,1}}},
  {2->2,{{1,2,3,4,5},{1,2}}},
  {3->3,{{1,2,3,4,5},{1,2}}},
  {3->1,{{2,3,4,5,1},{1,1}}}
 }
}
}
