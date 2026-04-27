{
(*
Train track graph for n=4 monogons and 4 edges.

Initial train track:

multigon 0 is a monogon connected to 1 edge
multigon 1 is a monogon connected to 1 edge
multigon 2 is a monogon connected to 1 edge
multigon 3 is a monogon connected to 1 edge
multigon 4 is a trigon with edge sequence 2 1 1 

Subgraph has 3 vertices.
*)
{
 {
  {"1111 1311 2311 1111 3312 1111 3322 1111",{1->6,2->5,3->5,4->7,5->6,5->7,6->7}},
  {"1111 1112 1122 1311 2311 1111 3311 1111",{1->4,2->6,3->7,4->5,5->6,5->7,6->7}},
  {"1111 1122 1112 1311 2311 1111 3311 1111",{1->4,2->6,3->7,4->5,5->6,5->7,6->7}}
 },
 {
  {1->2,{{3,4,1,2},{2,3}}},
  {1->3,{{4,3,1,2},{2,4}}},
  {2->1,{{4,2,1,3},{2,1}}},
  {2->2,{{1,2,3,4},{1,2}}},
  {3->3,{{1,2,3,4},{1,2}}},
  {3->1,{{2,3,4,1},{1,1}}}
 }
}
}
