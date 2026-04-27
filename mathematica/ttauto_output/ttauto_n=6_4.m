{
(*
Train track graph for n=6 monogons and 6 edges.

Initial train track:

multigon 0 is a monogon connected to 1 edge
multigon 1 is a monogon connected to 1 edge
multigon 2 is a monogon connected to 1 edge
multigon 3 is a monogon connected to 1 edge
multigon 4 is a monogon connected to 1 edge
multigon 5 is a monogon connected to 1 edge
multigon 6 is a pentagon with edge sequence 2 1 1 1 1 

Subgraph has 3 vertices.
*)
{
 {
  {"1111 1511 2511 1111 3511 1111 4511 1111 5512 1111 5522 1111",{1->8,2->7,3->7,4->9,5->10,6->11,7->8,7->11,8->9,9->10,10->11}},
  {"1111 1112 1122 1511 2511 1111 3511 1111 4511 1111 5511 1111",{1->6,2->8,3->9,4->10,5->11,6->7,7->8,7->11,8->9,9->10,10->11}},
  {"1111 1122 1112 1511 2511 1111 3511 1111 4511 1111 5511 1111",{1->6,2->8,3->9,4->10,5->11,6->7,7->8,7->11,8->9,9->10,10->11}}
 },
 {
  {1->2,{{5,6,1,2,3,4},{2,5}}},
  {1->3,{{6,5,1,2,3,4},{2,6}}},
  {2->1,{{4,5,6,2,1,3},{4,1}}},
  {2->2,{{1,2,3,4,5,6},{1,2}}},
  {3->3,{{1,2,3,4,5,6},{1,2}}},
  {3->1,{{2,3,4,5,6,1},{1,1}}}
 }
}
}
