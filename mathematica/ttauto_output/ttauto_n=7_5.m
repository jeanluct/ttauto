{
(*
Train track graph for n=7 monogons and 7 edges.

Initial train track:

multigon 0 is a monogon connected to 1 edge
multigon 1 is a monogon connected to 1 edge
multigon 2 is a monogon connected to 1 edge
multigon 3 is a monogon connected to 1 edge
multigon 4 is a monogon connected to 1 edge
multigon 5 is a monogon connected to 1 edge
multigon 6 is a monogon connected to 1 edge
multigon 7 is a hexagon with edge sequence 2 1 1 1 1 1 

Subgraph has 3 vertices.
*)
{
 {
  {"1111 1611 2611 1111 3611 1111 4611 1111 5611 1111 6612 1111 6622 1111",{1->9,2->8,3->8,4->10,5->11,6->12,7->13,8->9,8->13,9->10,10->11,11->12,12->13}},
  {"1111 1112 1122 1611 2611 1111 3611 1111 4611 1111 5611 1111 6611 1111",{1->7,2->9,3->10,4->11,5->12,6->13,7->8,8->9,8->13,9->10,10->11,11->12,12->13}},
  {"1111 1122 1112 1611 2611 1111 3611 1111 4611 1111 5611 1111 6611 1111",{1->7,2->9,3->10,4->11,5->12,6->13,7->8,8->9,8->13,9->10,10->11,11->12,12->13}}
 },
 {
  {1->2,{{6,7,1,2,3,4,5},{2,6}}},
  {1->3,{{7,6,1,2,3,4,5},{2,7}}},
  {2->1,{{4,5,6,7,2,1,3},{5,1}}},
  {2->2,{{1,2,3,4,5,6,7},{1,2}}},
  {3->3,{{1,2,3,4,5,6,7},{1,2}}},
  {3->1,{{2,3,4,5,6,7,1},{1,1}}}
 }
}
}
