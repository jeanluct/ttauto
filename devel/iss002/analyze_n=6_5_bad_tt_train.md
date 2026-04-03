From `n=6_5_bad_tt.train`:

 1 ->  6
 2 ->  1
 3 ->  2
 4 ->  4
 5 ->  5
 6 ->  7
 7 ->  3

 8 -> -13
 9 ->  13   8 -3 9
10 ->  10 -14 -6 14
11 -> -10  -4 -9
12 ->  14 -10 -4 -9 3 -8
13 ->  12
14 -> -11

Derivative map:

  8 -> -13
  9 ->  13
 10 ->  10
 11 -> -10
 12 ->  14
 13 ->  12
 14 -> -11
-8  ->  13
-9  -> -9
-10 -> -14
-11 ->  9
-12 ->  8
-13 -> -12
-14 ->  11

Example:

Edges at vertex 1: 1 -1 12

Gates are {-1}, {12}, {1} which are single-element equivalence classes.


Edges at vertex 5: -3 -8 9 3

Gate {-8, 9}: because -8 -> 13 and 9 -> 13 (k=1)


Edges at vertex 2: -4 10 5
Gates are {-4} {10} {5}
Infinitesimal edges join -4 to 10 (BUT NOT 5)


Edges at vertex 6: -9 4 -5
Gates are {-9} {4} {-5}
Infinitesimal edges join -9 to 4 (BUT NOT 5)


So the violation of B&H 3.3.2 is for the gates at vertices 2 and 6.




