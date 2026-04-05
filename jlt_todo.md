## Misc TODOs

freeauto_fold ?  (like mathmatrix_permplus1)  I wonder how much it would speed
things up in practice.

iss009: In traintrack.hpp, I'm worried it might be better to optionally find
the train track map.  I suspect computing the ttmap is expensive, so should be
done after the matrices are computed, using the path.  As things stand it
sounds like the train track map is always computed.  (Now issue #9)

Misusing infinitesimal edges I think.  Should use peripheral.  Will figure
this out as part of issue 2.

## DONE

matrix transposes

name? freeword freeauto
underscore
separate files

test_freeword / test_freeauto (no tt dep)

ttmap.hpp?

ttbuild -> traintrack_build

division between ttauto and traintracks

error messages scoping + namespaces

Is this comment correct in traintracks_util: fold_transition_matrix?
  // - Applying f1 then f2 corresponds to left-multiplication:
  //     TM(f2 followed by f1) = TM(f2) * TM(f1).

Add // comments describing functions and methods.

Is there anywhere where move semantics could help?  It might not be necessary
since we are mostly copying pointers.  The various swap functions?

auto keyword

Can you precede macros (including guard macros) by TTAUTO_ or TRAINTRACKS_ depending on the namespace/component.  Some macros already do this, but some should probably be rename from TTAUTO_ to TRAINTRACKS_ .

Which macros are enabled?  Code compiles with/without?

Verify content of #if 0/1 blocks.

The code a the start of traintrack/map.hpp looks like it could be in
mathmatrix_permplus1.

change year 2014 --> 2026 (Erwan?)

Why is there another perm+1 check in map.hpp: Check that the transition matrix
is permutation+1 or identity.  Seems redundant, again.

Why are we creating full mathmatrix objects for the folds?
fold_transition_matrix could return a mathmatrix_permplus1 object rather than
a jlt::mathmatrix<int> object.

Move printMathematicaForm out of class.

printMathematicaForm for freeauto.

### Move freeword/auto to jlt lib

Let's move freeword.hpp and freeauto.hpp to the jlt lib in extern.
* Read extern/jlt/AGENTS.md to know the local style.
* Move the hpp files to extern/jlt/jlt and update them to conform to the local style.  Examine other files in the same folder for guidance.  Make sure error handling is consistent with the jlt lib.
* Move test_free*.cpp to extern/jlt/examples.  Again update them to the local style.
* Update the scons build system in extern/jlt/examples.
* Update the includes in the ttauto project to <jlt/freeword.hpp>, etc.
* After you've done this, next step will be to add comprehensive tests for freeword and freeauto in the extern/jlt/tests folder.
