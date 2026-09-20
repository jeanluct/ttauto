## Remaining questions about issue 2

~Symmetries ok?~

Remove ~ entries in tables

Add comments to the code.  Big audit of the code.

## Ham Song

Ask AI: incorrect?  What follows from it?  Counterexample.

## Ask for a detailed latex file explaning the current state

On issue 2 branch.

## Plotting of train tracks.

Ask the AI for best option.  Goal is to make PDF, ultimately.  Vector graphics
preferred.

This is tough.  After a failed attempt, current suggestion is to use existing
graph plotting.

## Compile out or turn off the 5th digit of coding?

Will this break things?  The `ttauto` paper right now only refers to 4 digits.

Update: there is already a flag in the traintrack class:
```
  static const bool label_multiprongs = true;
```
I checked that everything still works fine when set to false, except for one
unit test that tries to call set_label.  I moved the bool to the public
section and skip the test of set to false.

Does the string input respect this flag?  I don't think so.  It seems like
only the print_coding method suppresses the middle label digit.

Ask the AI to help improve this.

## Canonicalization issue 12

Move canonicalization to separate file? [DONE]
It's such a crucial part of the codebase that it's good to isolate it.
[NOW DONE BETTER IN MASTER; WILL BE A CHALLENGE TO MERGE BACK INTO ISS012]

Consistent 0-indexing in test file

## Zero-indexing?

Get rid of any reference to 1-indexed quantities?  It would actually be kind
of ugly, since the coding has 2 "m of n" units.  Saying "0 of 2" and "1 of 2"
instead of "1 of 2" and "2 of 2" seems weird.

## Parallel path search?

No real reason not to.  It's embarassingly parallel.  Only the display would
have to be modified in some way.
