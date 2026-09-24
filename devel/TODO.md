## Remaining questions about issue 2

~Symmetries ok?~

Remove ~ entries in tables

Add comments to the code.  Big audit of the code.

## Ask for a detailed latex file explaning the current state

On issue 2 branch.

## Plotting of train tracks.

Ask the AI for best option.  Goal is to make PDF, ultimately.  Vector graphics
preferred.

This is tough.  After a failed attempt, current suggestion is to use existing
graph plotting.

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
