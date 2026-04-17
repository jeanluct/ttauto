# ttauto

*ttauto* is a C++ library for building train track automata for homeomorphisms of punctured discs.  It was written by [Jean-Luc Thiffeault][1] and [Erwan Lanneau][2].

### external code

*ttauto* uses Timothy A. Davis's [CSparse][3] and the [jlt][4] library.  These are included with the source code in the folder extern, along with their license information.

### license

*ttauto* is released under the [GNU General Public License v3][5].  See [COPYING][6] and [LICENSE][7].

### documentation

There is currently no real documentation for *ttauto*.  See the [examples folder][8] for some basic examples.  The most complete program is [ttauto.cpp][9], an interactive program.  The programs are readily compiled with the [SCONS][10] build tool.

The project targets C++17 by default.  For compatibility builds that avoid
shared-pointer ownership in core train-track structures, you can compile with
`-DTRAINTRACKS_NO_SHARED_PTR` (for example via `scons win32=1`).

To compile everything, invoke `scons` in the root folder.  To build only a subset, use `scons -u` in a subfolder (such as [examples][8]).

### support

The development of *ttauto* was supported by the [US National Science Foundation][11], under grants [DMS-0806821][12] and [CMMI-1233935][13].

[1]: https://www.math.wisc.edu/~jeanluc/
[2]: https://www-fourier.ujf-grenoble.fr/~lanneau/
[3]: https://www.suitesparse.com
[4]: https://github.com/jeanluct/jlt
[5]: https://www.gnu.org/licenses/gpl-3.0.html
[6]: https://github.com/jeanluct/ttauto/raw/master/COPYING
[7]: https://github.com/jeanluct/ttauto/raw/master/LICENSE
[8]: https://github.com/jeanluct/ttauto/raw/master/examples
[9]: https://github.com/jeanluct/ttauto/raw/master/examples/ttauto.cpp
[10]: https://www.scons.org
[11]: https://www.nsf.gov
[12]: https://www.nsf.gov/awardsearch/showAward?AWD_ID=0806821
[13]: https://www.nsf.gov/awardsearch/showAward?AWD_ID=1233935
