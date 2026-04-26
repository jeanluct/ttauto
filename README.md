# ttauto

*ttauto* is a C++ library for building train track automata for homeomorphisms of punctured discs.  It was written by [Jean-Luc Thiffeault][1] and [Erwan Lanneau][2].

### external code

*ttauto* uses Timothy A. Davis's [CSparse][3] and the [jlt][4] library.  These are included with the source code in the folder extern, along with their license information.

### license

*ttauto* is released under the [GNU General Public License v3][5].  See [COPYING][6] and [LICENSE][7].

### documentation

There is currently no real documentation for *ttauto*.  See the [examples folder][8] for some basic examples.  The most complete program is [ttauto.cpp][9], an interactive program.

### build (cmake)

The project targets C++17 by default and now builds with CMake.

From the repository root:

```bash
cmake -S . -B build
cmake --build build -j
```

Default output behavior is in-place:

- examples are written to `examples/`
- test programs are written to `tests/`
- the static library is written to `lib/`

Both `examples/*.cpp` and `tests/*.cpp` are auto-discovered by CMake.  Adding a new `.cpp` file in either folder is enough for it to compile on the next build, with no CMake file edits required.

### build (legacy scons)

SCons files are still present temporarily during migration:

```bash
scons
```

The CMake workflow above is the primary path.

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
