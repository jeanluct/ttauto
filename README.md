# ttauto

*ttauto* is a C++ library for building train track automata for homeomorphisms of punctured discs.  It was written by [Jean-Luc Thiffeault][1] and [Erwan Lanneau][2].

### contributors

*ttauto* uses Timothy A. Davis's [CSparse][3].

### license

*ttauto* is released under the [GNU General Public License v3][4].  See [COPYING][5] and [LICENSE][6].

### documentation

There is currently no real documentation for *ttauto*.  See the [examples folder][7] for some basic examples.  The most complete program is [ttauto.cpp][8], an interactive program.  The programs must be compiled with the [SCONS][9] build tool.

To compile everything, invoke `scons` in the root folder.  To build only a subset, invoke `scons -u` in a subfolder (such as [examples][7]).

### support

The development of *ttauto* was supported by the [US National Science Foundation][10], under grants [DMS-0806821][11] and [CMMI-1233935][12].

[1]: http://www.math.wisc.edu/~jeanluc/
[2]: http://www-fourier.ujf-grenoble.fr/~lanneau/
[3]: http://www.suitesparse.com
[4]: http://www.gnu.org/licenses/gpl-3.0.html
[5]: http://github.com/jeanluct/ttauto/raw/master/COPYING
[6]: http://github.com/jeanluct/ttauto/raw/master/LICENSE
[7]: http://github.com/jeanluct/ttauto/raw/master/examples
[8]: http://github.com/jeanluct/ttauto/raw/master/examples/ttauto.cpp
[9]: http://www.scons.org
[10]: http://www.nsf.gov
[11]: http://www.nsf.gov/awardsearch/showAward?AWD_ID=0806821
[12]: http://www.nsf.gov/awardsearch/showAward?AWD_ID=1233935
