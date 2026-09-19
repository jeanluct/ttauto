// <LICENSE
//   ttauto: a C++ library for building train track automata
//
//   https://github.com/jeanluct/ttauto
//
//   Copyright (C) 2010-2026  Jean-Luc Thiffeault   <jeanluc@math.wisc.edu>
//                            Erwan Lanneau <erwan.lanneau@ujf-grenoble.fr>
//
//   This file is part of ttauto.
//
//   ttauto is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   ttauto is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with ttauto.  If not, see <http://www.gnu.org/licenses/>.
// LICENSE>

// Minimal check helper for testsuite programs.
//
// The testsuite is built in Release mode with -UNDEBUG (see CMakeLists.txt)
// so plain assert() is live as well; prefer CHECK() in new tests because it
// prints the failing expression with its location and exits nonzero without
// aborting, which reads better under ctest --output-on-failure.

#ifndef TTAUTO_TESTSUITE_CHECK_HPP
#define TTAUTO_TESTSUITE_CHECK_HPP

#include <cstdlib>
#include <iostream>

#define CHECK(cond)                                                     \
  do {                                                                  \
    if (!(cond))                                                        \
      {                                                                 \
        std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK failed: "  \
                  << #cond << "\n";                                     \
        std::exit(1);                                                   \
      }                                                                 \
  } while (0)

#define CHECK_MSG(cond, msg)                                            \
  do {                                                                  \
    if (!(cond))                                                        \
      {                                                                 \
        std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK failed: "  \
                  << #cond << " (" << msg << ")\n";                     \
        std::exit(1);                                                   \
      }                                                                 \
  } while (0)

#endif // TTAUTO_TESTSUITE_CHECK_HPP
