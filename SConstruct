# <LICENSE
#   ttauto: a C++ library for building train track automata
#
#   https://github.com/jeanluct/ttauto
#
#   Copyright (C) 2010-2014  Jean-Luc Thiffeault   <jeanluc@math.wisc.edu>
#                            Erwan Lanneau <erwan.lanneau@ujf-grenoble.fr>
#
#   This file is part of ttauto.
#
#   ttauto is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   ttauto is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with ttauto.  If not, see <http://www.gnu.org/licenses/>.
# LICENSE>

import subprocess
from distutils.version import StrictVersion

# A '#' means relative to the root directory.
env = Environment(CC = 'gcc',
                  CCFLAGS = ['-Wall','-O3','-ffast-math'],
                  LIBS = ['csparse', 'ttauto'],
                  LIBPATH = ['#lib','#extern/CSparse'],
                  CPPPATH = ['#include','#extern/CSparse','#extern/jlt'])

# Use profile=1 on the command-line to turn on profiling.
profile = ARGUMENTS.get('profile', 0)
if int(profile):
   env.PrependUnique(CCFLAGS = ['-pg'])
   env.PrependUnique(LINKFLAGS = ['-pg'])

# Use static=1 on the command-line for static linking.
static = ARGUMENTS.get('static', 0)
if int(static):
   env.PrependUnique(LINKFLAGS = ['-static'])

# Use command 'scons win32=1 <program>.exe' to cross-compile for WIN32.
# Make sure to run 'scons -c' and 'scons win32=1 -c' beforehand.
# On Linux, a cross-compiler such as mingw32 needs to be installed.
win32 = ARGUMENTS.get('win32', 0)
if int(win32):
   env.AppendUnique(CXXFLAGS = ['-DTTAUTO_NO_SHARED_PTR'])
   env.Tool('crossmingw', toolpath = ['./devel'])

GCC_VERSION = StrictVersion(
   subprocess.check_output([env['CXX'],'-dumpversion']))

if GCC_VERSION < StrictVersion('4.5'):
   # For hash_set/hash_map with.
   # Eventually replaced by unordered_set/unordered_map.
   env.AppendUnique(CXXFLAGS = ['-DTTAUTO_OLD_HASH'])
   if StrictVersion('4.3') <= GCC_VERSION:
      # For hash_set/hash_map with gcc >= 4.3.3.
      env.PrependUnique(CXXFLAGS = ['-Wno-deprecated'])
else:
   env.PrependUnique(CXXFLAGS = ['-std=c++0x'])

env.SConscript(dirs = ['lib','examples','tests','extern/CSparse'],
               exports = 'env')
