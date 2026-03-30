# <LICENSE
#   ttauto: a C++ library for building train track automata
#
#   https://github.com/jeanluct/ttauto
#
#   Copyright (C) 2010-2026  Jean-Luc Thiffeault   <jeanluc@math.wisc.edu>
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

# Choose a particular compiler by specifying the version here.
#GCC_version = StrictVersion('9.3.0');

cc = 'gcc'
cxx = 'g++'

if 'GCC_version' not in locals():
   # If we didn't define a GCC version number, find out what it is.
   GCC_version = StrictVersion(
	       subprocess
	       .check_output([cxx,'-dumpfullversion'])
	       .decode('utf-8'))
else:
   # If we define a version number, call the appropriate compiler.
   cc = cc + '-' + str(GCC_version)
   cxx = cxx + '-' + str(GCC_version)

# A '#' means relative to the root directory.
jltdir = '#extern/jlt'
jltincdir = jltdir
csparsedir = jltdir + '/extern/CSparse'
csparselibdir = csparsedir + '/build'  # make sure to link to static lib below
csparseincdir = csparsedir + '/Include'
env = Environment(CC = cc, CXX = cxx,
                  CCFLAGS = ['-Wall','-O3','-ffast-math'],
                  LIBS = ['ttauto',File(csparselibdir + '/libcsparse.a')],
                  LIBPATH = ['#lib'],
                  CPPPATH = ['#include',csparseincdir,jltincdir])

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
   env.AppendUnique(CXXFLAGS = ['-DTRAINTRACKS_NO_SHARED_PTR'])
   env.Tool('crossmingw', toolpath = ['./devel'])

# Require a modern C++ standard across all builds.
env.PrependUnique(CXXFLAGS = ['-std=c++17'])

# Keep toolchain expectations explicit.
if GCC_version < StrictVersion('7.0'):
   print("Error: GCC >= 7.0 is required (C++17 baseline).")
   Exit(1)

env.SConscript(dirs = ['lib','examples','tests',csparsedir],
               exports = 'env',
               must_exist = False)
