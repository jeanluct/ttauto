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

TTAUTOROOTDIR = Dir('.').abspath
TTAUTOEXTDIR = TTAUTOROOTDIR + '/extern'

env = Environment(CC = 'gcc',
                  CCFLAGS = ['-Wall','-O3','-ffast-math'],
		  LIBS = ['csparse', 'ttauto'],
                  LIBPATH = [
		  	  TTAUTOEXTDIR + '/CSparse',
		   	  TTAUTOROOTDIR + '/lib'],
                  CPPPATH = [
		  	  TTAUTOEXTDIR + '/CSparse',
			  TTAUTOROOTDIR + '/include',
			  TTAUTOEXTDIR + '/jlt'])

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
   env.AppendUnique(CCFLAGS = ['-DTTAUTO_NO_BOOST'])
   env.Tool('crossmingw', toolpath = ['./devel'])

# For hash_set/hash_map with gcc >= 4.3.3.
# Eventually replace by unordered_set/unordered_map.
env.PrependUnique(CXXFLAGS = ['-Wno-deprecated'])

SConscript('extern/CSparse/SConscript', exports = 'env')
SConscript('lib/SConscript', exports = 'env')

Export('env')
