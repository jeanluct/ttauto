env = Environment(CC = 'gcc',
                  CCFLAGS = ['-Wall','-O3','-ffast-math'],
		  LIBS = ['csparse', 'ttauto'],
                  LIBPATH = ['../CSparse', '../lib' ],
                  CPPPATH = ['../CSparse', '../include', '..'])

# Use profile=1 on the command-line to turn on profiling.
profile = ARGUMENTS.get('profile', 0)
if int(profile):
   env.PrependUnique(CCFLAGS = ['-pg'])
   env.PrependUnique(LINKFLAGS = ['-pg'])

# Use static=1 on the command-line for static linking.
static = ARGUMENTS.get('static', 0)
if int(static):
   env.PrependUnique(LINKFLAGS = ['-static'])

# Use command 'scons win=1 <program>.exe' to cross-compile for WIN32.
# Make sure to run 'scons -c' and 'scons win=1 -c' beforehand.
# On Linux, a cross-compiler such as mingw32 needs to be installed.
win = ARGUMENTS.get('win', 0)
if int(win):
   env.AppendUnique(CCFLAGS = ['-DTTAUTO_NO_BOOST'])
   env.Tool('crossmingw', toolpath = ['.'])

# For hash_set/hash_map with gcc >= 4.3.3.
# Eventually replace by unordered_set/unordered_map.
env.PrependUnique(CXXFLAGS = ['-Wno-deprecated'])

SConscript('CSparse/SConscript', exports = 'env')
SConscript('lib/SConscript', exports = 'env')

Export('env')
