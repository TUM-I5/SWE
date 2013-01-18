#! /usr/bin/python

# @file
# This file is part of SWE.
#
# @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
# @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
#
# @section LICENSE
#
# SWE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SWE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SWE.  If not, see <http://www.gnu.org/licenses/>.
#
#
# @section DESCRIPTION
#
# Builds the SWE code with several options.
#

# print the welcome message
print '****************************************'
print '** Welcome to the build script of SWE **'
print '****************************************'
print 'SWE Copyright (C) 2012'
print ''
print '  Technische Universitaet Muenchen'
print '  Department of Informatics'
print '  Chair of Scientific Computing'
print '  http://www5.in.tum.de/SWE'
print ''
print 'SWE comes with ABSOLUTELY NO WARRANTY.'
print 'SWE is free software, and you are welcome to redistribute it'
print 'under certain conditions.'
print 'Details can be found in the file \'gpl.txt\'.'
print ''

import os

#
# set possible variables
#
vars = Variables()

# read parameters from a file if given
vars.AddVariables(
  PathVariable( 'buildVariablesFile', 'location of the python file, which contains the build variables', None, PathVariable.PathIsFile )
)
env = Environment(variables=vars)
if 'buildVariablesFile' in env:
  vars = Variables(env['buildVariablesFile'])

# SWE specific variables
vars.AddVariables(
  PathVariable( 'buildDir', 'where to build the code', 'build', PathVariable.PathIsDirCreate ),

  EnumVariable( 'compiler', 'used compiler', 'gnu',
                allowed_values=('gnu', 'intel')
              ),

  EnumVariable( 'compileMode', 'mode of the compilation', 'release',
                allowed_values=('debug', 'release')
              ),

  EnumVariable( 'parallelization', 'level of parallelization', 'none',
                allowed_values=('none', 'cuda', 'mpi_with_cuda', 'mpi')
              ),

  EnumVariable( 'computeCapability', 'optional architecture/compute capability of the CUDA card', 'sm_20',
                allowed_values=('sm_10', 'sm_11', 'sm_12','sm_13',
                                'sm_20', 'sm_21', 'sm_22', 'sm_23' )
              ),

  BoolVariable( 'openGL', 'compile with OpenGL visualization', False),

  BoolVariable( 'openGL_instr', 'add instructions to openGL version (requires SDL_ttf)', False ),

  BoolVariable( 'writeNetCDF', 'write output in the netCDF-format', False ),

  BoolVariable( 'asagi', 'use ASAGI', False ),

  PathVariable( 'asagiInputDir', 'location of netcdf input files', '', PathVariable.PathAccept ),

  EnumVariable( 'solver', 'Riemann solver', 'augrie',
                allowed_values=('rusanov', 'fwave', 'augrie', 'hybrid', 'fwavevec')
              ),
                  
  BoolVariable( 'vectorize', 'add pragmas to help vectorization (release only)', True ),
                  
  BoolVariable( 'showVectorization', 'show loop vectorization (Intel compiler only)', False ),

  EnumVariable( 'platform', 'compile for a specific platform (Intel compiler only', 'default',
                allowed_values=('default', 'mic' )
              ),

  BoolVariable( 'xmlRuntime', 'use a xml-file for runtime parameters', False )
)

# external variables
vars.AddVariables(
  PathVariable( 'compilerPath', 'location of the C++ compiler', None, PathVariable.PathIsFile ),
  PathVariable( 'linkerPath', 'location of the C++ linker', None, PathVariable.PathIsFile ),
  PathVariable( 'cudaToolkitDir', 'location of the CUDA toolkit', None ),
  PathVariable( 'cudaSDKDir', 'location of the CUDA SDK', None),
  PathVariable( 'libSDLDir', 'location of libSDL', None),
  PathVariable( 'netCDFDir', 'location of netCDF', None),
  PathVariable( 'asagiDir', 'location of ASAGI', None),
  PathVariable( 'libxmlDir', 'location of libxml2', None)
)

# set environment
env = Environment(ENV = {'PATH': os.environ['PATH']},
        variables=vars)

# generate help text
Help(vars.GenerateHelpText(env))

# handle unknown, maybe misspelled variables
unknownVariables = vars.UnknownVariables()

# remove the buildVariablesFile from the list of unknown variables (used before)
if 'buildVariablesFile' in unknownVariables:
  unknownVariables.pop('buildVariablesFile')

# exit in the case of unknown variables
if unknownVariables:
  print "*** The following build variables are unknown:", unknownVariables.keys()
  Exit(1)


# valid solver for CUDA?
if env['parallelization'] in ['cuda', 'mpi_with_cuda'] and env['solver'] != 'rusanov' and env['solver'] != 'fwave':
  print '** The "'+env['solver']+'" solver is not supported in CUDA.'
  Exit(1)

# CUDA parallelization for openGL
if env['parallelization'] != 'cuda' and env['openGL'] == True:
  print '** The parallelization "'+env['parallelization']+'" does not support OpenGL visualization (CUDA only).'
  Exit(1)

#
# precompiler, compiler and linker flags
#

# Set Intel compiler if selected
# If compiler is set to gnu, we use the default compiler
if env['parallelization'] in ['mpi', 'mpi_with_cuda']:
  # TODO when compilerPath/linkerPath is set, this is overwritten 
  env['CXX'] = 'mpiCC'
  if env['compiler'] == 'intel':
    # Environment variables to switch compiler for different MPI libraries
    envVars = ['OMPI_CXX', 'MPICH_CXX']
    for var in envVars:
      env['ENV'][var] = 'icpc'
else:
  if env['compiler'] == 'intel':
    env['CXX'] = 'icpc'

# eclipse specific flag
env.Append(CCFLAGS=['-fmessage-length=0'])

# xml parameters for the compiler TODO

# set (pre-)compiler flags for the compile modes
if env['compileMode'] == 'debug':
  env.Append(CPPDEFINES=['DEBUG'])

  if env['compiler'] == 'gnu':
    env.Append(CCFLAGS=['-O0','-g3','-Wall'])

  elif env['compiler'] == 'intel':
    env.Append(CCFLAGS=['-O0','-g'])

elif env['compileMode'] == 'release':
  env.Append(CPPDEFINES=['NDEBUG'])

  if env['compiler'] == 'gnu':
    env.Append(CCFLAGS=['-O3','-mtune=native'])

  elif env['compiler'] == 'intel':
    env.Append(CCFLAGS=['-O2'])
    
# Other compiler flags (for all compilers)
env.Append(CCFLAGS=['-fstrict-aliasing', '-fargument-noalias'])

# Vectorization?
if env['compileMode'] == 'release' and env['vectorize']:
  env.Append(CPPDEFINES=['VECTORIZE'])
  if env['compiler'] == 'intel':
    env.Append(CCFLAGS=['-xHost'])
if env['compiler'] == 'intel' and env['showVectorization']:
  env.Append(CCFLAGS=['-vec-report2'])
  
# Platform
if env['compiler'] == 'intel' and env['platform'] == 'mic':
  env.Append(CCFLAGS=['-mmic'])
  env.Append(LINKFLAGS=['-mmic'])
  # Add Intel specific libraries
  env.Append(LIBS=['svml', 'imf', 'intlc'])
  
# Add source directory to include path (important for subdirectories)
env.Append(CPPPATH=['.'])

# set the precompiler variables for the solver
if env['solver'] == 'fwave':
  env.Append(CPPDEFINES=['WAVE_PROPAGATION_SOLVER=1'])
elif env['solver'] == 'augrie':
  env.Append(CPPDEFINES=['WAVE_PROPAGATION_SOLVER=2'])
elif env['solver'] == 'hybrid':
  env.Append(CPPDEFINES=['WAVE_PROPAGATION_SOLVER=0'])
elif env['solver'] == 'fwavevec':
  env.Append(CPPDEFINES=['WAVE_PROPAGATION_SOLVER=4'])

# set the precompiler flags for CUDA
if env['parallelization'] in ['cuda', 'mpi_with_cuda']:
  env.Append(CPPDEFINES=['CUDA'])
  
  # set the directories for the CudaTool
  if 'cudaToolkitDir' in env:
    env['CUDA_TOOLKIT_PATH'] = env['cudaToolkitDir']
    env.Append(RPATH=[os.path.join(env['cudaToolkitDir'], 'lib64')])

  env.Tool('CudaTool', toolpath = ['.'])
  
  # set precompiler flag for nvcc
  env.Append(NVCCFLAGS=' -DCUDA')

  # set the compute capability of the cuda compiler (needs to be set after the CudaTool
  env.Append(NVCCFLAGS=' --gpu-architecture=')
  env.Append(NVCCFLAGS=env['computeCapability'])
  
  # compile explicitly with 64-bit on Mac OS X
  if env['PLATFORM'] == 'darwin':
    env.Append(NVCCFLAGS=' -m64')

# set the precompiler flags for MPI (CUDA)
if env['parallelization'] in ['mpi_with_cuda']:
  env.Append(NVCCFLAGS=' -DUSEMPI')

# set the precompiler flags for MPI (C++)
if env['parallelization'] in ['mpi_with_cuda', 'mpi']:
  env.Append(CPPDEFINES=['USEMPI'])
  if 'compilerPath' in env:
    env['CXX'] = env['compilerPath']
  else:
    env['CXX'] = 'mpiCC'
  if 'linkerPath' in env:
    env['LINKERFORPROGRAMS'] = env['linkerPath']
  env['LINKERFORPROGRAMS'] = 'mpiCC'

if env['openGL'] == True:
  env.Append(LIBS=['SDL', 'GL', 'GLU'])
  if env['openGL_instr'] == True:
    # We assume that SDL_ttf is in the same directory as SDL
    env.Append(LIBS=['SDL_ttf'])
    env.Append(CPPDEFINES=['USESDLTTF'])

# set the compiler flags for libSDL
if 'libSDLDir' in env:
  env.Append(CPPPATH=[env['libSDLDir']+'/include/SDL'])
  env.Append(LIBPATH=[env['libSDLDir']+'/lib'])

# set the precompiler flags and includes for netCDF
if env['writeNetCDF'] == True:
  env.Append(CPPDEFINES=['WRITENETCDF'])
  env.Append(LIBS=['netcdf'])
  # set netCDF location
  if 'netCDFDir' in env:
    env.Append(CPPPATH=[env['netCDFDir']+'/include'])
    env.Append(LIBPATH=[os.path.join(env['netCDFDir'], 'lib')])
    env.Append(RPATH=[os.path.join(env['netCDFDir'], 'lib')])

# set the precompiler flags, includes and libraries for ASAGI
if env['asagi'] == True:
  env.Append(CPPDEFINES=['ASAGI'])
  if env['parallelization'] == 'none' or env['parallelization'] == 'cuda':
    env.Append(CPPDEFINES=['ASAGI_NOMPI'])
    env.Append(LIBS=['asagi_nompi'])
  else:
    env.Append(LIBS=['asagi'])
  env.Append(LIBS=['netcdf_c++4'])
  if 'asagiDir' in env:
    env.Append(CPPPATH=[env['asagiDir']+'/include'])
    env.Append(LIBPATH=[env['asagiDir']+'/lib'])
    env.Append(RPATH=[os.path.join(env['asagiDir'], 'lib')])
  if 'netCDFDir' in env:
    env.Append(LIBPATH=[env['netCDFDir']+'/lib'])
    env.Append(RPATH=[os.path.join(env['netCDFDir'], 'lib')])
  if 'asagiInputDir' in env:
    env.Append(CPPFLAGS=['\'-DASAGI_INPUT_DIR="'+env['asagiInputDir']+'"\''])

# xml runtime parameters
if env['xmlRuntime'] == True: #TODO
  print 'xml runtime parameters are not implemented so far.'
  Exit(1)
  env.Append(CPPDEFINES=['READXML'])
  #set xmllib2 location
  if 'libxmlDir' in env:
    env.Append(CPPPATH=[env['libxmlDir']+'/include/libxml2'])
    env.Append(LIBPATH=[env['libxmlDir']+'/lib'])

#
# setup the program name and the build directory
#
program_name = 'SWE'

# compiler
program_name += '_'+env['compiler']

# compile mode
program_name += '_'+env['compileMode']

# parallelization
program_name += '_'+env['parallelization']

# parallelization
program_name += '_'+env['solver']

# build directory
build_dir = env['buildDir']+'/build_'+program_name

# get the src-code files
env.src_files = []
Export('env')
SConscript('src/SConscript', variant_dir=build_dir, duplicate=0)
Import('env')

# build the program
env.Program('build/'+program_name, env.src_files)
