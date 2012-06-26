#! /usr/bin/python

# @file
# This file is part of SWE.
#
# @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
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
print 'This program comes with ABSOLUTELY NO WARRANTY.'
print 'This is free software, and you are welcome to redistribute it'
print 'under certain conditions.'
print 'Details can be found in the file \'gpl.txt\'.'
print ''

#
# set possible variables
#
vars = Variables()

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

  EnumVariable( 'computeCapability', 'optional architecture/compute capability of the CUDA card', 'sm_21',
                allowed_values=('sm_20', 'sm_21', 'sm_22', 'sm_23')
              ),

  BoolVariable( 'writeNetCDF', 'write output in the netCDF-format', False ),

  BoolVariable( 'asagi', 'use ASAGI', False ),

  EnumVariable( 'solver', 'Riemann solver', 'augrie',
                allowed_values=('rusanov', 'fwave', 'augrie', 'hybrid')
              )
)

# external variables
vars.AddVariables(
  PathVariable( 'cudaToolkitDir', 'location of the CUDA toolkit', None ),
  PathVariable( 'cudaSDKDir', 'location of the CUDA SDK', None),
  PathVariable( 'netCDFDir', 'location of netCDF', None),
  PathVariable( 'asagiDir', 'location of ASAGI', None)
)

# set environment
env = Environment(variables=vars)

# generate help text
Help(vars.GenerateHelpText(env))

# valid solver for CUDA?
if env['parallelization'] == 'cuda' and env['solver'] != 'rusanov' and env['solver'] != 'fwave':
  print '** The "'+env['solver']+'" solver is not supported in CUDA.'
  Exit(1)

#
# precompiler, compiler and linker flags
#

# eclipse specific flag
env.Append(CCFLAGS=['-fmessage-length=0'])

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

# set the precompiler variables for the solver
if env['solver'] == 'fwave':
  env.Append(CPPDEFINES=['WAVE_PROPAGATION_SOLVER=1'])
elif env['solver'] == 'augrie':
  env.Append(CPPDEFINES=['WAVE_PROPAGATION_SOLVER=2'])
elif env['solver'] == 'hybrid':
  env.Append(CPPDEFINES=['WAVE_PROPAGATION_SOLVER=0'])

# set the precompiler flags for CUDA
if env['parallelization'] == 'cuda':
  env.Append(CPPDEFINES=['CUDA'])
  env.Append(CPPDEFINES=['NOMPI'])
  
  # set the directories for the CudaTool
  if 'cudaToolkitDir' in env:
    env['CUDA_TOOLKIT_PATH'] = env['cudaToolkitDir']
  if 'cudaSDKDir' in env:
    env['CUDA_SDK_PATH'] = env['cudaSDKDir']

  env.Tool('CudaTool', toolpath = ['.'])

  # set the compute capability of the cuda compiler (needs to be set after the CudaTool
  env.Append(NVCCFLAGS=' --gpu-architecture=')
  env.Append(NVCCFLAGS=env['computeCapability'])

# set the precompiler flags and includes for netCDF
if int(env['writeNetCDF']):
  env.Append(CPPDEFINES=['WRITENETCDF'])
  # set netCDF location
  if 'netCDFDir' in env:
    env.Append(CPPPATH=[env['netCDFDir']+'/include'])

# set the precompiler flags, includes and libraries for ASAGI
if int(env['asagi']):
  env.Append(CPPDEFINES=['ASAGI'])
  if env['parallelization'] == 'none' or env['parallelization'] == 'cuda':
    env.Append(CPPDEFINES=['ASAGI_NOMPI'])
  env.Append(LIBS=['netcdf_c++4'])
  env.Append(LIBS=['asagi'])
  if 'asagiDir' in env:
    env.Append(CPPPATH=[env['asagiDir']+'/include'])
    env.Append(LIBPATH=[env['asagiDir']])
  if 'netCDFDir' in env:
    env.Append(LIBPATH=[env['netCDFDir']+'/lib'])


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
