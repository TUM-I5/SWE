#!/usr/bin/python

# @file
# This file is part of SWE.
#
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
# Builds the SWE code with all available option combinations.
#

# print the welcome message
print '*********************************************'
print '** Welcome to the build test script of SWE **'
print '*********************************************'
print 'SWE Copyright (C) 2012-2013'
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
print 'This script will build all available build option combinations.'
print 'You may wan to comment out some of the options if you cannot'
print 'build them, e.g. you do have Intel compiler available.'

import collections
import itertools
import subprocess
import sys

#
# possible variables
#
vars = collections.OrderedDict({
    'compiler' : ['gnu', 'intel'],
    'compileMode' : ['debug', 'release'],
    'parallelization' : ['none', 'cuda', 'mpi_with_cuda', 'mpi'],
    'openGL' : ['on', 'off'],
    'openGL_instr' : ['on', 'off'],
    'writeNetCDF' : ['on', 'off'],
    'asagi' : ['on', 'off'],
#    'solver' : ['rusanov', 'fwave', 'augrie', 'hybrid', 'fwavevec'],
    'solver' : ['fwave', 'augrie', 'hybrid', 'fwavevec'],
    'vectorize' : ['on', 'off'],
    'libSDLDir' : ['/work/local'],
})

# Run scons with all combinations
for options in itertools.product(*vars.values()):
    cmd = ['scons']
    for i, var in enumerate(vars):
        cmd.append(var+'='+options[i])
        
    # Run scons
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    _, errmsg = p.communicate()
    
    if p.returncode == 3:
        # These are errors that happen due to unsupported option combinations
        # At the moment we just ignore them.
        pass
    elif p.returncode != 0:
        print ''
        print 'Failed to build SWE with the following options:'
        for i, var in enumerate(vars):
            print '\t'+var+'='+options[i]
        print 'The error message is:'
        print '>>>>>>>>>>>>>>>>>>>>>'
        print errmsg.rstrip()
        print '>>>>>>>>>>>>>>>>>>>>>'
        print 'Run the following command, to build SWE with this options manually:'
        print '\t'+' '.join(cmd)
        
        sys.exit(p.returncode)
        
print ''
print 'Done. Build everything successfully.'
print 'WARNING: Some combinations result in the same executable, so you may not'
print 'find all versions in your build directory!'