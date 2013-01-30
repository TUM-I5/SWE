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
# Example build parameters for OpenGL visualizations with ASAGI
#

# Build options
parallelization='cuda'
solver='fwave'
asagi='yes'
openGL='yes'
openGL_instr='yes'

# Hardware settings
computeCapability='sm_21'

# Directory containing ASAGI input files (needs to be set)
asagiInputDir=''

# Library paths (only required of not installed in default path)
#asagiDir=''
#netcdfDir=''
#libSDLDir=''
#cudaToolkitDir=''
