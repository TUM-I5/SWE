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
# Example build parameters for the NPS GPU-cluster.
#

#build options
parallelization='mpi_with_cuda'
computeCapability='sm_20'
solver='fwave'
asagi='yes'
writeNetCDF='yes'

# libraries (machine dependent)
compilerPath='/usr/lib64/openmpi/bin/mpiCC'
linkerPath='/usr/lib64/openmpi/bin/mpiCC'
cudaSDKDir='/'
asagiDir='/home/prof2/software/asagi/mpi'
netCDFDir='/home/prof2/software/netcdf'
