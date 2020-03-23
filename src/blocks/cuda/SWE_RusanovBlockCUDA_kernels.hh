/**
 * @file
 * This file is part of SWE.
 *
 * @author Michael Bader, Kaveh Rahnema, Tobias Schnabel
 *
 * @section LICENSE
 *
 * SWE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SWE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SWE.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * @section DESCRIPTION
 *
 * TODO
 */

#ifndef __SWE_RUSANOVBLOCKCUDAKERNELS_HH
#define __SWE_RUSANOVBLOCKCUDAKERNELS_HH

//******************************************************************
// kernels to implement Euler time-stepping
//******************************************************************

__global__
void kernelComputeFluxesF(float* hd, float* hud, float* hvd,
                          float* Fhd, float* Fhud, float* Fhvd,
                          int ny, float g, float llf, int istart);

__global__
void kernelComputeFluxesG(float* hd, float* hud, float* hvd,
                          float* Ghd, float* Ghud, float* Ghvd,
                          int ny, float g, float llf, int jstart);

__global__
void kernelComputeBathymetrySources(float* hd, float* bd, float* Bxd, float* Byd, 
                                    int ny, float g);

__global__
void kernelEulerTimestep(float* hd, float* hud, float* hvd,
                         float* Fhd, float* Fhud, float* Fhvd,
                         float* Ghd, float* Ghud, float* Ghvd,
			 float* Bxd, float* Byd,
			 float* maxhd, float* maxvd,
                         int nx, int ny, float dt, float dxi, float dyi);


__global__ 
void kernelMaximum(float* maxhd, float* maxvd, int start, int size);

#endif



