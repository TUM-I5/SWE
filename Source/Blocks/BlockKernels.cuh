/**
 * @file
 * This file is part of SWE.
 *
 * @author Michael Bader (bader AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Univ.-Prof._Dr._Michael_Bader)
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

#ifndef __SWE_BLOCKCUDAKERNELS_HH
#define __SWE_BLOCKCUDAKERNELS_HH

namespace cuda
{
    // declaration of CUDA kernels
    __global__ void kernelHdBufferEdges(RealType* hd_, int nx, int ny);

    __global__ void kernelMaximum(RealType* maxhd, RealType* maxvd, int start, int size);

    __global__ void kernelLeftBoundary(
        RealType* hd_, RealType* hud_, RealType* hvd_, int nx, int ny, BoundaryType bound
    );
    __global__ void kernelRightBoundary(
        RealType* hd_, RealType* hud_, RealType* hvd_, int nx, int ny, BoundaryType bound
    );
    __global__ void kernelBottomBoundary(
        RealType* hd_, RealType* hud_, RealType* hvd_, int nx, int ny, BoundaryType bound
    );
    __global__ void kernelTopBoundary(
        RealType* hd_, RealType* hud_, RealType* hvd_, int nx, int ny, BoundaryType bound
    );
    __global__ void kernelBottomGhostBoundary(
        RealType* hd_, RealType* hud_, RealType* hvd_, RealType* bottomGhostLayer, int nx, int ny
    );
    __global__ void kernelTopGhostBoundary(
        RealType* hd_, RealType* hud_, RealType* hvd_, RealType* topGhostLayer, int nx, int ny
    );
    __global__ void kernelBottomCopyLayer(
        RealType* hd_, RealType* hud_, RealType* hvd_, RealType* bottomCopyLayer, int nx, int ny
    );
    __global__ void kernelTopCopyLayer(
        RealType* hd_, RealType* hud_, RealType* hvd_, RealType* topCopyLayer, int nx, int ny
    );
} // namespace cuda

#endif
