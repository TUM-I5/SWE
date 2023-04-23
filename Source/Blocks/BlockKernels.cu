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

#include "Block.cuh"
#include "BlockKernels.cuh"

namespace cuda
{
    /**
        Sets corner values of hd_ (only needed for visualization)
        @param hd_		h-values on device
    */
    __global__ void kernelHdBufferEdges(RealType* hd_, int nx, int ny)
    {
        hd_[getCellCoord(0, 0, ny)]           = hd_[getCellCoord(1, 1, ny)];
        hd_[getCellCoord(0, ny + 1, ny)]      = hd_[getCellCoord(1, ny, ny)];
        hd_[getCellCoord(nx + 1, 0, ny)]      = hd_[getCellCoord(nx, 1, ny)];
        hd_[getCellCoord(nx + 1, ny + 1, ny)] = hd_[getCellCoord(nx, ny, ny)];

        // Corresponding C-Code:
        // h[0][0] = h[1][1];
        // h[0][ny+1] = h[1][ny];
        // h[nx+1][0] = h[nx][1];
        // h[nx+1][ny+1] = h[nx][ny];
    }

    //******************************************************************
    // kernels to implement boundary conditions
    //******************************************************************

    /**
     * CUDA kernel to set left boundary layer for conditions BoundaryType::Wall & BoundaryType::Outflow
     * blockIdx.y and threadIdx.y loop over the boundary elements
     * SWE_Block size ny is assumed to be a multiple of the kTileSize
     */
    __global__ void kernelLeftBoundary(
        RealType* hd_, RealType* hud_, RealType* hvd_, int nx, int ny, BoundaryType bound
    )
    {
        int j     = 1 + kTileSize * blockIdx.y + threadIdx.y;
        int ghost = getCellCoord(0, j, ny);
        int inner = getCellCoord(1, j, ny);

        // consider only BoundaryType::Wall & BoundaryType::Outflow boundary conditions
        hd_[ghost]  = hd_[inner];
        hud_[ghost] = (bound == BoundaryType::Wall) ? -hud_[inner] : hud_[inner];
        hvd_[ghost] = hvd_[inner];
    }

    /**
     * CUDA kernel to set right boundary layer for conditions BoundaryType::Wall & BoundaryType::Outflow
     * blockIdx.y and threadIdx.y loop over the boundary elements
     * SWE_Block size ny is assumed to be a multiple of the kTileSize
     */
    __global__ void kernelRightBoundary(
        RealType* hd_, RealType* hud_, RealType* hvd_, int nx, int ny, BoundaryType bound
    )
    {
        int j     = 1 + kTileSize * blockIdx.y + threadIdx.y;
        int ghost = getCellCoord(nx + 1, j, ny);
        int inner = getCellCoord(nx, j, ny);

        // consider only BoundaryType::Wall & BoundaryType::Outflow boundary conditions
        hd_[ghost]  = hd_[inner];
        hud_[ghost] = (bound == BoundaryType::Wall) ? -hud_[inner] : hud_[inner];
        hvd_[ghost] = hvd_[inner];
    }

    /**
     * CUDA kernel to set bottom boundary layer for conditions BoundaryType::Wall & BoundaryType::Outflow
     * blockIdx.x and threadIdx.x loop over the boundary elements
     * SWE_Block size ny is assumed to be a multiple of the kTileSize
     */
    __global__ void kernelBottomBoundary(
        RealType* hd_, RealType* hud_, RealType* hvd_, int nx, int ny, BoundaryType bound
    )
    {
        int i     = 1 + kTileSize * blockIdx.x + threadIdx.x;
        int ghost = getCellCoord(i, 0, ny);
        int inner = getCellCoord(i, 1, ny);

        // consider only BoundaryType::Wall & BoundaryType::Outflow boundary conditions
        hd_[ghost]  = hd_[inner];
        hud_[ghost] = hud_[inner];
        hvd_[ghost] = (bound == BoundaryType::Wall) ? -hvd_[inner] : hvd_[inner];
    }

    /**
     * CUDA kernel to set bottom boundary layer for conditions BoundaryType::Wall & BoundaryType::Outflow
     * blockIdx.x and threadIdx.x loop over the boundary elements
     */
    __global__ void kernelTopBoundary(RealType* hd_, RealType* hud_, RealType* hvd_, int nx, int ny, BoundaryType bound)
    {
        int i     = 1 + kTileSize * blockIdx.x + threadIdx.x;
        int ghost = getCellCoord(i, ny + 1, ny);
        int inner = getCellCoord(i, ny, ny);

        // consider only BoundaryType::Wall & BoundaryType::Outflow boundary conditions
        hd_[ghost]  = hd_[inner];
        hud_[ghost] = hud_[inner];
        hvd_[ghost] = (bound == BoundaryType::Wall) ? -hvd_[inner] : hvd_[inner];
    }

    /**
     * CUDA kernel to set bottom boundary layer according to the external
     * ghost layer status (conditions BoundaryType::Passive and BoundaryType::Connect)
     * blockIdx.x and threadIdx.x loop over the boundary elements.
     * Note that diagonal elements are currently not copied!
     * SWE_Block size ny is assumed to be a multiple of the kTileSize
     */
    __global__ void kernelBottomGhostBoundary(
        RealType* hd_, RealType* hud_, RealType* hvd_, RealType* bottomGhostLayer, int nx, int ny
    )
    {
        int i     = 1 + kTileSize * blockIdx.x + threadIdx.x;
        int ghost = getCellCoord(i, 0, ny);

        hd_[ghost]  = bottomGhostLayer[i];
        hud_[ghost] = bottomGhostLayer[(nx + 2) + i];
        hvd_[ghost] = bottomGhostLayer[2 * (nx + 2) + i];
    }

    /**
     * CUDA kernel to set top boundary layer according to the external
     * ghost layer status (conditions BoundaryType::Passive and BoundaryType::Connect)
     * blockIdx.x and threadIdx.x loop over the boundary elements
     * Note that diagonal elements are currently not copied!
     * SWE_Block size ny is assumed to be a multiple of the kTileSize
     */
    __global__ void kernelTopGhostBoundary(
        RealType* hd_, RealType* hud_, RealType* hvd_, RealType* topGhostLayer, int nx, int ny
    )
    {
        int i     = 1 + kTileSize * blockIdx.x + threadIdx.x;
        int ghost = getCellCoord(i, ny + 1, ny);

        hd_[ghost]  = topGhostLayer[i];
        hud_[ghost] = topGhostLayer[(nx + 2) + i];
        hvd_[ghost] = topGhostLayer[2 * (nx + 2) + i];
    }

    /**
     * CUDA kernel to update bottom copy layer according
     * (for boundary conditions BoundaryType::Passive and BoundaryType::Connect)
     * blockIdx.x and threadIdx.x loop over the boundary elements.
     * Note that diagonal elements are currently not copied!
     * SWE_Block size ny is assumed to be a multiple of the kTileSize
     */
    __global__ void kernelBottomCopyLayer(
        RealType* hd_, RealType* hud_, RealType* hvd_, RealType* bottomCopyLayer, int nx, int ny
    )
    {
        int i    = 1 + kTileSize * blockIdx.x + threadIdx.x;
        int copy = getCellCoord(i, 1, ny);

        bottomCopyLayer[i]                = hd_[copy];
        bottomCopyLayer[(nx + 2) + i]     = hud_[copy];
        bottomCopyLayer[2 * (nx + 2) + i] = hvd_[copy];
    }

    /**
     * CUDA kernel to set top boundary layer according to the external
     * ghost layer status (conditions BoundaryType::Passive and BoundaryType::Connect)
     * blockIdx.x and threadIdx.x loop over the boundary elements
     * Note that diagonal elements are currently not copied!
     * SWE_Block size ny is assumed to be a multiple of the kTileSize
     */
    __global__ void kernelTopCopyLayer(
        RealType* hd_, RealType* hud_, RealType* hvd_, RealType* topCopyLayer, int nx, int ny
    )
    {
        int i    = 1 + kTileSize * blockIdx.x + threadIdx.x;
        int copy = getCellCoord(i, ny, ny);

        topCopyLayer[i]                = hd_[copy];
        topCopyLayer[(nx + 2) + i]     = hud_[copy];
        topCopyLayer[2 * (nx + 2) + i] = hvd_[copy];
    }

} // namespace cuda