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

#include "SWE_BlockCUDA.hh"
#include "SWE_BlockCUDA_kernels.hh"

/**
    Sets corner values of hd (only needed for visualization)
	@param hd		h-values on device
*/	
__global__
void kernelHdBufferEdges(float* hd, int nx, int ny)
{ 
  hd[getCellCoord(0   ,0   ,ny)] = hd[getCellCoord(1 ,1 ,ny)];
  hd[getCellCoord(0   ,ny+1,ny)] = hd[getCellCoord(1 ,ny,ny)];
  hd[getCellCoord(nx+1,0   ,ny)] = hd[getCellCoord(nx,1 ,ny)];
  hd[getCellCoord(nx+1,ny+1,ny)] = hd[getCellCoord(nx,ny,ny)];
	
  //Corresponding C-Code:
  //h[0][0] = h[1][1];
  //h[0][ny+1] = h[1][ny];
  //h[nx+1][0] = h[nx][1];
  //h[nx+1][ny+1] = h[nx][ny];
}


//******************************************************************
// kernels to implement boundary conditions
//******************************************************************

/**
 * CUDA kernel to set left boundary layer for conditions WALL & OUTFLOW
 * blockIdx.y and threadIdx.y loop over the boundary elements
 * SWE_Block size ny is assumed to be a multiple of the TILE_SIZE
 */
__global__
void kernelLeftBoundary(float* hd, float* hud, float* hvd,
                        int nx, int ny, BoundaryType bound)
{
  int j = 1 + TILE_SIZE*blockIdx.y + threadIdx.y;
  int ghost = getCellCoord(0,j,ny);
  int inner = getCellCoord(1,j,ny);
  
  // consider only WALL & OUTFLOW boundary conditions
  hd[ghost] = hd[inner];
  hud[ghost] = (bound==WALL) ? -hud[inner] : hud[inner];
  hvd[ghost] = hvd[inner];
}

/**
 * CUDA kernel to set right boundary layer for conditions WALL & OUTFLOW
 * blockIdx.y and threadIdx.y loop over the boundary elements
 * SWE_Block size ny is assumed to be a multiple of the TILE_SIZE
 */
__global__
void kernelRightBoundary(float* hd, float* hud, float* hvd,
                         int nx, int ny, BoundaryType bound)
{
  int j = 1 + TILE_SIZE*blockIdx.y + threadIdx.y;
  int ghost = getCellCoord(nx+1,j,ny);
  int inner = getCellCoord(nx  ,j,ny);
  
  // consider only WALL & OUTFLOW boundary conditions
  hd[ghost] = hd[inner];
  hud[ghost] = (bound==WALL) ? -hud[inner] : hud[inner];
  hvd[ghost] = hvd[inner];
}


/**
 * CUDA kernel to set bottom boundary layer for conditions WALL & OUTFLOW
 * blockIdx.x and threadIdx.x loop over the boundary elements
 * SWE_Block size ny is assumed to be a multiple of the TILE_SIZE
 */
__global__
void kernelBottomBoundary(float* hd, float* hud, float* hvd,
                          int nx, int ny, BoundaryType bound)
{
  int i = 1 + TILE_SIZE*blockIdx.x + threadIdx.x;
  int ghost = getCellCoord(i,0,ny);
  int inner = getCellCoord(i,1,ny);
  
  // consider only WALL & OUTFLOW boundary conditions
  hd[ghost] = hd[inner];
  hud[ghost] = hud[inner];
  hvd[ghost] = (bound==WALL) ? -hvd[inner] : hvd[inner]; 
}

/**
 * CUDA kernel to set bottom boundary layer for conditions WALL & OUTFLOW
 * blockIdx.x and threadIdx.x loop over the boundary elements
 */
__global__
void kernelTopBoundary(float* hd, float* hud, float* hvd,
                       int nx, int ny, BoundaryType bound)
{
  int i = 1 + TILE_SIZE*blockIdx.x + threadIdx.x;
  int ghost = getCellCoord(i,ny+1,ny);
  int inner = getCellCoord(i,ny  ,ny);
  
  // consider only WALL & OUTFLOW boundary conditions
  hd[ghost] = hd[inner];
  hud[ghost] = hud[inner];
  hvd[ghost] = (bound==WALL) ? -hvd[inner] : hvd[inner]; 
}

/**
 * CUDA kernel to set bottom boundary layer according to the external 
 * ghost layer status (conditions PASSIVE and CONNECT)
 * blockIdx.x and threadIdx.x loop over the boundary elements.
 * Note that diagonal elements are currently not copied!
 * SWE_Block size ny is assumed to be a multiple of the TILE_SIZE
 */
__global__
void kernelBottomGhostBoundary(float* hd, float* hud, float* hvd,
                               float* bottomGhostLayer, int nx, int ny)
{
  int i = 1 + TILE_SIZE*blockIdx.x + threadIdx.x;
  int ghost = getCellCoord(i,0,ny);

  hd[ghost]  = bottomGhostLayer[i];
  hud[ghost] = bottomGhostLayer[(nx+2)+i];
  hvd[ghost] = bottomGhostLayer[2*(nx+2)+i];
}

/**
 * CUDA kernel to set top boundary layer according to the external 
 * ghost layer status (conditions PASSIVE and CONNECT)
 * blockIdx.x and threadIdx.x loop over the boundary elements
 * Note that diagonal elements are currently not copied!
 * SWE_Block size ny is assumed to be a multiple of the TILE_SIZE
 */
__global__
void kernelTopGhostBoundary(float* hd, float* hud, float* hvd,
                            float* topGhostLayer, int nx, int ny)
{
  int i = 1 + TILE_SIZE*blockIdx.x + threadIdx.x;
  int ghost = getCellCoord(i,ny+1,ny);
  
  hd[ghost]  = topGhostLayer[i];
  hud[ghost] = topGhostLayer[(nx+2)+i];
  hvd[ghost] = topGhostLayer[2*(nx+2)+i];
}

/**
 * CUDA kernel to update bottom copy layer according 
 * (for boundary conditions PASSIVE and CONNECT)
 * blockIdx.x and threadIdx.x loop over the boundary elements.
 * Note that diagonal elements are currently not copied!
 * SWE_Block size ny is assumed to be a multiple of the TILE_SIZE
 */
__global__
void kernelBottomCopyLayer(float* hd, float* hud, float* hvd,
                           float* bottomCopyLayer, int nx, int ny)
{
  int i = 1 + TILE_SIZE*blockIdx.x + threadIdx.x;
  int copy = getCellCoord(i,1,ny);

  bottomCopyLayer[i]          = hd[copy];  
  bottomCopyLayer[(nx+2)+i]   = hud[copy]; 
  bottomCopyLayer[2*(nx+2)+i] = hvd[copy]; 
}

/**
 * CUDA kernel to set top boundary layer according to the external 
 * ghost layer status (conditions PASSIVE and CONNECT)
 * blockIdx.x and threadIdx.x loop over the boundary elements
 * Note that diagonal elements are currently not copied!
 * SWE_Block size ny is assumed to be a multiple of the TILE_SIZE
 */
__global__
void kernelTopCopyLayer(float* hd, float* hud, float* hvd,
                        float* topCopyLayer, int nx, int ny)
{
  int i = 1 + TILE_SIZE*blockIdx.x + threadIdx.x;
  int copy = getCellCoord(i,ny,ny);
  
  topCopyLayer[i]          = hd[copy];  
  topCopyLayer[(nx+2)+i]   = hud[copy]; 
  topCopyLayer[2*(nx+2)+i] = hvd[copy]; 
}


