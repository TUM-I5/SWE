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
#include "SWE_RusanovBlockCUDA_kernels.hh"

//******************************************************************
// kernels to implement Euler time-stepping
//******************************************************************

inline __device__
float computeFlux(float fLow, float fHigh, float xiLow, float xiHigh, float llf) {
  // local Lax-Friedrich
  return 0.5f*(fLow+fHigh) - 0.5f*llf*(xiHigh-xiLow);
}

/**
 * computes the flux vector components Fhd, Fhud and Fhvd for a single 
 * edge by calling the function computeFlux 
 */
__global__
void kernelComputeFluxesF(float* hd, float* hud, float* hvd,
                          float* Fhd, float* Fhud, float* Fhvd,
                          int ny, float g, float llf, int istart)
{
   int i = istart + TILE_SIZE*blockIdx.x + threadIdx.x;
   int j = 1 + TILE_SIZE*blockIdx.y + threadIdx.y;
   int iL = getCellCoord(i,j,ny);	// index of left cell
   int iR = getCellCoord(i+1,j,ny);	// index of right cell
   int iEdge = getEdgeCoord(i,j,ny);	// index of current edge

   float upwind = max( fabs(hud[iL]/hd[iL]), fabs(hud[iR]/hd[iR]) ); 
   Fhd[iEdge] = computeFlux( hud[iL], hud[iR], hd[iL], hd[iR], upwind );
   Fhud[iEdge] = computeFlux( hud[iL]*hud[iL]/hd[iL] + 0.5f*g*hd[iL]*hd[iL],
            		      hud[iR]*hud[iR]/hd[iR] + 0.5f*g*hd[iR]*hd[iR],
	        	      hud[iL], 
	        	      hud[iR], 
                              llf );
   Fhvd[iEdge] = computeFlux( hud[iL]*hvd[iL]/hd[iL],hud[iR]*hvd[iR]/hd[iR], 
                              hvd[iL], hvd[iR], 
                              llf );
}

/**
 * computes the flux vector components Ghd, Ghud and Ghvd for a single 
 * edge by calling the function computeFlux 
 */
__global__
void kernelComputeFluxesG(float* hd, float* hud, float* hvd,
                          float* Ghd, float* Ghud, float* Ghvd,
                          int ny, float g, float llf, int jstart)
{
   int i = 1 + TILE_SIZE*blockIdx.x + threadIdx.x;
   int j = jstart + TILE_SIZE*blockIdx.y + threadIdx.y;
   int iB = getCellCoord(i,j  ,ny);
   int iT = getCellCoord(i,j+1,ny);
   int iEdge = getEdgeCoord(i,j,ny);

   float upwind = max( fabs(hvd[iB]/hd[iB]), fabs(hvd[iT]/hd[iT]) ); 
   Ghd[iEdge] = computeFlux( hvd[iB], hvd[iT], hd[iB], hd[iT], upwind );
   Ghud[iEdge] = computeFlux( hud[iB]*hvd[iB]/hd[iB],hud[iT]*hvd[iT]/hd[iT], 
                              hud[iB], hud[iT], 
                              llf );
   Ghvd[iEdge] = computeFlux( hvd[iB]*hvd[iB]/hd[iB] + 0.5f*g*hd[iB]*hd[iB],
                              hvd[iT]*hvd[iT]/hd[iT] + 0.5f*g*hd[iT]*hd[iT],
			      hvd[iB], hvd[iT], 
                              llf );
}

/**
 * computes the bathymetry source terms for the hu and hv equation for 
 * a given cell in the resp. array elements Bxd and Byd
 */
__global__
void kernelComputeBathymetrySources(float* hd, float* bd, float* Bxd, float* Byd, 
                                    int ny, float g)
{
// Note: different index ranges for h and b vs. Bxd, Byd: 
//       [0..nx+]x[0..ny+1] vs. [1..nx]x[1..ny]
// Note: indices for Bxd, Byd shifted to start with 0
   int i = TILE_SIZE*blockIdx.x + threadIdx.x;
   int j = TILE_SIZE*blockIdx.y + threadIdx.y;
   
   // compute indices of involved array elements
   int ij = getBathyCoord(i,j,ny);
   int left  = getCellCoord(i  ,j+1,ny); // index of left cell (arrays hd,bd)
   int right = getCellCoord(i+2,j+1,ny); // index of right cell (array hd,bb)

   Bxd[ij] = g * 0.5f*(hd[right] + hd[left]) * 0.5f*(bd[right] - bd[left]);

   int bot = getCellCoord(i+1,j,ny);   // index of left cell (arrays hd,bd)
   int top = getCellCoord(i+1,j+2,ny); // index of right cell (array hd,bb)

   Byd[ij] = g * 0.5f*(hd[top] + hd[bot]) * 0.5f*(bd[top] - bd[bot]);

}


/**
 * CUDA kernel for Euler time step
 */
__global__
void kernelEulerTimestep(float* hd, float* hud, float* hvd,
                         float* Fhd, float* Fhud, float* Fhvd,
                         float* Ghd, float* Ghud, float* Ghvd,
			 float* Bxd, float* Byd,
			 float* maxhd, float* maxvd,
                         int nx, int ny, float dt, float dxi, float dyi)
{

   __shared__ float Fds[TILE_SIZE+1][TILE_SIZE+1];
   __shared__ float Gds[TILE_SIZE+1][TILE_SIZE+1];

   int tx = threadIdx.x;
   int ty = threadIdx.y;
   
   int i = 1 + TILE_SIZE*blockIdx.x + tx;
   int j = 1 + TILE_SIZE*blockIdx.y + ty;
   int iElem = getCellCoord(i,j,ny);   // index of current cell
   int iEdge = getEdgeCoord(i,j,ny);   // index of right/top Edge
   int iLeft = getEdgeCoord(i-1,j,ny); // index of left Edge
   int iBot  = getEdgeCoord(i,j-1,ny); // index of bottom Edge
   
   float h;
   float hu;
   float hv;

   // copy flux unknowns from global into local memory
   // -> for fluxes corresponding to variable h
   Fds[tx+1][ty] = Fhd[iEdge];
   Gds[tx][ty+1] = Ghd[iEdge];
   if (tx==0) Fds[tx][ty] = Fhd[iLeft];
   if (ty==0) Gds[tx][ty] = Ghd[iBot];
   __syncthreads();

   // compute new value of h from fluxes
   h = hd[iElem] - dt *( (Fds[tx+1][ty]-Fds[tx][ty])*dxi 
           	        +(Gds[tx][ty+1]-Gds[tx][ty])*dyi );
   __syncthreads();

   // copy flux unknowns from global into local memory
   // -> for fluxes corresponding to variable hu
   Fds[tx+1][ty] = Fhud[iEdge];
   Gds[tx][ty+1] = Ghud[iEdge];
   if (tx==0) Fds[tx][ty] = Fhud[iLeft];
   if (ty==0) Gds[tx][ty] = Ghud[iBot];
   __syncthreads();

   // compute new value of hu from fluxes
   hu = hud[iElem] - dt *( (Fds[tx+1][ty]-Fds[tx][ty])*dxi 
           	          +(Gds[tx][ty+1]-Gds[tx][ty])*dyi 
			  + Bxd[getBathyCoord(i-1,j-1,ny)]*dxi );
   __syncthreads();

   // copy flux unknowns from global into local memory
   // -> for fluxes corresponding to variable hv
   Fds[tx+1][ty] = Fhvd[iEdge];
   Gds[tx][ty+1] = Ghvd[iEdge];
   if (tx==0) Fds[tx][ty] = Fhvd[iLeft];
   if (ty==0) Gds[tx][ty] = Ghvd[iBot];
   __syncthreads();

   // compute new value of hv from fluxes
   hv = hvd[iElem] - dt *( (Fds[tx+1][ty]-Fds[tx][ty])*dxi 
           	          +(Gds[tx][ty+1]-Gds[tx][ty])*dyi 
			  + Byd[getBathyCoord(i-1,j-1,ny)]*dyi );
   __syncthreads();

   /* precompute maxmimal height and velocity per thread block 
    * (for computation of allowed time step size)
    */

   // compute absolute values of h and absolute velocity
   hd[iElem] = h; Fds[tx][ty] = h;
   hud[iElem] = hu; hu = (h>0.0) ? fabs(hu/h) : 0.0;
   hvd[iElem] = hv; hv = (h>0.0) ? fabs(hv/h) : 0.0; 
   Gds[tx][ty] = (hu>hv) ? hu : hv;
   
   // parallel reduction on thread block:
   // determine maximum wave height and velocity
   // step 1: reduction in ty-direction
   for (i=TILE_SIZE>>1; i>0; i>>=1) {
      __syncthreads();
      if (ty < i) {
         if( Fds[tx][ty] < Fds[tx][ty+i]) Fds[tx][ty] = Fds[tx][ty+i];
         if( Gds[tx][ty] < Gds[tx][ty+i]) Gds[tx][ty] = Gds[tx][ty+i];
      };
   };
   // step 2: reduction in ty-direction
   for (i=TILE_SIZE>>1; i>0; i>>=1) {
      __syncthreads();
      if ((tx < i) && (ty==0)) {
         if( Fds[tx][ty] < Fds[tx+i][ty]) Fds[tx][ty] = Fds[tx+i][ty];
         if( Gds[tx][ty] < Gds[tx+i][ty]) Gds[tx][ty] = Gds[tx+i][ty];
      };
   };
   // save maxima in array maxhd and maxvd
   
   if ((tx == 0) && (ty==0)) {
      j = blockIdx.x*(nx/TILE_SIZE)+blockIdx.y;
      maxhd[j] = Fds[0][0];
      maxvd[j] = Gds[0][0];
   };
   
}


//******************************************************************
// kernels to implement boundary conditions
//******************************************************************


/**
 * CUDA kernel for maximum reduction
 * required to compute maximum water height and velocities to determine 
 * allow time step
 */
__global__ 
void kernelMaximum(float* maxhd, float* maxvd, int start, int size) {
  int tx = start+threadIdx.x;
  for (int i=size>>1; i>0; i>>=1) {
     __syncthreads();
     if (tx < i) {
        if( maxhd[tx] < maxhd[tx+i] ) maxhd[tx] = maxhd[tx+i];
        if( maxvd[tx] < maxvd[tx+i] ) maxvd[tx] = maxvd[tx+i];
     };
  };
}



