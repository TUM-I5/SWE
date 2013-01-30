/**
 * @file
 * This file is part of SWE.
 *
 * @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
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
 * SWE_Block in CUDA, which uses solvers in the wave propagation formulation.
 */

#include "SWE_WavePropagationBlockCuda.hh"
#include "SWE_BlockCUDA.hh"

#include "SWE_WavePropagationBlockCuda_kernels.hh"

#include "tools/Logger.hh"

#include <cassert>

// CUDA-C includes
#include <cuda.h>
#include <cuda_runtime_api.h>

// Thrust library (used for the final maximum reduction in the method computeNumericalFluxes(...))
#include <thrust/device_vector.h>


/**
 * Constructor of a SWE_WavePropagationBlockCuda.
 *
 * Allocates the variables for the simulation:
 *   Please note: The definition of indices changed in contrast to the CPU-Implementation.
 *
 *   unknowns hd,hud,hvd,bd stored on the CUDA device are defined for grid indices [0,..,nx+1]*[0,..,ny+1] (-> Abstract class SWE_BlockCUDA)
 *     -> computational domain is [1,..,nx]*[1,..,ny]
 *     -> plus ghost cell layer
 *
 *   net-updates are defined for edges with indices [0,..,nx]*[0,..,ny] for horizontal and vertical edges for simplicity (one layer is not necessary).
 *
 *   A left/right net update with index (i-1,j) is located on the edge between
 *   cells with index (i-1,j) and (i,j):
 * <pre>
 *   *********************
 *   *         *         *
 *   * (i-1,j) *  (i,j)  *
 *   *         *         *
 *   *********************
 *
 *             *
 *            ***
 *           *****
 *             *
 *             *
 *   NetUpdatesLeft(i-1,j)
 *             or
 *   NetUpdatesRight(i-1,j)
 * </pre>
 *
 *   A below/above net update with index (i, j-1) is located on the edge between
 *   cells with index (i, j-1) and (i,j):
 * <pre>
 *   ***********
 *   *         *
 *   * (i, j)  *   *
 *   *         *  **  NetUpdatesBelow(i,j-1)
 *   *********** *****         or
 *   *         *  **  NetUpdatesAbove(i,j-1)
 *   * (i,j-1) *   *
 *   *         *
 *   ***********
 * </pre>
 * @param i_offsetX spatial offset of the block in x-direction.
 * @param i_offsetY spatial offset of the offset in y-direction.
 * @param i_cudaDevice ID of the CUDA-device, which should be used.
 */
SWE_WavePropagationBlockCuda::SWE_WavePropagationBlockCuda(
		int l_nx, int l_ny,
		float l_dx, float l_dy)
	: SWE_BlockCUDA(l_nx, l_ny, l_dx, l_dy)
{
  // compute the size of one 1D net-update array.
  int sizeOfNetUpdates = (nx+1)*(ny+1)*sizeof(float);

  // allocate CUDA memory for the net-updates
  cudaMalloc((void**)&hNetUpdatesLeftD, sizeOfNetUpdates);
  checkCUDAError("allocate device memory for hNetUpdatesLeftD");

  cudaMalloc((void**)&hNetUpdatesRightD, sizeOfNetUpdates);
  checkCUDAError("allocate device memory for hNetUpdatesRightD");

  cudaMalloc((void**)&huNetUpdatesLeftD, sizeOfNetUpdates);
  checkCUDAError("allocate device memory for huNetUpdatesLeftD");

  cudaMalloc((void**)&huNetUpdatesRightD, sizeOfNetUpdates);
  checkCUDAError("allocate device memory for huNetUpdatesRightD");

  cudaMalloc((void**)&hNetUpdatesBelowD, sizeOfNetUpdates);
  checkCUDAError("allocate device memory for hNetUpdatesBelowD");

  cudaMalloc((void**)&hNetUpdatesAboveD, sizeOfNetUpdates);
  checkCUDAError("allocate device memory for hNetUpdatesAboveD");

  cudaMalloc((void**)&hvNetUpdatesBelowD, sizeOfNetUpdates);
  checkCUDAError("allocate device memory for hvNetUpdatesBelowD");

  cudaMalloc((void**)&hvNetUpdatesAboveD, sizeOfNetUpdates);
  checkCUDAError("allocate device memory for hNetUpdatesAboveD");
}

/**
 * Destructor of a SWE_WavePropagationBlockCuda.
 *
 * Frees all of the memory, which was allocated within the constructor.
 * Resets the CUDA device: Useful if error occured and printf is used on the device (buffer).
 */
SWE_WavePropagationBlockCuda::~SWE_WavePropagationBlockCuda() {
  // free the net-updates memory
  cudaFree(hNetUpdatesLeftD);
  cudaFree(hNetUpdatesRightD);
  cudaFree(huNetUpdatesLeftD);
  cudaFree(huNetUpdatesRightD);

  cudaFree(hNetUpdatesBelowD);
  cudaFree(hNetUpdatesAboveD);
  cudaFree(hvNetUpdatesBelowD);
  cudaFree(hvNetUpdatesAboveD);
}

/**
 * Compute a single global time step of a given time step width.
 * Remark: The user has to take care about the time step width. No additional check is done. The time step width typically available
 *         after the computation of the numerical fluxes (hidden in this method).
 *
 * First the net-updates are computed.
 * Then the cells are updated with the net-updates and the given time step width.
 *
 * @param i_dT time step width in seconds.
 */
__host__
void SWE_WavePropagationBlockCuda::simulateTimestep(float i_dT) {
 // Compute the numerical fluxes/net-updates in the wave propagation formulation.
 computeNumericalFluxes();

 // Update the unknowns with the net-updates.
 updateUnknowns(i_dT);
}

/**
 * perform forward-Euler time steps, starting with simulation time tStart,:
 * until simulation time tEnd is reached; 
 * device-global variables hd, hud, hvd are updated;
 * unknowns h, hu, hv in main memory are not updated.
 * Ghost layers and bathymetry sources are updated between timesteps.
 * intended as main simulation loop between two checkpoints
 */
__host__
float SWE_WavePropagationBlockCuda::simulate(float tStart, float tEnd) {
  float t = tStart;
  do {
     // set values in ghost cells:
     setGhostLayer();
     
     // Compute the numerical fluxes/net-updates in the wave propagation formulation.
     computeNumericalFluxes();

     // Update the unknowns with the net-updates.
     updateUnknowns(maxTimestep);
	 
	 t += maxTimestep;
//     cout << "Simulation at time " << t << endl << flush;
  } while(t < tEnd);

  return t;
}


/**
 * Compute the numerical fluxes (net-update formulation here) on all edges.
 *
 * The maximum wave speed is computed within the net-updates kernel for each CUDA-block.
 * To finalize the method the Thrust-library is called, which does the reduction over all blockwise maxima.
 *   In the wave speed reduction step the actual cell width in x- and y-direction is not taken into account.
 *
 *   TODO: A splitting or direct computation of the time step width might increase the total time step size.
 *     Example:
 *       dx = 11, dy = 6;
 *       max wave speed in x-direction: 10
 *       max wave speed in y-direction: 5.5
 *       max wave speed in both direction: 10
 *
 *       => maximum time step (current implementation): min(11/10, 6/10) = 0.6
 *       => maximum time step (splitting the dimensions): min(11/10, 6/5.5) = 1.09..
 */
void SWE_WavePropagationBlockCuda::computeNumericalFluxes() {
  /*
   * Initialization.
   */

  /**
   * <b>Row-major vs column-major</b>
   *
   * C/C++ arrays are row-major whereas warps are constructed in column-major order from threads/blocks. To get coalesced
   * memory access in CUDA, we can use a 2-dimensional CUDA structure but we have to switch x and y inside a block.
   *
   * This means, we have to switch threadIdx.x <-> threadIdx.y as well as blockDim.x <-> blockDim.y.
   * Important: blockDim has to be switched for the kernel call as well!
   */

  //! definition of one CUDA-block. Typical size are 8*8 or 16*16
  dim3 dimBlock(TILE_SIZE, TILE_SIZE);

  /**
   * Definition of the "main" CUDA-grid.
   * This grid covers only edges 0..#(edges in x-direction)-2 and 0..#(edges in y-direction)-2.
   *
   * An example with a computational domain of size
   *  nx = 24, ny = 16
   * with a 1 cell ghost layer would result in a grid with
   *  (nx+2)*(ny+2) = (26*18)
   * cells and
   *  (nx+1)*(ny+1) = (25*17)
   * edges.
   *
   * The CUDA-blocks (here 8*8) mentioned above would cover all edges except
   * the ones lying between the computational domain and the right/top ghost layer:
   * <pre>
   *                                                          *
   *                                                         **        top ghost layer,
   *                                                        ********   cell ids
   *                        *******************************  **        = (*, ny+1)
   *                        *         *         *         *   *
   *                        *   8*8   *   8*8   *   8*8   *
   *                        *  block  *  block  *  block  *
   *                        *         *         *         *
   *                        *******************************
   *                        *         *         *         *
   *                        *   8*8   *   8*8   *   8*8   *
   *                    *   *  block  *  block  *  block  *
   *     bottom         **  *         *         *         *
   *     ghost     ******** *******************************
   *     layer,         **
   *     cell ids       *   *                              *
   *     =(*,0)            ***                            ***
   *                        *                              *
   *                        *                              *
   *                  left ghost layer,             right ghost layer,
   *                  cell ids = (0,*)             cell ids = (nx+1, *)
   * </pre>
   */
  dim3 dimGrid(nx/TILE_SIZE,ny/TILE_SIZE);

  // assert a valid tile size
  assert(nx%TILE_SIZE==0);
  assert(ny%TILE_SIZE==0);

  // "2D array" which holds the blockwise maximum wave speeds
  float* l_maximumWaveSpeedsD;

  // size of the maximum wave speed array (dimension of the grid + ghost layers, without the top right block), sizeof(float) not included
  int l_sizeMaxWaveSpeeds = ((dimGrid.x+1)*(dimGrid.y+1)-1);
  cudaMalloc((void**)&l_maximumWaveSpeedsD, (l_sizeMaxWaveSpeeds*sizeof(float)) );


  /*
   * Compute the net updates for the 'main part and the two 'boundary' parts.
   */
  // compute the net-updates for the "main" part.
  computeNetUpdatesKernel<<<dimGrid,dimBlock>>>( hd, hud, hvd, bd,
                                                 hNetUpdatesLeftD,  hNetUpdatesRightD,
                                                 huNetUpdatesLeftD, huNetUpdatesRightD,
                                                 hNetUpdatesBelowD, hNetUpdatesAboveD,
                                                 hvNetUpdatesBelowD, hvNetUpdatesAboveD,
                                                 l_maximumWaveSpeedsD,
                                                 nx,ny
                                               );

  // compute the "remaining" net updates (edges "simulation domain"/"top ghost layer" and "simulation domain"/"right ghost layer")
  // edges between cell nx and ghost layer nx+1
  dim3 dimRightBlock(TILE_SIZE, 1);
  dim3 dimRightGrid(1, ny/TILE_SIZE);
  computeNetUpdatesKernel<<<dimRightGrid, dimRightBlock>>>( hd, hud, hvd, bd,
                                                            hNetUpdatesLeftD,  hNetUpdatesRightD,
                                                            huNetUpdatesLeftD, huNetUpdatesRightD,
                                                            hNetUpdatesBelowD, hNetUpdatesAboveD,
                                                            hvNetUpdatesBelowD, hvNetUpdatesAboveD,
                                                            l_maximumWaveSpeedsD,
                                                            nx, ny,
                                                            nx, 0,
                                                            dimGrid.x, 0);

  // edges between cell ny and ghost layer ny+1
  dim3 dimTopBlock(1, TILE_SIZE);
  dim3 dimTopGrid(nx/TILE_SIZE, 1);
  computeNetUpdatesKernel<<<dimTopGrid, dimTopBlock>>>( hd, hud, hvd, bd,
                                                        hNetUpdatesLeftD,  hNetUpdatesRightD,
                                                        huNetUpdatesLeftD, huNetUpdatesRightD,
                                                        hNetUpdatesBelowD, hNetUpdatesAboveD,
                                                        hvNetUpdatesBelowD, hvNetUpdatesAboveD,
                                                        l_maximumWaveSpeedsD,
                                                        nx, ny,
                                                        0, ny,
                                                        0, dimGrid.y);

  /*
   * Finalize (max reduction of the maximumWaveSpeeds-array.)
   *
   * The Thrust library is used in this step.
   * An optional kernel could be written for the maximum reduction.
   */
  // Thrust pointer to the device array
  thrust::device_ptr<float> l_thrustDevicePointer(l_maximumWaveSpeedsD);

  // use Thrusts max_element-function for the maximum reduction
  thrust::device_ptr<float> l_thrustDevicePointerMax = thrust::max_element(l_thrustDevicePointer, l_thrustDevicePointer+l_sizeMaxWaveSpeeds);

  // get the result from the device
  float l_maximumWaveSpeed = l_thrustDevicePointerMax[0];

  // free the max wave speeds array on the device
  cudaFree(l_maximumWaveSpeedsD);

  // set the maximum time step for this SWE_WavePropagationBlockCuda
  maxTimestep = std::min( dx/l_maximumWaveSpeed, dy/l_maximumWaveSpeed );

  // CFL = 0.5
  maxTimestep *= (float)0.4;
}

/**
 * Update the cells with a given global time step.
 *
 * @param i_deltaT time step size.
 */
void SWE_WavePropagationBlockCuda::updateUnknowns(const float i_deltaT) {
  //! definition of one CUDA-block. Typical size are 8*8 or 16*16
  dim3 dimBlock(TILE_SIZE,TILE_SIZE);

  //! definition of the CUDA-grid.
  dim3 dimGrid(nx/TILE_SIZE,ny/TILE_SIZE);

  // assert a valid tile size
  assert(nx%TILE_SIZE==0);
  assert(ny%TILE_SIZE==0);

  // compute the update width.
  float l_updateWidthX = i_deltaT / dx;
  float l_updateWidthY = i_deltaT / dy;

  // update the unknowns (global time step)
  updateUnknownsKernel<<<dimGrid,dimBlock>>>( hNetUpdatesLeftD, hNetUpdatesRightD,
                                              huNetUpdatesLeftD, huNetUpdatesRightD,
                                              hNetUpdatesBelowD, hNetUpdatesAboveD,
                                              hvNetUpdatesBelowD, hvNetUpdatesAboveD,
                                              hd, hud, hvd,
                                              l_updateWidthX, l_updateWidthY,
                                              nx, ny);

  // synchronize the copy layer for MPI communication
  #ifdef USEMPI
  synchCopyLayerBeforeRead();
  #endif
}
