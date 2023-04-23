/**
 * @file
 * This file is part of SWE.
 *
 * @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de,
 * http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
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

#include <cassert>

#include "Block.cuh"
#include "WavePropagationBlock.cuh"
#include "WavePropagationBlockKernels.cuh"

#include "Tools/Logger.hpp"

// CUDA-C includes
#include <cuda.h>
#include <cuda_runtime_api.h>

// Thrust library (used for the final maximum reduction in the method computeNumericalFluxes(...))
#include <thrust/device_vector.h>
#include <thrust/extrema.h>

namespace cuda
{
    /**
     * Constructor of a WavePropagationBlock.
     *
     * Allocates the variables for the simulation:
     *   Please note: The definition of indices changed in contrast to the CPU-Implementation.
     *
     *   unknowns H_,HU_,HV_,B_ stored on the CUDA device are defined for grid indices [0,..,nx_+1]*[0,..,ny_+1] (->
     * Abstract class Block)
     *     -> computational domain is [1,..,nx_]*[1,..,ny_]
     *     -> plus ghost cell layer
     *
     *   net-updates are defined for edges with indices [0,..,nx_]*[0,..,ny_] for horizontal and vertical edges for
     * simplicity (one layer is not necessary).
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
    WavePropagationBlock::WavePropagationBlock(int nx, int ny, RealType dx, RealType dy):
        Block(nx, ny, dx, dy)
    {
        // assert a valid tile size
        ASSERTION(
            nx_ % kTileSize == 0,
            "Decomposition of nx_= " << nx_ << " to kTileSize= " << kTileSize << " is " << nx_ % kTileSize
        );
        ASSERTION(
            ny_ % kTileSize == 0,
            "Decomposition of ny_= " << ny_ << " to kTileSize= " << kTileSize << " is " << ny_ % kTileSize
        );

        // compute the size of one 1D net-update array.
        int sizeOfNetUpdates = (nx_ + 1) * (ny_ + 1) * sizeof(RealType);

        // allocate CUDA memory for the net-updates
        CUDA_MALLOC((void**)&hNetUpdatesLeftD, sizeOfNetUpdates);
        CUDA_MALLOC((void**)&hNetUpdatesRightD, sizeOfNetUpdates);
        CUDA_MALLOC((void**)&huNetUpdatesLeftD, sizeOfNetUpdates);
        CUDA_MALLOC((void**)&huNetUpdatesRightD, sizeOfNetUpdates);
        CUDA_MALLOC((void**)&hNetUpdatesBelowD, sizeOfNetUpdates);
        CUDA_MALLOC((void**)&hNetUpdatesAboveD, sizeOfNetUpdates);
        CUDA_MALLOC((void**)&hvNetUpdatesBelowD, sizeOfNetUpdates);
        CUDA_MALLOC((void**)&hvNetUpdatesAboveD, sizeOfNetUpdates);

        // "2D array" which holds the blockwise maximum wave speeds
        dim3 dimGrid(nx_ / kTileSize, ny_ / kTileSize);

        // size of the maximum wave speed array (dimension of the grid + ghost layers, without the top right block),
        // sizeof(RealType) not included
        int l_sizeMaxWaveSpeeds = ((dimGrid.x + 1) * (dimGrid.y + 1) - 1);
        CUDA_MALLOC((void**)&l_maximumWaveSpeedsD, (l_sizeMaxWaveSpeeds * sizeof(RealType)));
    }

    /**
     * Destructor of a WavePropagationBlock.
     *
     * Frees all of the memory, which was allocated within the constructor.
     * Resets the CUDA device: Useful if error occured and printf is used on the device (buffer).
     */
    WavePropagationBlock::~WavePropagationBlock()
    {
        // free the net-updates memory
        CUDA_FREE(hNetUpdatesLeftD);
        CUDA_FREE(hNetUpdatesRightD);
        CUDA_FREE(huNetUpdatesLeftD);
        CUDA_FREE(huNetUpdatesRightD);

        CUDA_FREE(hNetUpdatesBelowD);
        CUDA_FREE(hNetUpdatesAboveD);
        CUDA_FREE(hvNetUpdatesBelowD);
        CUDA_FREE(hvNetUpdatesAboveD);

        // free the max wave speeds array on the device
        CUDA_FREE(l_maximumWaveSpeedsD);
    }

    /**
     * Compute a single global time step of a given time step width.
     * Remark: The user has to take care about the time step width. No additional check is done. The time step width
     * typically available after the computation of the numerical fluxes (hidden in this method).
     *
     * First the net-updates are computed.
     * Then the cells are updated with the net-updates and the given time step width.
     *
     * @param i_dT time step width in seconds.
     */
    __host__ void WavePropagationBlock::simulateTimestep(RealType i_dT)
    {
        // Compute the numerical fluxes/net-updates in the wave propagation formulation.
        computeNumericalFluxes();

        // Update the unknowns with the net-updates.
        updateUnknowns(i_dT);
    }

    /**
     * perform forward-Euler time steps, starting with simulation time tStart,:
     * until simulation time tEnd is reached;
     * device-global variables H_, HU_, HV_ are updated;
     * unknowns h, hu, hv in main memory are not updated.
     * Ghost layers and bathymetry sources are updated between timesteps.
     * intended as main simulation loop between two checkpoints
     */
    __host__ RealType WavePropagationBlock::simulate(RealType tStart, RealType tEnd)
    {
        RealType t = tStart;
        do
        {
            // set values in ghost cells:
            setGhostLayer();

            // Compute the numerical fluxes/net-updates in the wave propagation formulation.
            computeNumericalFluxes();

            // Update the unknowns with the net-updates.
            updateUnknowns(maxTimeStep_);

            t += maxTimeStep_;
            //     cout << "Simulation at time " << t << endl << flush;
        } while (t < tEnd);

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
    void WavePropagationBlock::computeNumericalFluxes()
    {
        PROFILER_INSTANCE(0);

        /*
         * Initialization.
         */

        /**
         * <b>Row-major vs column-major</b>
         *
         * C/C++ arrays are row-major whereas warps are constructed in column-major order from threads/blocks. To get
         * coalesced memory access in CUDA, we can use a 2-dimensional CUDA structure but we have to switch x and y
         * inside a block.
         *
         * This means, we have to switch threadIdx.x <-> threadIdx.y as well as blockDim.x <-> blockDim.y.
         * Important: blockDim has to be switched for the kernel call as well!
         */

        // ---------------------------------------------------------------------------------------------------------------
        // Grid indexing Ex: nx = 2, ny = 2
        //
        //        i    i+1   i+2   i+3   i+4
        //
        //  j     + --- + --- + --- + --- +
        //        |  g  |  g  |  g  |  g  |
        //  j+1   + --- + --- + --- + --- +
        //        |  g  |  c  |  c  |  g  |
        //  j+2   + --- + --- + --- + --- +    ny+2
        //        |  g  |  c  |  c  |  g  |
        //  j+3   + --- + --- + --- + --- +
        //        |  g  |  g  |  g  |  g  |
        //  j+4   + --- + --- + --- + --- +
        //
        //                  nx+2
        //
        // g = ghost cell
        // c = computational cell
        // CUDA Kernel layout options:
        //
        //         row-based decomposition           column-based decomposition           block-based decomposition
        //           kernel size = 2 * 4                 kernel size = 4 * 2                kernel size = 2 * 2
        //        + === + === + === + === +         ++ --- + --- ++ --- + --- ++         + === + === ++ === + === +
        //        |  g  |  g  |  g  |  g  |         ||  g  |  g  ||  g  |  g  ||         |  g  |  g  ||  g  |  g  |
        //        + --- + --- + --- + --- +         ++ --- + --- ++ --- + --- ++         + --- + --- ++ --- + --- +
        //        |  g  |  c  |  c  |  g  |         ||  g  |  c  ||  c  |  g  ||         |  g  |  c  ||  c  |  g  |
        //        + === + === + === + === +         ++ --- + --- ++ --- + --- ++         + === + === ++ === + === +
        //        |  g  |  c  |  c  |  g  |         ||  g  |  c  ||  c  |  g  ||         |  g  |  c  ||  c  |  g  |
        //        + --- + --- + --- + --- +         ++ --- + --- ++ --- + --- ++         + --- + --- ++ --- + --- +
        //        |  g  |  g  |  g  |  g  |         ||  g  |  g  ||  g  |  g  ||         |  g  |  g  ||  g  |  g  |
        //        + === + === + === + === +         ++ --- + --- ++ --- + --- ++         + === + === ++ === + === +
        // REVIEW Based on the Nsight Compute profiler, we will decide which decomposition is the best for our
        // application.
        // ---------------------------------------------------------------------------------------------------------------

        /**
         * Definition of the "main" CUDA-grid.
         * This grid covers only edges 0..#(edges in x-direction)-2 and 0..#(edges in y-direction)-2.
         *
         * An example with a computational domain of size
         *  nx_ = 24, ny_ = 16
         * with a 1 cell ghost layer would result in a grid with
         *  (nx_+2)*(ny_+2) = (26*18)
         * cells and
         *  (nx_+1)*(ny_+1) = (25*17)
         * edges.
         *
         * The CUDA-blocks (here 8*8) mentioned above would cover all edges except
         * the ones lying between the computational domain and the right/top ghost layer:
         * <pre>
         *                                                          *
         *                                                         **        top ghost layer,
         *                                                        ********   cell ids
         *                        *******************************  **        = (*, ny_+1)
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
         *                  cell ids = (0,*)             cell ids = (nx_+1, *)
         * </pre>
         */
        dim3 dimGrid(nx_ / kTileSize, ny_ / kTileSize);

        //! definition of one CUDA-block. Typical size are 8*8 or 16*16
        dim3 dimBlock(kTileSize , kTileSize);
        /*
         * Compute the net updates for the 'main part and the two 'boundary' parts.
         */
        // compute the net-updates for the "main" part.
        computeNetUpdatesKernel<<<dimGrid, dimBlock>>>(
            H_,
            HU_,
            HV_,
            B_,
            hNetUpdatesLeftD,
            hNetUpdatesRightD,
            huNetUpdatesLeftD,
            huNetUpdatesRightD,
            hNetUpdatesBelowD,
            hNetUpdatesAboveD,
            hvNetUpdatesBelowD,
            hvNetUpdatesAboveD,
            l_maximumWaveSpeedsD,
            nx_,
            ny_
        );

        // compute the "remaining" net updates (edges "simulation domain"/"top ghost layer" and "simulation
        // domain"/"right ghost layer") edges between cell nx_ and ghost layer nx_+1
        dim3 dimRightBlock(kTileSize, 1);
        dim3 dimRightGrid(1, ny_ / kTileSize);
        computeNetUpdatesKernel<<<dimRightGrid, dimRightBlock>>>(
            H_,
            HU_,
            HV_,
            B_,
            hNetUpdatesLeftD,
            hNetUpdatesRightD,
            huNetUpdatesLeftD,
            huNetUpdatesRightD,
            hNetUpdatesBelowD,
            hNetUpdatesAboveD,
            hvNetUpdatesBelowD,
            hvNetUpdatesAboveD,
            l_maximumWaveSpeedsD,
            nx_,
            ny_,
            nx_,
            0,
            dimGrid.x,
            0
        );

        // edges between cell ny_ and ghost layer ny_+1
        dim3 dimTopBlock(1, kTileSize);
        dim3 dimTopGrid(nx_ / kTileSize, 1);
        computeNetUpdatesKernel<<<dimTopGrid, dimTopBlock>>>(
            H_,
            HU_,
            HV_,
            B_,
            hNetUpdatesLeftD,
            hNetUpdatesRightD,
            huNetUpdatesLeftD,
            huNetUpdatesRightD,
            hNetUpdatesBelowD,
            hNetUpdatesAboveD,
            hvNetUpdatesBelowD,
            hvNetUpdatesAboveD,
            l_maximumWaveSpeedsD,
            nx_,
            ny_,
            0,
            ny_,
            0,
            dimGrid.y
        );

        /*
         * Finalize (max reduction of the maximumWaveSpeeds-array.)
         *
         * The Thrust library is used in this step.
         * An optional kernel could be written for the maximum reduction.
         */
        // Thrust pointer to the device array
        thrust::device_ptr<RealType> l_thrustDevicePointer(l_maximumWaveSpeedsD);

        int l_sizeMaxWaveSpeeds = ((dimGrid.x + 1) * (dimGrid.y + 1) - 1);
        // use Thrusts max_element-function for the maximum reduction
        thrust::device_ptr<RealType> l_thrustDevicePointerMax = thrust::max_element(
            l_thrustDevicePointer, l_thrustDevicePointer + l_sizeMaxWaveSpeeds
        );

        // get the result from the device
        RealType l_maximumWaveSpeed = l_thrustDevicePointerMax[0];

        // set the maximum time step for this WavePropagationBlock
        maxTimeStep_ = std::min(dx_ / l_maximumWaveSpeed, dy_ / l_maximumWaveSpeed);

        // CFL = 0.5
        maxTimeStep_ *= (RealType)0.4;
    }

    /**
     * Update the cells with a given global time step.
     *
     * @param i_deltaT time step size.
     */
    void WavePropagationBlock::updateUnknowns(const RealType i_deltaT)
    {
        PROFILER_INSTANCE(0);

        //! definition of one CUDA-block. Typical size are 8*8 or 16*16
        dim3 dimBlock(kTileSize, kTileSize);

        //! definition of the CUDA-grid.
        dim3 dimGrid(nx_ / kTileSize, ny_ / kTileSize);

        // assert a valid tile size
        assert(nx_ % kTileSize == 0);
        assert(ny_ % kTileSize == 0);

        // compute the update width.
        RealType l_updateWidthX = i_deltaT / dx_;
        RealType l_updateWidthY = i_deltaT / dy_;

        // update the unknowns (global time step)
        updateUnknownsKernel<<<dimGrid, dimBlock>>>(
            hNetUpdatesLeftD,
            hNetUpdatesRightD,
            huNetUpdatesLeftD,
            huNetUpdatesRightD,
            hNetUpdatesBelowD,
            hNetUpdatesAboveD,
            hvNetUpdatesBelowD,
            hvNetUpdatesAboveD,
            H_,
            HU_,
            HV_,
            l_updateWidthX,
            l_updateWidthY,
            nx_,
            ny_
        );

// FIXME[epic=SWE,seq=53] mpi for cuda not communicating ghost layer correctly
#ifdef ENABLE_MPI
#if !defined(ENABLE_CUDA_UMA_ALLOCATOR)
        // synchronize the copy layer for MPI communication
        synchBeforeRead();
#endif // !defined(ENABLE_CUDA_UMA_ALLOCATOR)

#endif // ENABLE_MPI
    }
} // namespace cuda
