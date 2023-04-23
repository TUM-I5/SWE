/**
 * @file
 * This file is part of SWE.
 *
 * @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
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
 * CUDA Kernels for a SWE_Block, which uses solvers in the wave propagation formulation.
 */

#ifndef SWEWAVEPROPAGATIONBLOCKCUDAKERNELS_HH_
#define SWEWAVEPROPAGATIONBLOCKCUDAKERNELS_HH_

namespace cuda
{

    // CUDA-kernel which computes the net-updates
    __global__ void computeNetUpdatesKernel(
        const RealType* i_h,
        const RealType* i_hu,
        const RealType* i_hv,
        const RealType* i_b,
        RealType*       o_hNetUpdatesLeftD,
        RealType*       o_hNetUpdatesRightD,
        RealType*       o_huNetUpdatesLeftD,
        RealType*       o_huNetUpdatesRightD,
        RealType*       o_hNetUpdatesBelowD,
        RealType*       o_hNetUpdatesAboveD,
        RealType*       o_hvNetUpdatesBelowD,
        RealType*       o_hvNetUpdatesAboveD,
        RealType*       o_maximumWaveSpeeds,
        const int       i_nx,
        const int       i_ny,
        const int       i_offsetX      = 0,
        const int       i_offsetY      = 0,
        const int       i_blockOffSetX = 0,
        const int       i_blockOffSetY = 0
    );

    // CUDA-kernel which updates the unknowns
    __global__ void updateUnknownsKernel(
        const RealType* i_hNetUpdatesLeftD,
        const RealType* i_hNetUpdatesRightD,
        const RealType* i_huNetUpdatesLeftD,
        const RealType* i_huNetUpdatesRightD,
        const RealType* i_hNetUpdatesBelowD,
        const RealType* i_hNetUpdatesAboveD,
        const RealType* i_hvNetUpdatesBelowD,
        const RealType* i_hvNetUpdatesAboveD,
        RealType*       io_h,
        RealType*       io_hu,
        RealType*       io_hv,
        const RealType  i_updateWidthX,
        const RealType  i_updateWidthY,
        const int       i_nx,
        const int       i_ny
    );

    // CUDA-kernel which computes the 1D position in an array from a given 2D index
    __device__ inline int computeOneDPositionKernel(const int i_i, const int i_j, const int i_nx);
} // namespace cuda

#endif /* SWEWAVEPROPAGATIONBLOCKCUDAKERNELS_HH_ */
