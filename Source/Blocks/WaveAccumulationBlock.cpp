/**
 * @file
 * This file is part of SWE.
 *
 * @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de,
 * http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
 * @author Michael Bader (bader AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Michael_Bader)
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
 * Blocks::Block, which uses solvers in the wave propagation formulation.
 */

#include "WaveAccumulationBlock.hpp"

#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

#include <iostream>

Blocks::WaveAccumulationBlock::WaveAccumulationBlock(int nx, int ny, RealType dx, RealType dy):
  Block(nx, ny, dx, dy),
  hNetUpdates_(nx + 2, ny + 2),
  huNetUpdates_(nx + 2, ny + 2),
  hvNetUpdates_(nx + 2, ny + 2) {}

void Blocks::WaveAccumulationBlock::computeNumericalFluxes() {
  RealType dxInv = RealType(1.0) / dx_;
  RealType dyInv = RealType(1.0) / dy_;

  // Maximum (linearized) wave speed within one iteration
  RealType maxWaveSpeed = RealType(0.0);

#ifdef ENABLE_OPENMP
#pragma omp parallel
  {
    // Thread-local maximum wave speed:
    RealType maxWaveSpeedLocal = RealType(0.0);

    // Use OpenMP for the outer loop
#pragma omp for
#endif
    // Compute the net-updates for the vertical edges
    for (int i = 1; i < nx_ + 2; i++) {
      const int nyEnd = ny_ + 1; // Compiler might refuse to vectorize j-loop without this ...
#if defined(ENABLE_VECTORIZATION) && defined(ENABLE_OPENMP) // Vectorize the inner loop
#pragma omp simd reduction(max : maxWaveSpeedLocal)
#elif defined(ENABLE_OPENMP)
#pragma omp reduction(max : maxWaveSpeedLocal)
#elif defined(ENABLE_VECTORIZATION)
#pragma omp simd reduction(max : maxWaveSpeed)
#endif
      for (int j = 1; j < nyEnd; j++) {
        RealType maxEdgeSpeed = RealType(0.0);
        RealType hNetUpLeft = 0.0, hNetUpRight = 0.0;
        RealType huNetUpLeft = 0.0, huNetUpRight = 0.0;
        wavePropagationSolver_.computeNetUpdates(
          h_[i - 1][j],
          h_[i][j],
          hu_[i - 1][j],
          hu_[i][j],
          b_[i - 1][j],
          b_[i][j],
          hNetUpLeft,
          hNetUpRight,
          huNetUpLeft,
          huNetUpRight,
          maxEdgeSpeed
        );

        // Accumulate net updates to cell-wise net updates for h and hu
        hNetUpdates_[i - 1][j] += dxInv * hNetUpLeft;
        huNetUpdates_[i - 1][j] += dxInv * huNetUpLeft;
        hNetUpdates_[i][j] += dxInv * hNetUpRight;
        huNetUpdates_[i][j] += dxInv * huNetUpRight;

#ifdef ENABLE_OPENMP
        // Update the thread-local maximum wave speed
        maxWaveSpeedLocal = std::max(maxWaveSpeed, maxEdgeSpeed);
#else
      // Update the maximum wave speed
      maxWaveSpeed = std::max(maxWaveSpeed, maxEdgeSpeed);
#endif
      }
    }

#ifdef ENABLE_OPENMP
#pragma omp for
#endif
    // Compute the net-updates for the horizontal edges
    for (int i = 1; i < nx_ + 1; i++) {
      const int nyEnd = ny_ + 2;
#if defined(ENABLE_VECTORIZATION) && defined(ENABLE_OPENMP)
#pragma omp simd reduction(max : maxWaveSpeedLocal)
#elif defined(ENABLE_OPENMP)
#pragma omp reduction(max : maxWaveSpeedLocal)
#elif defined(ENABLE_VECTORIZATION)
#pragma omp simd reduction(max : maxWaveSpeed)
#endif
      for (int j = 1; j < nyEnd; j++) {
        RealType maxEdgeSpeed = 0.0;
        RealType hNetUpDow = 0.0, hNetUpUpw = 0.0;
        RealType hvNetUpDow = 0.0, hvNetUpUpw = 0.0;
        wavePropagationSolver_.computeNetUpdates(
          h_[i][j - 1],
          h_[i][j],
          hv_[i][j - 1],
          hv_[i][j],
          b_[i][j - 1],
          b_[i][j],
          hNetUpDow,
          hNetUpUpw,
          hvNetUpDow,
          hvNetUpUpw,
          maxEdgeSpeed
        );

        // Accumulate net updates to cell-wise net updates for h and hu
        hNetUpdates_[i][j - 1] += dyInv * hNetUpDow;
        hvNetUpdates_[i][j - 1] += dyInv * hvNetUpDow;
        hNetUpdates_[i][j] += dyInv * hNetUpUpw;
        hvNetUpdates_[i][j] += dyInv * hvNetUpUpw;

#ifdef ENABLE_OPENMP
        // Update the thread-local maximum wave speed
        maxWaveSpeedLocal = std::max(maxWaveSpeed, maxEdgeSpeed);
#else
      // Update the maximum wave speed
      maxWaveSpeed = std::max(maxWaveSpeed, maxEdgeSpeed);
#endif
      }
    }

#ifdef ENABLE_OPENMP
#pragma omp critical
    { maxWaveSpeed = std::max(maxWaveSpeedLocal, maxWaveSpeed); }

  } // #pragma omp parallel
#endif

  if (maxWaveSpeed > 0.00001) {
    // Compute the time step width
    maxTimeStep_ = std::min(dx_ / maxWaveSpeed, dy_ / maxWaveSpeed);

    // Reduce maximum time step size by "safety factor"
    maxTimeStep_ *= RealType(0.4); // CFL-number = 0.5
  } else {
    // Might happen in dry cells
    maxTimeStep_ = std::numeric_limits<float>::max();
  }
}

void Blocks::WaveAccumulationBlock::updateUnknowns(RealType dt) {
  // Update cell averages with the net-updates

#ifdef SWE_USE_OPENMP
#pragma omp parallel for
#endif
  for (int i = 1; i < nx_ + 1; i++) {
    const int nyEnd = ny_ + 1;
#ifdef ENABLE_VECTORIZATION
    // Tell the compiler that he can safely ignore all dependencies in this loop
#pragma omp simd
#endif
    for (int j = 1; j < nyEnd; j++) {
      h_[i][j] -= dt * hNetUpdates_[i][j];
      hu_[i][j] -= dt * huNetUpdates_[i][j];
      hv_[i][j] -= dt * hvNetUpdates_[i][j];

      hNetUpdates_[i][j]  = RealType(0.0);
      huNetUpdates_[i][j] = RealType(0.0);
      hvNetUpdates_[i][j] = RealType(0.0);

      if (h_[i][j] < 0.1) {                    // dryTol
        hu_[i][j] = hv_[i][j] = RealType(0.0); // No water, no speed!
      }

      if (h_[i][j] < 0) {
#ifndef NDEBUG
        // Only print this warning when debug is enabled
        // Otherwise we cannot vectorize this loop
        if (h_[i][j] < -0.1) {
          std::cerr << "Warning, negative height: (i,j)=(" << i << "," << j << ")=" << h_[i][j] << std::endl;
          std::cerr << "         b: " << b_[i][j] << std::endl;
        }
#endif

        // Zero (small) negative depths
        h_[i][j] = RealType(0.0);
      }
    }
  }
}
