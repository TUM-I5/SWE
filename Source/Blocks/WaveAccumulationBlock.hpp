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

#pragma once

#include "Block.hpp"

// *** Blocks::WaveAccumulationBlock only supports the following wave propagation solvers:
//  - Approximate Augmented Riemann solver (functional implementation: AugRieFunSolver)
//  - F-Wave (vectorized implementation: FWaveVecSolver)
#if defined(WITH_SOLVER_AUGRIE)
#include "Solvers/AugRieFunSolver.hpp"
#elif defined(WITH_SOLVER_HLLE)
#include "Solvers/HLLEFunSolver.hpp"
#elif defined(WITH_SOLVER_FWAVE)
  #if defined(ENABLE_VECTORIZATION) && defined(ENABLE_VECTORIZATION_WITH_SIMD)
    #include "Solvers/FWaveSIMDsolver.hpp"
  #else
    #include "Solvers/FWaveSolver.hpp"
  #endif
#else
#error Chosen wave propagation solver not supported by WaveAccumulationBlock
#endif

namespace Blocks {

  /**
   * Blocks::WaveAccumulationBlock is an implementation of the Blocks::Block abstract class.
   * It uses a wave propagation solver which is defined with the pre-compiler flag WITH_SOLVER_.
   *
   * Possible wave propagation solvers are:
   *  F-Wave, Approximate Augmented Riemann, Hybrid (f-wave + augmented).
   *  (details can be found in the corresponding source files)
   */
  class WaveAccumulationBlock: public Block {
#ifdef WITH_SOLVER_AUGRIE
    //! Approximate Augmented Riemann solver
    Solvers::AugRieFunSolver<RealType> wavePropagationSolver_;
#elif defined(WITH_SOLVER_HLLE)
    //! Vectorized FWave solver
    Solvers::HLLEFunSolver<RealType> wavePropagationSolver_;
#elif defined(WITH_SOLVER_FWAVE)
    //! Approximate Augmented Riemann solver
    #if defined(ENABLE_VECTORIZATION) && defined(ENABLE_VECTORIZATION_WITH_SIMD)
      Solvers::FWaveSIMDsolver<RealType> wavePropagationSolver_;
    #else
      Solvers::FWaveSolver<RealType> wavePropagationSolver_;
    #endif
#endif

    //! net-updates for the heights of the cells (for accumulation)
    Tools::Float2D<RealType> hNetUpdates_;

    //! net-updates for the x-momentums of the cells (for accumulation)
    Tools::Float2D<RealType> huNetUpdates_;

    //! net-updates for the y-momentums of the cells (for accumulation)
    Tools::Float2D<RealType> hvNetUpdates_;

  public:
    /**
     * Constructor of a Blocks::WaveAccumulationBlock.
     *
     * Allocates the variables for the simulation:
     *   unknowns h,hu,hv,b are defined on grid indices [0,..,nx+1]*[0,..,ny+1] (-> Abstract class Blocks::Block)
     *     -> computational domain is [1,..,nx]*[1,..,ny]
     *     -> plus ghost cell layer
     *
     * Similar, all net-updates are defined as cell-local variables with indices [0,..,nx+1]*[0,..,ny+1],
     * however, only values on [1,..,nx]*[1,..,ny] are used (i.e., ghost layers are not accessed).
     * Net updates are intended to hold the accumulated(!) net updates computed on the edges.
     */
    WaveAccumulationBlock(int nx, int ny, RealType dx, RealType dy);
    ~WaveAccumulationBlock() override = default;

    /**
     * Compute net updates for the block.
     * The member variable #maxTimestep will be updated with the
     * maximum allowed time step size
     */
    void computeNumericalFluxes() override;

    /**
     * Updates the unknowns with the already computed net-updates.
     *
     * @param dt time step width used in the update.
     */
    void updateUnknowns(RealType dt) override;
  };

} // namespace Blocks