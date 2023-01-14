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
 * Implementation of Blocks::Block that uses solvers in the wave propagation formulation.
 */

#pragma once

#include "Block.hpp"

#ifdef WITH_SOLVER_HYBRID
#include "HybridSolver.hpp"
#elif defined(WITH_SOLVER_FWAVE)
#include "FWaveSolver.hpp"
#elif defined(WITH_SOLVER_AUGRIE)
#include "AugRieSolver.hpp"
#else
#warning WavePropagationBlock should only be used with Riemann solvers (FWave, AugRie or Hybrid)
#endif

namespace Blocks {

  /**
   * Blocks::WavePropagationBlock is an implementation of the Blocks::Block abstract class.
   * It uses a wave propagation solver which is defined with the pre-compiler flag WITH_SOLVER_ (see above).
   *
   * Possible wave propagation solvers are:
   *  F-Wave, Approximate Augmented Riemann, Hybrid (f-wave + augmented).
   *  (details can be found in the corresponding source files)
   */
  class WavePropagationBlock: public Block {
  private:
#ifdef WITH_SOLVER_HYBRID
    //! Hybrid solver (f-wave + augmented)
    Solvers::HybridSolver<RealType> wavePropagationSolver_;
#elif defined(WITH_SOLVER_FWAVE)
    //! F-wave Riemann solver
    Solvers::FWaveSolver<RealType> wavePropagationSolver_;
#elif defined(WITH_SOLVER_AUGRIE)
    //! Approximate Augmented Riemann solver
    Solvers::AugRieSolver<RealType> wavePropagationSolver_;
#endif

    //! net-updates for the heights of the cells on the left sides of the vertical edges.
    Tools::Float2D<RealType> hNetUpdatesLeft_;
    //! net-updates for the heights of the cells on the right sides of the vertical edges.
    Tools::Float2D<RealType> hNetUpdatesRight_;

    //! net-updates for the x-momentums of the cells on the left sides of the vertical edges.
    Tools::Float2D<RealType> huNetUpdatesLeft_;
    //! net-updates for the x-momentums of the cells on the right sides of the vertical edges.
    Tools::Float2D<RealType> huNetUpdatesRight_;

    //! net-updates for the heights of the cells below the horizontal edges.
    Tools::Float2D<RealType> hNetUpdatesBelow_;
    //! net-updates for the heights of the cells above the horizontal edges.
    Tools::Float2D<RealType> hNetUpdatesAbove_;

    //! net-updates for the y-momentums of the cells below the horizontal edges.
    Tools::Float2D<RealType> hvNetUpdatesBelow_;
    //! net-updates for the y-momentums of the cells above the horizontal edges.
    Tools::Float2D<RealType> hvNetUpdatesAbove_;

  public:
    /**
     * Constructor of a Blocks::WavePropagationBlock.
     *
     * Allocates the variables for the simulation:
     *   unknowns h,hu,hv,b are defined on grid indices [0,..,nx+1]*[0,..,ny+1] (-> Abstract class Blocks::Block)
     *     -> computational domain is [1,..,nx]*[1,..,ny]
     *     -> plus ghost cell layer
     *
     *   net-updates are defined for edges with indices [0,..,nx]*[0,..,ny-1]
     *   or [0,..,nx-1]*[0,..,ny] (for horizontal/vertical edges)
     *
     *   A left/right net update with index (i-1,j-1) is located on the edge between
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
     *   NetUpdatesLeft(i-1,j-1)
     *             or
     *   NetUpdatesRight(i-1,j-1)
     * </pre>
     *
     *   A below/above net update with index (i-1, j-1) is located on the edge between
     *   cells with index (i, j-1) and (i,j):
     * <pre>
     *   ***********
     *   *         *
     *   * (i, j)  *   *
     *   *         *  **  NetUpdatesBelow(i-1,j-1)
     *   *********** *****         or
     *   *         *  **  NetUpdatesAbove(i-1,j-1)
     *   * (i,j-1) *   *
     *   *         *
     *   ***********
     * </pre>
     */
    WavePropagationBlock(int nx, int ny, RealType dx, RealType dy);
    WavePropagationBlock(
      int nx, int ny, RealType dx, RealType dy,
      Tools::Float2D<RealType>& h,
      Tools::Float2D<RealType>& hu,
      Tools::Float2D<RealType>& hv
    );
    ~WavePropagationBlock() override = default;

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