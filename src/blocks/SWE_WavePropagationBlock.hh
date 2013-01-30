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
 * SWE_Block, which uses solvers in the wave propagation formulation.
 */

#ifndef SWEWAVEPROPAGATIONBLOCK_HH_
#define SWEWAVEPROPAGATIONBLOCK_HH_

#include "blocks/SWE_Block.hh"
#ifdef DYNAMIC_DISPLACEMENTS
#include "scenarios/Asagi.hpp"
#endif
#include "tools/help.hh"

#include <string>

//which wave propagation solver should be used
//  0: Hybrid
//  1: f-Wave
//  2: Approximate Augmented Riemann solver
//  3: Approximate Augmented Riemann solver which uses the underlying Fortran routines of GeoClaw directly.
#if WAVE_PROPAGATION_SOLVER==1
#include "solvers/FWave.hpp"
#elif WAVE_PROPAGATION_SOLVER==2
#include "solvers/AugRie.hpp"
#elif WAVE_PROPAGATION_SOLVER==3
#include "solvers/AugRieGeoClaw.hpp"
#elif WAVE_PROPAGATION_SOLVER==4
#include "solvers/FWaveVec.hpp"
#else
#include "solvers/Hybrid.hpp"
#endif

/**
 * SWE_WavePropagationBlock is an implementation of the SWE_Block abstract class.
 * It uses a wave propagation solver which is defined with the pre-compiler flag WAVE_PROPAGATION_SOLVER (see above).
 *
 * Possible wave propagation solvers are:
 *  F-Wave, Apprximate Augmented Riemann, Hybrid (f-wave + augmented).
 *  (details can be found in the corresponding source files)
 */
class SWE_WavePropagationBlock: public SWE_Block {
//private:
  //OpenMp: Every task defines it own solver -> SWE_WavePropagationBlock.cpp
#ifndef LOOP_OPENMP
    //specify the wave propagation solver
#if WAVE_PROPAGATION_SOLVER==1
    //! F-wave Riemann solver
    solver::FWave<float> wavePropagationSolver;
#elif WAVE_PROPAGATION_SOLVER==2
    //! Approximate Augmented Riemann solver
    solver::AugRie<float> wavePropagationSolver;
#elif WAVE_PROPAGATION_SOLVER==3
    //! Approximate Augmented Riemann solver
    solver::AugRieGeoClaw<double> wavePropagationSolver;
#elif WAVE_PROPAGATION_SOLVER==4
    //! Vectorized FWave solver
    solver::FWaveVec<float> wavePropagationSolver;
#else
    //! Hybrid solver (f-wave + augmented)
    solver::Hybrid<float> wavePropagationSolver;
#endif
#endif

    //! net-updates for the heights of the cells on the left sides of the vertical edges.
    Float2D hNetUpdatesLeft;
    //! net-updates for the heights of the cells on the right sides of the vertical edges.
    Float2D hNetUpdatesRight;

    //! net-updates for the x-momentums of the cells on the left sides of the vertical edges.
    Float2D huNetUpdatesLeft;
    //! net-updates for the x-momentums of the cells on the right sides of the vertical edges.
    Float2D huNetUpdatesRight;

    #if WAVE_PROPAGATION_SOLVER==3
    //! net-updates for the y-momentums of the cells on the left sides of the vertical edges.
    Float2D hvNetUpdatesLeft;
    //! net-updates for the y-momentums of the cells on the right sides of the vertical edges.
    Float2D hvNetUpdatesRight;
    #endif

    //! net-updates for the heights of the cells below the horizontal edges.
    Float2D hNetUpdatesBelow;
    //! net-updates for the heights of the cells above the horizontal edges.
    Float2D hNetUpdatesAbove;

    //! net-updates for the y-momentums of the cells below the horizontal edges.
    Float2D hvNetUpdatesBelow;
    //! net-updates for the y-momentums of the cells above the horizontal edges.
    Float2D hvNetUpdatesAbove;

    #if WAVE_PROPAGATION_SOLVER==3
    //! net-updates for the x-momentums of the cells below the horizontal edges.
    Float2D huNetUpdatesBelow;
    //! net-updates for the x-momentums of the cells below the horizontal edges.
    Float2D huNetUpdatesAbove;
    #endif

  public:
    //constructor of a SWE_WavePropagationBlock.
    SWE_WavePropagationBlock(int l_nx, int l_ny,
    					float l_dx, float l_dy);

    //executes a single timestep.
    virtual void simulateTimestep(float dt);

    //computes the net-updates for the block
    void computeNumericalFluxes();

    //update the cells
    void updateUnknowns(float dt);

    //runs the simulation until i_tEnd is reached.
    float simulate(float i_tStart, float i_tEnd);

    //updates the bathymetry with the current displacment values
#ifdef DYNAMIC_DISPLACEMENTS
    bool updateBathymetryWithDynamicDisplacement(scenarios::Asagi &i_asagiScenario, float time);
#endif

    /**
     * Destructor of a SWE_WavePropagationBlock.
     *
     * In the case of a hybrid solver (NDEBUG not defined) information about the used solvers will be printed.
     */
    virtual ~SWE_WavePropagationBlock() {}
};

#endif /* SWEWAVEPROPAGATIONBLOCK_HH_ */
