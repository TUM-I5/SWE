/**
 * @file
 * This file is part of SWE.
 *
 * @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
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
 * SWE_Block, which uses solvers in the wave propagation formulation.
 */

#ifndef SWE_WAVEACCUMULATION_BLOCK_HH_
#define SWE_WAVEACCUMULATION_BLOCK_HH_

#include "blocks/SWE_Block.hh"
#ifdef DYNAMIC_DISPLACEMENTS
#include "scenarios/Asagi.hpp"
#endif
#include "tools/help.hh"

#include <string>

// *** SWE_WaveAccumulationBlock only supports the following wave propagation solvers:
//  2: Approximate Augmented Riemann solver (functional implementation: AugRieFun)
//  4: f-Wave (vectorized implementation: FWaveVec) 
#if WAVE_PROPAGATION_SOLVER==2
#include "solvers/AugRieFun.hpp"
#elif WAVE_PROPAGATION_SOLVER==4
#include "solvers/FWaveVec.hpp"
#else
#warning chosen wave propagation solver not supported by SWE_WaveAccumulationBlock
#endif

/**
 * SWE_WaveAccumulationBlock is an implementation of the SWE_Block abstract class.
 * It uses a wave propagation solver which is defined with the pre-compiler flag WAVE_PROPAGATION_SOLVER (see above).
 *
 * Possible wave propagation solvers are:
 *  F-Wave, Apprximate Augmented Riemann, Hybrid (f-wave + augmented).
 *  (details can be found in the corresponding source files)
 */
class SWE_WaveAccumulationBlock: public SWE_Block {

#if WAVE_PROPAGATION_SOLVER==2
    //! Approximate Augmented Riemann solver
    solver::AugRieFun<float> wavePropagationSolver;
#elif WAVE_PROPAGATION_SOLVER==4
    //! Vectorized FWave solver
    solver::FWaveVec<float> wavePropagationSolver;
#endif

    //! net-updates for the heights of the cells (for accumulation)
    Float2D hNetUpdates;

    //! net-updates for the x-momentums of the cells (for accumulation)
    Float2D huNetUpdates;

    //! net-updates for the y-momentums of the cells (for accumulation)
    Float2D hvNetUpdates;

  public:
    //constructor of a SWE_WaveAccumulationBlock.
    SWE_WaveAccumulationBlock(int l_nx, int l_ny, float l_dx, float l_dy);
    //destructor of a SWE_WaveAccumulationBlock.
    virtual ~SWE_WaveAccumulationBlock() {}

    //computes the net-updates for the block
    void computeNumericalFluxes();

    //update the cells
    void updateUnknowns(float dt);

    //updates the bathymetry with the current displacment values
#ifdef DYNAMIC_DISPLACEMENTS
    bool updateBathymetryWithDynamicDisplacement(scenarios::Asagi &i_asagiScenario, float time);
#endif

};

#endif /* SWE_WAVEACCUMULATION_BLOCK_HH_ */
