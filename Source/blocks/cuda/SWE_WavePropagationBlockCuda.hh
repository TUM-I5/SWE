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

#ifndef SWEWAVEPROPAGATIONBLOCKCUDA_HH_
#define SWEWAVEPROPAGATIONBLOCKCUDA_HH_

#include <cassert>
#include "SWE_BlockCUDA.hh"

/**
 * SWE_WavePropagationBlockCuda is an implementation of the SWE_BlockCuda abstract class.
 * It uses a wave propagation solver which is defined with the pre-compiler flag WAVE_PROPAGATION_SOLVER (see above).
 *
 * Possible wave propagation solvers are:
 *  F-Wave, <strike>Approximate Augmented Riemann, Hybrid (f-wave + augmented).</strike>
 *  (details can be found in the corresponding source files)
 */
class SWE_WavePropagationBlockCuda: public SWE_BlockCUDA {
  //private:
    //! "2D array" which holds the net-updates for the water height (wave propagating to the left).
    float* hNetUpdatesLeftD;
    //! "2D array" which holds the net-updates for the water height (wave propagating to the right).
    float* hNetUpdatesRightD;

    //! "2D array" which holds the net-updates for the momentum in x-direction (wave propagating to the left).
    float* huNetUpdatesLeftD;
    //! "2D array" which holds the net-updates for the momentum in x-direction (wave propagating to the right).
    float* huNetUpdatesRightD;


    //! "2D array" which holds the net-updates for the water height (wave propagating to the top).
    float* hNetUpdatesBelowD;
    //! "2D array" which holds the net-updates for the water height (wave propagating to the bottom).
    float* hNetUpdatesAboveD;

    //! "2D array" which holds the net-updates for the momentum in y-direction (wave propagating to the top).
    float* hvNetUpdatesBelowD;
    //! "2D array" which holds the net-updates for the momentum in y-direction (wave propagating to the bottom).
    float* hvNetUpdatesAboveD;

  public:
    // constructor of SWE_WavePropagationBlockCuda
    SWE_WavePropagationBlockCuda(
    		int l_nx, int l_ny,
    		float l_dx, float l_dy);

    // destructor of SWE_WavePropagationBlockCuda
    ~SWE_WavePropagationBlockCuda();

    // compute a single time step (net-updates + update of the cells).
    void simulateTimestep( float i_dT );

    // simulate multiple time steps (start and end time provided as parameters)
    float simulate(float, float);

    // TODO: not implemented, max time step reduction is done in each call of computeNumericalFluxes(...)
    //void computeMaxTimestep() {
    //  assert(false);
    //};

    // compute the numerical fluxes (net-update formulation here).
    void computeNumericalFluxes();

    // compute the new cell values.
    void updateUnknowns(const float i_deltaT);
};

#endif /* SWEWAVEPROPAGATIONBLOCKCUDA_HH_ */
