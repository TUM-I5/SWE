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

#include <cassert>
#include <string>
#include <limits>
#include "SWE_Block.hh"
#include "SWE_WavePropagationBlock.hh"
#ifdef LOOP_OPENMP
#include <omp.h>
#endif

/**
 * Constructor of a SWE_WavePropagationBlock.
 *
 * Allocates the variables for the simulation:
 *   unknowns h,hu,hv,b are defined on grid indices [0,..,nx+1]*[0,..,ny+1] (-> Abstract class SWE_Block)
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
SWE_WavePropagationBlock::SWE_WavePropagationBlock():
  SWE_Block(),
  hNetUpdatesLeft  (nx+1, ny),
  hNetUpdatesRight (nx+1, ny),
  huNetUpdatesLeft (nx+1, ny),
  huNetUpdatesRight(nx+1, ny),
  #if WAVE_PROPAGATION_SOLVER==3
  hvNetUpdatesLeft (nx+1, ny),
  hvNetUpdatesRight(nx+1, ny),
  #endif

  hNetUpdatesBelow (nx, ny+1),
  hNetUpdatesAbove (nx, ny+1),
  #if WAVE_PROPAGATION_SOLVER==3
  huNetUpdatesBelow(nx, ny+1),
  huNetUpdatesAbove(nx, ny+1),
  #endif
  hvNetUpdatesBelow(nx, ny+1),
  hvNetUpdatesAbove(nx, ny+1)
{}

/**
 * Compute net updates for the block.
 * The member variable #maxTimestep will be updated with the 
 * maximum allowed time step size
 */
void SWE_WavePropagationBlock::computeNumericalFluxes() {
  //maximum (linearized) wave speed within one iteration
  float maxWaveSpeed = (float) 0.;

  //compute the net-updates for the vertical edges
  #ifdef LOOP_OPENMP
  #pragma omp parallel
  {
  float l_maxWaveSpeed = (float) 0.;
  solver::Hybrid<float> wavePropagationSolver;
  #pragma omp for
  #endif
  for(int i = 1; i < nx+2; i++) {
#if  WAVE_PROPAGATION_SOLVER==4
#ifdef VECTORIZE
	#pragma simd
#endif
#endif
    for(int j = 1; j < ny+1; j++) {
      float maxEdgeSpeed;
      #if WAVE_PROPAGATION_SOLVER!=3
      wavePropagationSolver.computeNetUpdates( h[i-1][j], h[i][j],
                                               hu[i-1][j], hu[i][j],
                                               b[i-1][j], b[i][j],
                                               hNetUpdatesLeft[i-1][j-1], hNetUpdatesRight[i-1][j-1],
                                               huNetUpdatesLeft[i-1][j-1], huNetUpdatesRight[i-1][j-1],
                                               maxEdgeSpeed );
     #else
     //TODO: implement again.
     assert(false);
     #endif

      #ifdef LOOP_OPENMP
      //update the maximum wave speed
      l_maxWaveSpeed = std::max(l_maxWaveSpeed, maxEdgeSpeed);
      #else
      //update the maximum wave speed
      maxWaveSpeed = std::max(maxWaveSpeed, maxEdgeSpeed);
      #endif
    }
  }

  //compute the net-updates for the horizontal edges
  #ifdef LOOP_OPENMP
  #pragma omp for
  #endif
  for(int i = 1; i < nx+1; i++) {
#if  WAVE_PROPAGATION_SOLVER==4
#ifdef VECTORIZE
	#pragma simd
#endif
#endif
    for(int j = 1; j < ny+2; j++) {
      float maxEdgeSpeed;
      #if WAVE_PROPAGATION_SOLVER!=3
      wavePropagationSolver.computeNetUpdates( h[i][j-1], h[i][j],
                                               hv[i][j-1], hv[i][j],
                                               b[i][j-1], b[i][j],
                                               hNetUpdatesBelow[i-1][j-1], hNetUpdatesAbove[i-1][j-1],
                                               hvNetUpdatesBelow[i-1][j-1], hvNetUpdatesAbove[i-1][j-1],
                                               maxEdgeSpeed );
      #else
      //TODO: implement again.
      assert(false);
      #endif

      #ifdef LOOP_OPENMP
      //update the maximum wave speed
      l_maxWaveSpeed = std::max(l_maxWaveSpeed, maxEdgeSpeed);
      #else
      //update the maximum wave speed
      maxWaveSpeed = std::max(maxWaveSpeed, maxEdgeSpeed);
      #endif
    }
  }
  #ifdef LOOP_OPENMP
  #pragma omp critical
  {
  maxWaveSpeed = std::max(l_maxWaveSpeed, maxWaveSpeed);
  }
  } // end of parallel for block
  #endif

  if(maxWaveSpeed > 0.00001) { //TODO zeroTol
    //compute the time step width
    //CFL-Codition
    //(max. wave speed) * dt / dx < .5
    // => dt = .5 * dx/(max wave speed)
    maxTimestep = std::min( dx/maxWaveSpeed, dy/maxWaveSpeed );

//    #if WAVE_PROPAGATION_SOLVER!=3
    maxTimestep *= (float) .4; //CFL-number = .5
//    #else
//    dt *= (float) .8; //CFL-number = 1. (wave limiters)
//    #endif
  }
  else
    maxTimestep = std::numeric_limits<float>::max(); //might happen in dry cells

}

/**
 * Updates the unknowns with the already computed net-updates.
 *
 * @param dt time step width used in the update.
 */
void SWE_WavePropagationBlock::updateUnknowns(float dt) {
  //update cell averages with the net-updates
  #ifdef LOOP_OPENMP
  #pragma omp parallel for
  #endif
  for(int i = 1; i < nx+1; i++) {
#ifdef VECTORIZE
	  // Tell the compiler that he can safely ignore all dependencies in this loop
	  #pragma ivdep
#endif
      for(int j = 1; j < ny+1; j++) {
        h[i][j] -=   dt/dx * (hNetUpdatesRight[i-1][j-1] + hNetUpdatesLeft[i][j-1])
                   + dt/dy * (hNetUpdatesAbove[i-1][j-1] + hNetUpdatesBelow[i-1][j]);
        hu[i][j] -= dt/dx * (huNetUpdatesRight[i-1][j-1] + huNetUpdatesLeft[i][j-1]);
        hv[i][j] -= dt/dy * (hvNetUpdatesAbove[i-1][j-1] + hvNetUpdatesBelow[i-1][j]);
        #if WAVE_PROPAGATION_SOLVER==3
        hv[i][j] -= dt/dx * (hvNetUpdatesRight[i-1][j-1] + hvNetUpdatesLeft[i][j-1]);
        hu[i][j] -= dt/dy * (huNetUpdatesAbove[i-1][j-1] + huNetUpdatesBelow[i-1][j]);
        #endif

        //TODO: dryTol
        if(h[i][j] < 0) {
#ifndef NDEBUG
          // Only print this warning when debug is enabled.
          // Otherwise we cannot vectorize this loop
          if(h[i][j] < -0.1) {
            std::cerr << "Warning, negative height: (i,j)=(" << i << "," << j << ")=" << h[i][j] << std::endl;
            std::cerr << "         b: " << b[i][j] << std::endl;
          }
#endif
          //zero (small) negative depths
          h[i][j] = hu[i][j] = hv[i][j] = 0.;
        }
        else if(h[i][j] < 0.1)
          hu[i][j] = hv[i][j] = 0.; //no water, no speed!
      }
  }
}

/**
 * Update the bathymetry values with the displacement corresponding to the current time step.
 *
 * @param i_asagiScenario the corresponding ASAGI-scenario
 */
#ifdef DYNAMIC_DISPLACEMENTS
bool SWE_WavePropagationBlock::updateBathymetryWithDynamicDisplacement(scenarios::Asagi &i_asagiScenario, const float i_time) {
  if (!i_asagiScenario.dynamicDisplacementAvailable(i_time))
    return false;

  // update the bathymetry
  for(int i=0; i<=nx+1; i++) {
    for(int j=0; j<=ny+1; j++) {
      b[i][j] = i_asagiScenario.getBathymetryAndDynamicDisplacement( offsetX + (i-0.5f)*dx,
                                                                     offsetY + (j-0.5f)*dy,
                                                                     i_time
                                                                   );
    }
  }
  return true;
}
#endif

/**
 * Executes a single timestep.
 *  * compute net updates for every edge
 *  * update cell values with the net updates
 *
 * @param dt	time step width of the update
 */
void SWE_WavePropagationBlock::simulateTimestep(float dt) {
  computeNumericalFluxes();
  updateUnknowns(dt);
}

/**
 * Runs the simulation until i_tEnd is reached.
 *
 * @param i_tStart time when the simulation starts
 * @param i_tEnd  time when the simulation should end
 * @return time we reached after the last update step, in general a bit later than i_tEnd
 */
float SWE_WavePropagationBlock::simulate(float i_tStart,float i_tEnd) {
  float t = i_tStart;
  do {
     //set values in ghost cells
     setGhostLayer();

     // compute net updates for every edge
     computeNumericalFluxes();
     //execute a wave propagation time step
     updateUnknowns(maxTimestep);
     t += maxTimestep;

     std::cout << "Simulation at time " << t << std::endl << std::flush;
  } while(t < i_tEnd);

  return t;
}

#if WAVE_PROPAGATION_SOLVER==0
#ifndef NDEBUG
#ifndef LOOP_OPENMP
/**
 * Gets the current solver statistics for a hybrid solver.
 *
 * @param o_counterFWave will be set to: times the f-Wave solver was used.
 * @param o_counterAugRie will be set to: times the Augmented Riemann solver was used.
 */
void SWE_WavePropagationBlock::getSolverStats( long &o_counterFWave, long &o_counterAugRie ) {
  wavePropagationSolver.getStats( o_counterFWave, o_counterAugRie );
}
#endif
#endif
#endif
