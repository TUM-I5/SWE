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

#include "SWE_WavePropagationBlockSIMD.hh"

#if WAVE_PROPAGATION_SOLVER==3
#error "Solver not implemented in SWE_WavePropagationBlockSIMD"
#endif // WAVE_PROPAGATION_SOLVER==3

#include <cassert>
#include <string>
#include <limits>

#ifdef LOOP_OPENMP
#include <omp.h>
#endif

/**
 * Constructor of a SWE_WavePropagationBlockSIMD.
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
SWE_WavePropagationBlockSIMD::SWE_WavePropagationBlockSIMD (int l_nx, int l_ny, float l_dx, float l_dy) :
	SWE_Block (l_nx, l_ny, l_dx, l_dy),
	hNetUpdatesLeft (nx + 1, ny),
	hNetUpdatesRight (nx + 1, ny),
	huNetUpdatesLeft (nx + 1, ny),
	huNetUpdatesRight (nx + 1, ny),

	hNetUpdatesBelow (nx, ny + 1),
	hNetUpdatesAbove (nx, ny + 1),
	hvNetUpdatesBelow (nx, ny + 1),
	hvNetUpdatesAbove (nx, ny + 1)
#ifdef COUNTFLOPS
	, // don't forget the comma, to extend the list above
	flops(0),
	time_needed(0.0)
#endif
{
}

/**
 * Compute net updates for the block.
 * The member variable #maxTimestep will be updated with the 
 * maximum allowed time step size
 */
void
SWE_WavePropagationBlockSIMD::computeNumericalFluxes ()
{
#ifdef COUNTFLOPS
#ifdef LOOP_OPENMP
	const double time_begin = omp_get_wtime();
#else
	const double time_begin = clock();
#endif
#endif
	//maximum (linearized) wave speed within one iteration
	float maxWaveSpeed = (float) 0.;

	// compute the loop limits
	const int end_ny_1_1 = ny + 1;
	const int end_ny_1_2 = ny + 2;

#if  WAVE_PROPAGATION_SOLVER==5
	// Note, that ny is used instead of (ny + 1). This is due to the fact, that in the loop below, j starts with 1!
	// So, for SSE, for instance, j takes values 1, 5, 9, 13, ...
	// Now, consider ny == 3
	//
	// Then the code below returns 0 as end of the vector loop, which is correct (!)
	//
	// if, instead (ny + 1) was used, the end of the vector loop is 4 (!)
	// that means, the solver accesses the elements (1, 2, 3, 4), even though only the elements (0, 1, 2, 3) exist (recall, ny == 3) (!)
	const int end_ny_V_1 = ny & (~(VECTOR_LENGTH - 1));
	const int end_ny_V_2 = (ny + 1) & (~(VECTOR_LENGTH - 1));
#endif /*  WAVE_PROPAGATION_SOLVER==5 */

	/***************************************************************************************
	 * compute the net-updates for the vertical edges
	 **************************************************************************************/

#ifdef LOOP_OPENMP
#pragma omp parallel
	{

		float l_maxWaveSpeed = (float) 0.;
//		solver::Hybrid<float> wavePropagationSolver;
#pragma message "augmented Riemann solver was hardcoded set for OpenMP!"
		solver::AugRie_SIMD wavePropagationSolver;

		// Use OpenMP for the outer loop
#pragma omp for
#endif // LOOP_OPENMP
		for (int i = 1; i < nx + 2; i++) {
			int j = 1;

#if  WAVE_PROPAGATION_SOLVER==5 and (not defined VECTOR_NOVEC)
			for (; j < end_ny_V_1; j += VECTOR_LENGTH) {
				float maxEdgeSpeed;

				wavePropagationSolver.computeNetUpdates_SIMD (
					&h[i - 1][j], &h[i][j],
					&hu[i - 1][j], &hu[i][j],
					&b[i - 1][j], &b[i][j],
					&hNetUpdatesLeft[i - 1][j - 1], &hNetUpdatesRight[i - 1][j - 1],
					&huNetUpdatesLeft[i - 1][j - 1], &huNetUpdatesRight[i - 1][j - 1],
					maxEdgeSpeed
				);

#ifdef LOOP_OPENMP
				//update the thread-local maximum wave speed
				l_maxWaveSpeed = std::max (l_maxWaveSpeed, maxEdgeSpeed);
#else // LOOP_OPENMP
				//update the maximum wave speed
				maxWaveSpeed = std::max (maxWaveSpeed, maxEdgeSpeed);
#endif // LOOP_OPENMP
			}
#endif /* WAVE_PROPAGATION_SOLVER==5 */
#if  WAVE_PROPAGATION_SOLVER==4 and defined VECTORIZE
			// Vectorization is currently only possible for the FWaveVec solver
			// Vectorize the inner loop
#pragma simd
#endif // WAVE_PROPAGATION_SOLVER==4 and defined VECTORIZE
			for (; j < end_ny_1_1; ++j) {
				float maxEdgeSpeed;

				wavePropagationSolver.computeNetUpdates (
					h[i - 1][j], h[i][j],
					hu[i - 1][j], hu[i][j],
					b[i - 1][j], b[i][j],
					hNetUpdatesLeft[i - 1][j - 1], hNetUpdatesRight[i - 1][j - 1],
					huNetUpdatesLeft[i - 1][j - 1], huNetUpdatesRight[i - 1][j - 1],
					maxEdgeSpeed
				);

#ifdef LOOP_OPENMP
				//update the thread-local maximum wave speed
				l_maxWaveSpeed = std::max (l_maxWaveSpeed, maxEdgeSpeed);
#else // LOOP_OPENMP
				//update the maximum wave speed
				maxWaveSpeed = std::max (maxWaveSpeed, maxEdgeSpeed);
#endif // LOOP_OPENMP
			}
			assert (j == ny + 1);
		}

	/***************************************************************************************
	 * compute the net-updates for the horizontal edges
	 **************************************************************************************/

#ifdef LOOP_OPENMP
		// Use OpenMP for the outer loop
#pragma omp for
#endif // LOOP_OPENMP
		for (int i = 1; i < nx + 1; i++) {
			int j = 1;
#if  WAVE_PROPAGATION_SOLVER==5 and (not defined VECTOR_NOVEC)
			for (; j < end_ny_V_2; j += VECTOR_LENGTH) {
				float maxEdgeSpeed;

				wavePropagationSolver.computeNetUpdates_SIMD(
					&h[i][j - 1], &h[i][j],
					&hv[i][j - 1], &hv[i][j],
					&b[i][j - 1], &b[i][j],
					&hNetUpdatesBelow[i - 1][j - 1], &hNetUpdatesAbove[i - 1][j - 1],
					&hvNetUpdatesBelow[i - 1][j - 1], &hvNetUpdatesAbove[i - 1][j - 1],
					maxEdgeSpeed
				);

#ifdef LOOP_OPENMP
				//update the thread-local maximum wave speed
				l_maxWaveSpeed = std::max (l_maxWaveSpeed, maxEdgeSpeed);
#else // LOOP_OPENMP
				//update the maximum wave speed
				maxWaveSpeed = std::max (maxWaveSpeed, maxEdgeSpeed);
#endif // LOOP_OPENMP
			}
#endif /* WAVE_PROPAGATION_SOLVER==5 */

#if  WAVE_PROPAGATION_SOLVER==4 and defined VECTORIZE
		// Vectorization is currently only possible for the FWaveVec solver
		// Vectorize the inner loop
#pragma simd
#endif // WAVE_PROPAGATION_SOLVER==4
			for (; j < end_ny_1_2; j++) {
				float maxEdgeSpeed;

				wavePropagationSolver.computeNetUpdates (
					h[i][j - 1], h[i][j],
					hv[i][j - 1], hv[i][j],
					b[i][j - 1], b[i][j],
					hNetUpdatesBelow[i - 1][j - 1], hNetUpdatesAbove[i - 1][j - 1],
					hvNetUpdatesBelow[i - 1][j - 1], hvNetUpdatesAbove[i - 1][j - 1],
					maxEdgeSpeed
				);

#ifdef LOOP_OPENMP
				//update the thread-local maximum wave speed
				l_maxWaveSpeed = std::max (l_maxWaveSpeed, maxEdgeSpeed);
#else // LOOP_OPENMP
				//update the maximum wave speed
				maxWaveSpeed = std::max (maxWaveSpeed, maxEdgeSpeed);
#endif // LOOP_OPENMP
			}
			assert (j = ny + 2);
		}

#ifdef LOOP_OPENMP
#pragma omp critical
		{
			maxWaveSpeed = std::max (l_maxWaveSpeed, maxWaveSpeed);
#ifdef COUNTFLOPS
			flops += wavePropagationSolver.flops;
#endif
		}

	} // #pragma omp parallel
#endif

	if (maxWaveSpeed > 0.00001) {
		//TODO zeroTol

		//compute the time step width
		//CFL-Codition
		//(max. wave speed) * dt / dx < .5
		// => dt = .5 * dx/(max wave speed)

		maxTimestep = std::min (dx / maxWaveSpeed, dy / maxWaveSpeed);

		maxTimestep *= (float) .4; //CFL-number = .5
	} else {
		//might happen in dry cells
		maxTimestep = std::numeric_limits<float>::max ();
	}
#ifdef COUNTFLOPS
#ifdef LOOP_OPENMP
	time_needed += omp_get_wtime() - time_begin;
#else
	time_needed += clock() - time_begin;
#endif
#endif
}

/**
 * Updates the unknowns with the already computed net-updates.
 *
 * @param dt time step width used in the update.
 */
void
SWE_WavePropagationBlockSIMD::updateUnknowns (float dt)
{
	//update cell averages with the net-updates
#ifdef LOOP_OPENMP
#pragma omp parallel for
#endif // LOOP_OPENMP
	for (int i = 1; i < nx + 1; i++) {

#ifdef VECTORIZE
		// Tell the compiler that he can safely ignore all dependencies in this loop
#pragma ivdep
#endif // VECTORIZE
		for (int j = 1; j < ny + 1; j++) {
			h[i][j] -= dt / dx * (hNetUpdatesRight[i - 1][j - 1] + hNetUpdatesLeft[i][j - 1]) + dt / dy * (hNetUpdatesAbove[i - 1][j - 1] + hNetUpdatesBelow[i - 1][j]);
			hu[i][j] -= dt / dx * (huNetUpdatesRight[i - 1][j - 1] + huNetUpdatesLeft[i][j - 1]);
			hv[i][j] -= dt / dy * (hvNetUpdatesAbove[i - 1][j - 1] + hvNetUpdatesBelow[i - 1][j]);

			if (h[i][j] < 0) {
				//TODO: dryTol
#ifndef NDEBUG
				// Only print this warning when debug is enabled
				// Otherwise we cannot vectorize this loop
				if (h[i][j] < -0.1) {
					std::cerr << "Warning, negative height: (i,j)=(" << i << "," << j << ")=" << h[i][j] << std::endl;
					std::cerr << "         b: " << b[i][j] << std::endl;
				}
#endif // NDEBUG
				//zero (small) negative depths
				h[i][j] = hu[i][j] = hv[i][j] = 0.;
			} else if (h[i][j] < 0.1)
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

bool
SWE_WavePropagationBlockSIMD::updateBathymetryWithDynamicDisplacement (scenarios::Asagi &i_asagiScenario, const float i_time)
{
	if (!i_asagiScenario.dynamicDisplacementAvailable (i_time))
		return false;

	// update the bathymetry
	for (int i = 0; i <= nx + 1; i++) {
		for (int j = 0; j <= ny + 1; j++) {
			b[i][j] = i_asagiScenario.getBathymetryAndDynamicDisplacement (
				offsetX + (i - 0.5f) * dx,
				offsetY + (j - 0.5f) * dy,
				i_time
			);
		}
	}

	setBoundaryBathymetry ();

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
void
SWE_WavePropagationBlockSIMD::simulateTimestep (float dt)
{
	computeNumericalFluxes ();
	updateUnknowns (dt);
}

/**
 * Runs the simulation until i_tEnd is reached.
 *
 * @param i_tStart time when the simulation starts
 * @param i_tEnd  time when the simulation should end
 * @return time we reached after the last update step, in general a bit later than i_tEnd
 */
float
SWE_WavePropagationBlockSIMD::simulate (float i_tStart, float i_tEnd)
{
	float t = i_tStart;
	do {
		//set values in ghost cells
		setGhostLayer ();

		// compute net updates for every edge
		computeNumericalFluxes ();
		//execute a wave propagation time step
		updateUnknowns (maxTimestep);
		t += maxTimestep;

		std::cout << "Simulation at time " << t << std::endl << std::flush;
	} while (t < i_tEnd);

	return t;
}
