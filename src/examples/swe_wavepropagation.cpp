/**
 * @file
 * This file is part of SWE.
 *
 * @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *         Michael Bader (bader AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Univ.-Prof._Dr._Michael_Bader)
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
 * Very basic setting of SWE, which uses a wave propagation solver and an artificial scenario on a single block.
 */

#include "tools/help.hh"
#include <cstdlib>
#include <string>

#include "../SWE_Block.hh"
#include "../SWE_WavePropagationBlock.hh"
#include "../scenarios/SWE_simple_scenarios.h"
#include "../tools/Logger.hpp"

#ifndef STATICLOGGER
#define STATICLOGGER
#include "tools/Logger.hpp"
static tools::Logger s_sweLogger;
#endif

/**
 * Main program for the simulation on a single SWE_WavePropagationBlock.
 */
int main( int argc, char** argv ) {
  /**
   * Initialization.
   */
  // check if the necessary command line input parameters are given
  if(argc != 4) {
    std::cout << "Aborting ... please provide proper input parameters." << std::endl
              << "Example: ./SWE_parallel 200 300 /work/openmp_out" << std::endl
              << "\tfor a single block of size 200 * 300" << std::endl;
    assert(false);
  }

  //! number of grid cells in x- and y-direction.
  int l_nX, l_nY;

  //! l_baseName of the plots.
  std::string l_baseName;

  // read command line parameters
  l_nY = l_nX = atoi(argv[1]);
  l_nY = atoi(argv[2]);
  l_baseName = std::string(argv[3]);

  //create a simple artificial scenario
  SWE_BathymetryDamBreakScenario l_scenario;

  //! number of checkpoints for visualization (at each checkpoint in time, an output file is written).
  int l_numberOfCheckPoints = 40;

  //! size of a single cell in x- and y-direction
  float l_dX, l_dY;

  // compute the size of a single cell
  l_dX = (l_scenario.getBoundaryPos(BND_RIGHT) - l_scenario.getBoundaryPos(BND_LEFT) )/l_nX;
  l_dY = (l_scenario.getBoundaryPos(BND_TOP) - l_scenario.getBoundaryPos(BND_BOTTOM) )/l_nY;

  // initialize the grid data and the corresponding static variables
  SWE_Block::initGridData(l_nX,l_nY,l_dX,l_dY);

  //! origin of the simulation domain in x- and y-direction
  float l_originX, l_originY;

  // get the origin from the scenario
  l_originX = l_scenario.getBoundaryPos(BND_LEFT);
  l_originY = l_scenario.getBoundaryPos(BND_BOTTOM);

  // create a single wave propagation block
  SWE_WavePropagationBlock l_wavePropgationBlock(l_originX, l_originY);

  // initialize the wave propgation block
  l_wavePropgationBlock.initScenario(l_scenario);


  //! time when the simulation ends.
  float l_endSimulation = l_scenario.endSimulation();

  //! checkpoints when output files are written.
  float* l_checkPoints = new float[l_numberOfCheckPoints+1];

  // compute the checkpoints in time
  for(int cp = 0; cp <= l_numberOfCheckPoints; cp++) {
     l_checkPoints[cp] = cp*(l_endSimulation/l_numberOfCheckPoints);
  }

  // write the output at time zero
  l_wavePropgationBlock.writeVTKFileXML(generateFileName(l_baseName,0,0,0), l_nX, l_nY);


  /**
   * Simulation.
   */
  // print a start message and reset the wall clock time
  s_sweLogger.printStartMessage();
  s_sweLogger.initWallClockTime(time(NULL));

  //! simulation time.
  float l_t = 0.0;

  // loop over checkpoints
  for(int c=1; c<=l_numberOfCheckPoints; c++) {
    // reset the cpu clock
    s_sweLogger.resetCpuClockToCurrentTime();

    // do time steps until next checkpoint is reached
    while( l_t < l_checkPoints[c] ) {
      // set values in ghost cells:
      l_wavePropgationBlock.setGhostLayer();

      // compute numerical flux on each edge
      l_wavePropgationBlock.computeNumericalFluxes();

      //! maximum allowed time step width.
      float l_maxTimeStepWidth = l_wavePropgationBlock.getMaxTimestep();

      // update the cell values
      l_wavePropgationBlock.updateUnknowns(l_maxTimeStepWidth);

      // update simulation time with time step width.
      l_t += l_maxTimeStepWidth;

      // print the current simulation time
      s_sweLogger.printSimulationTime(l_t);
    }

    // update the cpu time in the logger
    s_sweLogger.updateCpuTime();

    // print current simulation time
    s_sweLogger.printOutputTime(l_t);

    // write vtk output
    l_wavePropgationBlock.writeVTKFileXML(generateFileName(l_baseName,c,0,0), l_nX, l_nY);
  }

  /**
   * Finalize.
   */
  // write the statistics message
  s_sweLogger.printStatisticsMessage();

  // print the cpu time
  s_sweLogger.printCpuTime("CPU time");

  // print the wall clock time (includes plotting)
  s_sweLogger.printWallClockTime(time(NULL));

  return 0;
}
