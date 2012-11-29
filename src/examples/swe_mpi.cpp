/**
 * @file
 * This file is part of SWE.
 *
 * @author Michael Bader (bader AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Univ.-Prof._Dr._Michael_Bader)
 *         Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
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
 * Setting of SWE, which uses a wave propagation solver and an artificial or ASAGI scenario on multiple blocks.
 */

#include <mpi.h>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <string>

#include "../tools/help.hh"

#include "../SWE_Block.hh"

#ifndef CUDA
#include "../SWE_WavePropagationBlock.hh"
#else
#include "../SWE_WavePropagationBlockCuda.hh"
#endif

#ifdef WRITENETCDF
#include "writer/NetCdfWriter.hh"
#else
#include "writer/VtkWriter.hh"
#endif

#ifdef ASAGI
#include "../scenarios/SWE_AsagiScenario.hpp"
#else
#include "../scenarios/SWE_simple_scenarios.h"
#endif

#ifdef READXML
#include "../tools/CXMLConfig.hpp"
#endif


#include "../tools/Logger.hpp"
#include "tools/ProgressBar.hh"

/**
 * Compute the number of block rows from the total number of processes.
 *
 * The number of rows is determined as the square root of the
 * number of processes, if this is a square number;
 * otherwise, we use the largest number that is smaller than the square
 * root and still a divisor of the number of processes.
 *
 * @param numProcs number of process.
 * @return number of block rows
 */
int computeNumberOfBlockRows(int i_numberOfProcesses) {
  int l_numberOfRows = std::sqrt(i_numberOfProcesses);
  while (i_numberOfProcesses % l_numberOfRows != 0) l_numberOfRows--;
  return l_numberOfRows;
};

// Exchanges the left and right ghost layers.
void exchangeLeftRightGhostLayers( const int i_leftNeighborRank,  SWE_Block1D* o_leftInflow,  SWE_Block1D* i_leftOutflow,
                                   const int i_rightNeighborRank, SWE_Block1D* o_rightInflow, SWE_Block1D* i_rightOutflow,
                                   MPI_Datatype i_mpiCol);

// Exchanges the bottom and top ghist layers.
void exchangeBottomTopGhostLayers( const int i_bottomNeighborRank, SWE_Block1D* o_bottomNeighborInflow, SWE_Block1D* i_bottomNeighborOutflow,
                                   const int i_topNeighborRank,    SWE_Block1D* o_topNeighborInflow,    SWE_Block1D* i_topNeighborOutflow,
                                   const MPI_Datatype i_mpiRow);

/**
 * Main program for the simulation on a single SWE_WavePropagationBlock.
 */
int main( int argc, char** argv ) {
  /**
   * Initialization.
   */
  //! MPI Rank of a process.
  int l_mpiRank;
  //! number of MPI processes.
  int l_numberOfProcesses;

  // initialize MPI
  if ( MPI_Init(&argc,&argv) != MPI_SUCCESS ) {
    std::cerr << "MPI_Init failed." << std::endl;
  }

  // determine local MPI rank
  MPI_Comm_rank(MPI_COMM_WORLD,&l_mpiRank);
  // determine total number of processes
  MPI_Comm_size(MPI_COMM_WORLD,&l_numberOfProcesses);

  // initialize a logger for every MPI process
  tools::Logger::logger.setProcessRank(l_mpiRank);

  // print the welcome message
  tools::Logger::logger.printWelcomeMessage();

  // set current wall clock time within the solver
  tools::Logger::logger.initWallClockTime( MPI_Wtime() );
  //print the number of processes
  tools::Logger::logger.printNumberOfProcesses(l_numberOfProcesses);

  // check if the necessary command line input parameters are given
  #ifndef READXML
  if(argc != 4) {
    std::cout << "Aborting ... please provide proper input parameters." << std::endl
              << "Example: ./SWE_parallel 200 300 /work/openmp_out" << std::endl
              << "\tfor a multiple block grid of total size 200*300" << std::endl << std::flush;

    // abort MPI execution
    MPI_Abort(MPI_COMM_WORLD, -1);

    return 1;
  }
  #endif

  //! total number of grid cell in x- and y-direction.
  int l_nX, l_nY;

  //! l_baseName of the plots.
  std::string l_baseName;

  // read command line parameters
  #ifndef READXML
  l_nX = atoi(argv[1]);
  l_nY = atoi(argv[2]);

  l_baseName = std::string(argv[3]);
  #endif

  // read xml file
  #ifdef READXML
  assert(false); //TODO: not implemented.
  if(argc != 2) {
    l_sweLogger.printString("Aborting. Please provide a proper input file.");
    l_sweLogger.printString("Example: ./SWE_gnu_debug_none_augrie config.xml");
    return 1;
  }
  l_sweLogger.printString("Reading xml-file.");

  std::string l_xmlFile = std::string(argv[1]);
  l_sweLogger.printString(l_xmlFile);

  CXMLConfig l_xmlConfig;
  l_xmlConfig.loadConfig(l_xmlFile.c_str());
  #endif

  //! number of SWE_Blocks in x- and y-direction.
  int l_blocksX, l_blocksY;

  // determine the layout of MPI-ranks: use l_blocksX*l_blocksY grid blocks
  l_blocksY = computeNumberOfBlockRows(l_numberOfProcesses);
  l_blocksX = l_numberOfProcesses/l_blocksY;

  // print information about the grid
  tools::Logger::logger.printNumberOfCells(l_nX, l_nY);
  tools::Logger::logger.printNumberOfBlocks(l_blocksX, l_blocksY);


  //! local position of each MPI process in x- and y-direction.
  int l_blockPositionX, l_blockPositionY;

  // determine local block coordinates of each SWE_Block
  l_blockPositionX = l_mpiRank / l_blocksY;
  l_blockPositionY = l_mpiRank % l_blocksY;

  #ifdef ASAGI
  /*
   * Pixel node registration used [Cartesian grid]
   * Grid file format: nf = GMT netCDF format (float)  (COARDS-compliant)
   * x_min: -500000 x_max: 6500000 x_inc: 500 name: x nx: 14000
   * y_min: -2500000 y_max: 1500000 y_inc: 500 name: y ny: 8000
   * z_min: -6.48760175705 z_max: 16.1780223846 name: z
   * scale_factor: 1 add_offset: 0
   * mean: 0.00217145586762 stdev: 0.245563641735 rms: 0.245573241263
   */

  //simulation area
  float simulationArea[4];
  simulationArea[0] = -450000;
  simulationArea[1] = 6450000;
  simulationArea[2] = -2450000;
  simulationArea[3] = 1450000;

  SWE_AsagiScenario l_scenario( ASAGI_INPUT_DIR "tohoku_gebco_ucsb3_500m_hawaii_bath.nc",
		  	  	  	  	  	  	ASAGI_INPUT_DIR "tohoku_gebco_ucsb3_500m_hawaii_displ.nc",
                                (float) 14400., simulationArea);
  #else
  // create a simple artificial scenario
  SWE_BathymetryDamBreakScenario l_scenario;
  #endif

  //! number of checkpoints for visualization (at each checkpoint in time, an output file is written).
  int l_numberOfCheckPoints = 40;


  //! number of grid cells in x- and y-direction per process.
  int l_nXLocal, l_nYLocal;

  //! size of a single cell in x- and y-direction
  float l_dX, l_dY;

  // compute local number of cells for each SWE_Block
  l_nXLocal = (l_blockPositionX < l_blocksX-1) ? l_nX/l_blocksX : l_nX - (l_blocksX-1)*(l_nX/l_blocksX);
  l_nYLocal = (l_blockPositionY < l_blocksY-1) ? l_nY/l_blocksY : l_nY - (l_blocksY-1)*(l_nY/l_blocksY);

  // compute the size of a single cell
  l_dX = (l_scenario.getBoundaryPos(BND_RIGHT) - l_scenario.getBoundaryPos(BND_LEFT) )/l_nX;
  l_dY = (l_scenario.getBoundaryPos(BND_TOP) - l_scenario.getBoundaryPos(BND_BOTTOM) )/l_nY;

  // print information about the cell size and local number of cells
  tools::Logger::logger.printCellSize(l_dX, l_dY);
  tools::Logger::logger.printNumberOfCellsPerProcess(l_nXLocal, l_nYLocal);

  // initialize the grid data and the corresponding static variables
  SWE_Block::initGridData(l_nXLocal,l_nYLocal,l_dX,l_dY);

  //! origin of the simulation domain in x- and y-direction
  float l_originX, l_originY;

  // get the origin from the scenario
  l_originX = l_scenario.getBoundaryPos(BND_LEFT) + l_blockPositionX*l_nXLocal*l_dX;;
  l_originY = l_scenario.getBoundaryPos(BND_BOTTOM) + l_blockPositionY*l_nYLocal*l_dY;

  // create a single wave propagation block
  #ifndef CUDA
  SWE_WavePropagationBlock l_wavePropgationBlock;
  #else
  //! number of CUDA devices per node TODO: hardcoded
  int l_cudaDevicesPerNode = 7;

  //! the id of the node local GPU
  int l_cudaDeviceId = l_mpiRank % l_cudaDevicesPerNode;

  SWE_WavePropagationBlockCuda l_wavePropgationBlock(l_cudaDeviceId);
  #endif

  // initialize the wave propgation block
  l_wavePropgationBlock.initScenario(l_originX, l_originY, l_scenario, true);

  //! time when the simulation ends.
  float l_endSimulation = l_scenario.endSimulation();

  //! checkpoints when output files are written.
  float* l_checkPoints = new float[l_numberOfCheckPoints+1];

  // compute the checkpoints in time
  for(int cp = 0; cp <= l_numberOfCheckPoints; cp++) {
     l_checkPoints[cp] = cp*(l_endSimulation/l_numberOfCheckPoints);
  }

  /*
   * Connect SWE blocks at boundaries
   */
  // left and right boundaries
  tools::Logger::logger.printString("Connecting SWE blocks at left boundaries.");
  SWE_Block1D* l_leftInflow  = l_wavePropgationBlock.grabGhostLayer(BND_LEFT);
  SWE_Block1D* l_leftOutflow = l_wavePropgationBlock.registerCopyLayer(BND_LEFT);
  if (l_blockPositionX == 0)
    l_wavePropgationBlock.setBoundaryType(BND_LEFT, OUTFLOW);

  tools::Logger::logger.printString("Connecting SWE blocks at right boundaries.");
  SWE_Block1D* l_rightInflow  = l_wavePropgationBlock.grabGhostLayer(BND_RIGHT);
  SWE_Block1D* l_rightOutflow = l_wavePropgationBlock.registerCopyLayer(BND_RIGHT);
  if (l_blockPositionX == l_blocksX-1)
    l_wavePropgationBlock.setBoundaryType(BND_RIGHT, OUTFLOW);

  // bottom and top boundaries
  tools::Logger::logger.printString("Connecting SWE blocks at bottom boundaries.");
  SWE_Block1D* l_bottomInflow  = l_wavePropgationBlock.grabGhostLayer(BND_BOTTOM);
  SWE_Block1D* l_bottomOutflow = l_wavePropgationBlock.registerCopyLayer(BND_BOTTOM);
  if (l_blockPositionY == 0)
    l_wavePropgationBlock.setBoundaryType(BND_BOTTOM, OUTFLOW);

  tools::Logger::logger.printString("Connecting SWE blocks at top boundaries.");
  SWE_Block1D* l_topInflow  = l_wavePropgationBlock.grabGhostLayer(BND_TOP);
  SWE_Block1D* l_topOutflow = l_wavePropgationBlock.registerCopyLayer(BND_TOP);
  if (l_blockPositionY == l_blocksY-1)
    l_wavePropgationBlock.setBoundaryType(BND_TOP, OUTFLOW);

  /*
   * The grid is stored column wise in memory:
   *
   *        ************************** . . . **********
   *        *       *  ny+2 *2(ny+2)*         * (ny+1)*
   *        *  ny+1 * +ny+1 * +ny+1 *         * (ny+2)*
   *        *       *       *       *         * +ny+1 *
   *        ************************** . . . **********
   *        *       *       *       *         *       *
   *        .       .       .       .         .       .
   *        .       .       .       .         .       .
   *        .       .       .       .         .       .
   *        *       *       *       *         *       *
   *        ************************** . . . **********
   *        *       *  ny+2 *2(ny+2)*         * (ny+1)*
   *        *   1   *   +1  *   +1  *         * (ny+2)*
   *        *       *       *       *         *   +1  *
   *        ************************** . . . **********
   *        *       *  ny+2 *2(ny+2)*         * (ny+1)*
   *        *   0   *   +0  *   +0  *         * (ny+2)*
   *        *       *       *       *         *   +0  *
   *        ************************** . . . ***********
   *
   *
   *  -> The stride for a row is ny+2, because we have to jump over a whole column
   *     for every row-element. This holds only in the CPU-version, in CUDA a buffer is implemented.
   *     See SWE_BlockCUDA.hh/.cu for details.
   *  -> The stride for a column is 1, because we can access the elements linear in memory.
   */
  //! MPI row-vector: l_nXLocal+2 blocks, 1 element per block, stride of l_nYLocal+2
  MPI_Datatype l_mpiRow;
  #ifndef CUDA
  MPI_Type_vector(l_nXLocal+2, 1          , l_nYLocal+2, MPI_FLOAT, &l_mpiRow);
  #else
  MPI_Type_vector(1,           l_nXLocal+2, 1          , MPI_FLOAT, &l_mpiRow);
  #endif
  MPI_Type_commit(&l_mpiRow);

  //! MPI row-vector: 1 block, l_nYLocal+2 elements per block, stride of 1
  MPI_Datatype l_mpiCol;
  MPI_Type_vector(1,           l_nYLocal+2, 1,           MPI_FLOAT, &l_mpiCol);
  MPI_Type_commit(&l_mpiCol);

  //! MPI ranks of the neighbors
  int l_leftNeighborRank, l_rightNeighborRank, l_bottomNeighborRank, l_topNeighborRank;

  // compute MPI ranks of the neighbour processes
  l_leftNeighborRank   = (l_blockPositionX > 0) ? l_mpiRank-l_blocksY : MPI_PROC_NULL;
  l_rightNeighborRank  = (l_blockPositionX < l_blocksX-1) ? l_mpiRank+l_blocksY : MPI_PROC_NULL;
  l_bottomNeighborRank = (l_blockPositionY > 0) ? l_mpiRank-1 : MPI_PROC_NULL;
  l_topNeighborRank    = (l_blockPositionY < l_blocksY-1) ? l_mpiRank+1 : MPI_PROC_NULL;

  // print the MPI grid
  tools::Logger::logger.cout() << "neighbors: "
                     << l_leftNeighborRank << " (left), "
                     << l_rightNeighborRank << " (right), "
                     << l_bottomNeighborRank << " (bottom), "
                     << l_topNeighborRank << " (top)" << std::endl;

  // intially exchange ghost and copy layers
  exchangeLeftRightGhostLayers( l_leftNeighborRank,  l_leftInflow,  l_leftOutflow,
                  l_rightNeighborRank, l_rightInflow, l_rightOutflow,
                  l_mpiCol );

  exchangeBottomTopGhostLayers( l_bottomNeighborRank, l_bottomInflow, l_bottomOutflow,
                  l_topNeighborRank,    l_topInflow,    l_topOutflow,
                  l_mpiRow );

  // Init fancy progressbar
  tools::ProgressBar progressBar(l_endSimulation, l_mpiRank);

  // write the output at time zero
  tools::Logger::logger.printOutputTime(0);
  progressBar.update(0.);

  std::string l_fileName = generateBaseFileName(l_baseName,l_blockPositionX,l_blockPositionY);
  //boundary size of the ghost layers
  io::BoundarySize l_boundarySize = {1};
#ifdef WRITENETCDF
  //construct a NetCdfWriter
  io::NetCdfWriter l_writer( l_fileName,
		  l_wavePropgationBlock.getBathymetry(),
		  l_boundarySize,
		  l_nXLocal, l_nYLocal,
		  l_dX, l_dY,
          l_originX, l_originY );
#else
  // Construct a VtkWriter
  io::VtkWriter l_writer( l_fileName,
		  l_wavePropgationBlock.getBathymetry(),
		  l_boundarySize,
		  l_nXLocal, l_nYLocal,
		  l_dX, l_dY,
		  l_blockPositionX*l_nXLocal, l_blockPositionY*l_nYLocal );
#endif
  // Write zero time step
  l_writer.writeTimeStep( l_wavePropgationBlock.getWaterHeight(),
                          l_wavePropgationBlock.getDischarge_hu(),
                          l_wavePropgationBlock.getDischarge_hv(),
                          (float) 0.);
  /**
   * Simulation.
   */
  // print the start message and reset the wall clock time
  progressBar.clear();
  tools::Logger::logger.printStartMessage();
  tools::Logger::logger.initWallClockTime(time(NULL));

  //! simulation time.
  float l_t = 0.0;
  progressBar.update(l_t);

  // loop over checkpoints
  for(int c=1; c<=l_numberOfCheckPoints; c++) {
    //reset CPU-Communication clock
	  tools::Logger::logger.resetCpuCommunicationClockToCurrentTime();

    // do time steps until next checkpoint is reached
    while( l_t < l_checkPoints[c] ) {
      // exchange ghost and copy layers
      exchangeLeftRightGhostLayers( l_leftNeighborRank,  l_leftInflow,  l_leftOutflow,
                      l_rightNeighborRank, l_rightInflow, l_rightOutflow,
                      l_mpiCol );

      exchangeBottomTopGhostLayers( l_bottomNeighborRank, l_bottomInflow, l_bottomOutflow,
                      l_topNeighborRank,    l_topInflow,    l_topOutflow,
                      l_mpiRow );

      // reset the cpu clock
      tools::Logger::logger.resetCpuClockToCurrentTime();

      // set values in ghost cells
      l_wavePropgationBlock.setGhostLayer();

      // compute numerical flux on each edge
      l_wavePropgationBlock.computeNumericalFluxes();

      //! maximum allowed time step width within a block.
      float l_maxTimeStepWidth = l_wavePropgationBlock.getMaxTimestep();

      //! maximum allowed time steps of all blocks
      float l_maxTimeStepWidthGlobal;

      // determine smallest time step of all blocks
      MPI_Allreduce(&l_maxTimeStepWidth, &l_maxTimeStepWidthGlobal, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);

      // update the cell values
      l_wavePropgationBlock.updateUnknowns(l_maxTimeStepWidthGlobal);

      // update simulation time with time step width.
      l_t += l_maxTimeStepWidthGlobal;

      // print the current simulation time
      progressBar.clear();
      tools::Logger::logger.printSimulationTime(l_t);
      progressBar.update(l_t);
    }

    // update and reset the cpu time in the logger
    tools::Logger::logger.updateCpuTime();
    tools::Logger::logger.resetCpuClockToCurrentTime();

    // print current simulation time
    progressBar.clear();
    tools::Logger::logger.printOutputTime(l_t);
    progressBar.update(l_t);

    // write output
    l_writer.writeTimeStep( l_wavePropgationBlock.getWaterHeight(),
                            l_wavePropgationBlock.getDischarge_hu(),
                            l_wavePropgationBlock.getDischarge_hv(),
                            l_t);
  }

  /**
   * Finalize.
   */
#ifdef ASAGI
  // Free ASAGI resources
  l_scenario.deleteGrids();
#endif

  progressBar.clear();

  // write the statistics message
  tools::Logger::logger.printStatisticsMessage();

  // print the cpu time
  tools::Logger::logger.printCpuTime("CPU time");

  // print the wall clock time (includes plotting)
  tools::Logger::logger.printWallClockTime(time(NULL));

  // print the finish message
  tools::Logger::logger.printFinishMessage();

  // finalize MPI execution
  MPI_Finalize();

  return 0;
}


/**
 * Exchanges the left and right ghost layers with MPI's SendReceive.
 *
 * @param i_leftNeighborRank MPI rank of the  left neighbor.
 * @param o_leftInflow ghost layer, where the left neighbor writes into.
 * @param i_leftOutflow layer where the left neighbor reads from.
 * @param i_rightNeighborRank MPI rank of the right neighbor.
 * @param o_rightInflow ghost layer, where the right neighbor writes into.
 * @param i_rightOutflow layer, where the right neighbor reads form.
 * @param i_mpiCol MPI data type for the vertical gost layers.
 */
void exchangeLeftRightGhostLayers( const int i_leftNeighborRank,  SWE_Block1D* o_leftInflow,  SWE_Block1D* i_leftOutflow,
                                   const int i_rightNeighborRank, SWE_Block1D* o_rightInflow, SWE_Block1D* i_rightOutflow,
                                   MPI_Datatype i_mpiCol) {

  MPI_Status l_status;

  // send to left, receive from the right:
  MPI_Sendrecv( i_leftOutflow->h.elemVector(), 1, i_mpiCol, i_leftNeighborRank,  1,
                o_rightInflow->h.elemVector(), 1, i_mpiCol, i_rightNeighborRank, 1,
                MPI_COMM_WORLD, &l_status );

  MPI_Sendrecv( i_leftOutflow->hu.elemVector(), 1, i_mpiCol, i_leftNeighborRank,  2,
                o_rightInflow->hu.elemVector(), 1, i_mpiCol, i_rightNeighborRank, 2,
                MPI_COMM_WORLD, &l_status );

  MPI_Sendrecv( i_leftOutflow->hv.elemVector(), 1, i_mpiCol, i_leftNeighborRank,  3,
                o_rightInflow->hv.elemVector(), 1, i_mpiCol, i_rightNeighborRank, 3,
                MPI_COMM_WORLD, &l_status );

  // send to right, receive from the left:
  MPI_Sendrecv( i_rightOutflow->h.elemVector(), 1, i_mpiCol, i_rightNeighborRank, 4,
                o_leftInflow->h.elemVector(),   1, i_mpiCol, i_leftNeighborRank,  4,
                MPI_COMM_WORLD, &l_status );

  MPI_Sendrecv( i_rightOutflow->hu.elemVector(), 1, i_mpiCol, i_rightNeighborRank, 5,
                o_leftInflow->hu.elemVector(),   1, i_mpiCol, i_leftNeighborRank,  5,
                MPI_COMM_WORLD, &l_status);

  MPI_Sendrecv( i_rightOutflow->hv.elemVector(), 1, i_mpiCol, i_rightNeighborRank, 6,
                o_leftInflow->hv.elemVector(),   1, i_mpiCol, i_leftNeighborRank,  6,
                MPI_COMM_WORLD, &l_status );

}

/**
 * Exchanges the bottom and top ghost layers with MPI's SendReceive.
 *
 * @param i_bottomNeighborRank MPI rank of the bottom neighbor.
 * @param o_bottomNeighborInflow ghost layer, where the bottom neighbor writes into.
 * @param i_bottomNeighborOutflow host layer, where the bottom neighbor reads from.
 * @param i_topNeighborRank MPI rank of the top neighbor.
 * @param o_topNeighborInflow ghost layer, where the top neighbor writes into.
 * @param i_topNeighborOutflow ghost layer, where the top neighbor reads from.
 * @param i_mpiRow MPI data type for the horizontal ghost layers.
 */
void exchangeBottomTopGhostLayers( const int i_bottomNeighborRank, SWE_Block1D* o_bottomNeighborInflow, SWE_Block1D* i_bottomNeighborOutflow,
                                   const int i_topNeighborRank,    SWE_Block1D* o_topNeighborInflow,    SWE_Block1D* i_topNeighborOutflow,
                                   const MPI_Datatype i_mpiRow) {
  MPI_Status l_status;

  // send to bottom, receive from the top:
  MPI_Sendrecv( i_bottomNeighborOutflow->h.elemVector(), 1, i_mpiRow, i_bottomNeighborRank, 11,
                o_topNeighborInflow->h.elemVector(),     1, i_mpiRow, i_topNeighborRank,11,
                MPI_COMM_WORLD, &l_status );

  MPI_Sendrecv( i_bottomNeighborOutflow->hu.elemVector(), 1, i_mpiRow, i_bottomNeighborRank, 12,
                o_topNeighborInflow->hu.elemVector(),     1, i_mpiRow, i_topNeighborRank,    12,
                MPI_COMM_WORLD, &l_status );

  MPI_Sendrecv( i_bottomNeighborOutflow->hv.elemVector(), 1, i_mpiRow, i_bottomNeighborRank, 13,
                o_topNeighborInflow->hv.elemVector(),     1, i_mpiRow, i_topNeighborRank, 13,
                MPI_COMM_WORLD, &l_status);

  // send to top, receive from the bottom:
  MPI_Sendrecv( i_topNeighborOutflow->h.elemVector(),   1, i_mpiRow, i_topNeighborRank,    14,
                o_bottomNeighborInflow->h.elemVector(), 1, i_mpiRow, i_bottomNeighborRank, 14,
                MPI_COMM_WORLD, &l_status );

  MPI_Sendrecv( i_topNeighborOutflow->hu.elemVector(),   1, i_mpiRow, i_topNeighborRank, 15,
                o_bottomNeighborInflow->hu.elemVector(), 1, i_mpiRow, i_bottomNeighborRank, 15,
                MPI_COMM_WORLD, &l_status );

  MPI_Sendrecv( i_topNeighborOutflow->hv.elemVector(),   1, i_mpiRow, i_topNeighborRank,    16,
                o_bottomNeighborInflow->hv.elemVector(), 1, i_mpiRow, i_bottomNeighborRank, 16,
                MPI_COMM_WORLD, &l_status );

}
