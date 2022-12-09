/**
 * @file
 * This file is part of SWE.
 *
 * @author Michael Bader (bader AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Univ.-Prof._Dr._MichaeBader)
 * @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de,
 * http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
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

#include <cmath>
#include <csignal>
#include <fenv.h>
#include <GASPI.h>
#include <GASPI_Ext.h>

#include "Blocks/Block.hpp"
#include "Blocks/WavePropagationBlock.hpp"
#include "Scenarios/BathymetryDamBreakScenario.hpp"
#include "Scenarios/RadialDamBreakScenario.hpp"
#include "Scenarios/SeaAtRestScenario.hpp"
#include "Scenarios/SplashingConeScenario.hpp"
#include "Scenarios/SplashingPoolScenario.hpp"
#include "Tools/Args.hpp"
#include "Tools/Logger.hpp"
#include "Tools/ProgressBar.hpp"
#include "Writers/Writer.hpp"

#ifndef _MSC_VER
#pragma RealType_control(precise, on)
#pragma STDC FENV_ACCESS ON
#endif

#ifdef _MSC_VER
void fpExceptionHandler(const int signal, const int nSubCode);
#endif

void success_or_exit (const char *file, const int line, const gaspi_return_t ec)
{
  if (ec != GASPI_SUCCESS)
  {
    Tools::Logger::logger.getDefaultOutputStream() << "Assertion failed in " << file 
                               << "[" << line << "]: Return " << ec 
                               << ": " <<  gaspi_error_str (ec) << std::endl;
    exit(EXIT_FAILURE);
  }
}

#define ASSERT(ec)  success_or_exit (__FILE__, __LINE__, ec)

/**
 * Computes the number of block rows from the total number of processes.
 *
 * The number of rows is determined as the square root of the
 * number of processes, if this is a square number;
 * otherwise, we use the largest number that is smaller than the square
 * root and still a divisor of the number of processes.
 *
 * @param numberOfProcesses number of processes
 * @return number of block rows
 */
int  computeNumberOfBlockRows(int numberOfProcesses);

void exchangeLeftRightGhostLayers( 
  const int           leftNeighborRank,
  const int           nleftInflowOffset, 
  gaspi_offset_t*     leftInflowOffset,
  const int           nleftOutflowOffset,
  gaspi_offset_t*     leftOutflowOffset,
  const int           rightNeighborRank,
  const int           nrightInflowOffset,
  gaspi_offset_t*     rightInflowOffset,
  int                 nrightOutflowOffset,
  gaspi_offset_t*     rightOutflowOffset,
  gaspi_segment_id_t* segment_id,
  gaspi_size_t*       size,
  gaspi_rank_t        gpiRank
);

void exchangeBottomTopGhostLayers( 
  const int           bottomNeighborRank,
  const int           nbottomInflowOffset, 
  gaspi_offset_t*     bottomInflowOffset,
  const int           nbottomOutflowOffset,
  gaspi_offset_t*     bottomOutflowOffset,
  const int           topNeighborRank,
  const int           ntopInflowOffset,
  gaspi_offset_t*     topInflowOffset,
  int                 ntopOutflowOffset,
  gaspi_offset_t*     topOutflowOffset,
  gaspi_segment_id_t* segment_id,
  gaspi_size_t*       size,
  gaspi_rank_t        gpiRank
);

gaspi_offset_t* calculateOffsets(Tools::Float2D<RealType>& grid, BoundaryEdge edge, BoundaryType type);

int main(int argc, char** argv) {
  //! GPI Rank of a process.
  gaspi_rank_t gpiRank;
  //! number of GPI processes.
  gaspi_rank_t numberOfProcesses;

  // initialize GPI
  gaspi_config_t config;
  ASSERT(gaspi_config_get(&config));
  config.queue_size_max = 1024;
  config.rw_list_elem_max = 1024;
  ASSERT(gaspi_config_set(config));
  ASSERT(gaspi_proc_init(GASPI_BLOCK));

  // determine local GPI rank
  gaspi_proc_rank(&gpiRank);
  // determine total number of processes
  gaspi_proc_num(&numberOfProcesses);

  Tools::Logger::logger.setProcessRank(gpiRank);
  Tools::Logger::logger.printWelcomeMessage();
  Tools::Logger::logger.printNumberOfProcesses(numberOfProcesses);

#ifndef _MSC_VER
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#else
  unsigned int oldState = 0;
  _controlfp_s(&oldState, _MCW_EM, _MCW_EM);
  const unsigned int flags      = _EM_INVALID | _EM_ZERODIVIDE | _EM_OVERFLOW;
  const unsigned int enableBits = flags & _MCW_EM;
  _clearfp();
  _controlfp_s(0, ~enableBits, enableBits);
  // std::signal(SIGFPE, (void(__cdecl*)(int)) fpExceptionHandler); // Gives better printing, but hides the actual
  // location of the RealTypeing point error.
#endif

  Tools::Args args;
  args.addOption("grid-size-x", 'x', "Number of cells in x direction");
  args.addOption("grid-size-y", 'y', "Number of cells in y direction");
  args.addOption("output-basepath", 'o', "Output base file name");
  args.addOption("number-of-checkpoints", 'n', "Number of checkpoints to write output files");

  Tools::Args::Result ret = args.parse(argc, argv, gpiRank == 0);

  switch (ret) {
  case Tools::Args::Result::Error:
    return EXIT_FAILURE;
  case Tools::Args::Result::Help:
    ASSERT(gaspi_proc_term(GASPI_BLOCK));
    return EXIT_SUCCESS;
  default:
    break;
  }

  int         numberOfGridCellsX  = args.getArgument<int>("grid-size-x", 16);
  int         numberOfGridCellsY  = args.getArgument<int>("grid-size-y", 16);
  std::string baseName            = args.getArgument<std::string>("output-basepath", "SWE");
  int         numberOfCheckPoints = args.getArgument<int>(
    "number-of-checkpoints", 20
  ); //! Number of checkpoints for visualization (at each checkpoint in time, an output file is written).

  // Print information about the grid
  Tools::Logger::logger.printNumberOfCells(numberOfGridCellsX, numberOfGridCellsY);

  // Determine the layout of GPI-ranks: use numberOfBlocksY*numberOfBlocksX grid blocks
  int numberOfBlocksY = computeNumberOfBlockRows(numberOfProcesses);
  int numberOfBlocksX = numberOfProcesses / numberOfBlocksY;
  Tools::Logger::logger.printNumberOfBlocks(numberOfBlocksX, numberOfBlocksY);

  // Determine local block coordinates of each block
  int blockPositionX = gpiRank / numberOfBlocksY;
  int blockPositionY = gpiRank % numberOfBlocksY;

  // Number of grid cells in x- and y-direction per process
  // Compute local number of cells for each block
  int nXLocal  = (blockPositionX < numberOfBlocksX - 1)
                   ? numberOfGridCellsX / numberOfBlocksX
                   : numberOfGridCellsX - (numberOfBlocksX - 1) * (numberOfGridCellsX / numberOfBlocksX);
  int nYLocal  = (blockPositionY < numberOfBlocksY - 1)
                   ? numberOfGridCellsY / numberOfBlocksY
                   : numberOfGridCellsY - (numberOfBlocksY - 1) * (numberOfGridCellsY / numberOfBlocksY);
  int nXNormal = numberOfGridCellsX / numberOfBlocksX;
  int nYNormal = numberOfGridCellsY / numberOfBlocksY;

  Tools::Logger::logger.printNumberOfCellsPerProcess(nXLocal, nYLocal);

  // Create a simple artificial scenario
  Scenarios::RadialDamBreakScenario scenario;

  // Compute the size of a single cell
  RealType cellSizeX = (scenario.getBoundaryPos(BoundaryEdge::Right) - scenario.getBoundaryPos(BoundaryEdge::Left))
                       / numberOfGridCellsX;
  RealType cellSizeY = (scenario.getBoundaryPos(BoundaryEdge::Top) - scenario.getBoundaryPos(BoundaryEdge::Bottom))
                       / numberOfGridCellsY;
  Tools::Logger::logger.printCellSize(cellSizeX, cellSizeY);

  Tools::Float2D<RealType> h(nXLocal+2, nYLocal+2);
  Tools::Float2D<RealType> hu(nXLocal+2, nYLocal+2);
  Tools::Float2D<RealType> hv(nXLocal+2, nYLocal+2);

  auto waveBlock = Blocks::Block::getBlockInstance(nXLocal, nYLocal, cellSizeX, cellSizeY, h, hu, hv);
  //Bind the block to a gaspi segment
  RealType *height = h.getData();
  RealType *discharge_hu = hu.getData();
  RealType *discharge_hv = hv.getData();

  ASSERT(gaspi_segment_use(0, (gaspi_pointer_t) height,
                           (nXLocal+2) * (nYLocal+2) * sizeof(RealType),
                           GASPI_GROUP_ALL, GASPI_BLOCK, 0));
  ASSERT(gaspi_segment_use(1, (gaspi_pointer_t) discharge_hu,
                           (nXLocal+2) * (nYLocal+2) * sizeof(RealType),
                           GASPI_GROUP_ALL, GASPI_BLOCK, 0));
  ASSERT(gaspi_segment_use(2, (gaspi_pointer_t) discharge_hv,
                           (nXLocal+2) * (nYLocal+2) * sizeof(RealType),
                           GASPI_GROUP_ALL, GASPI_BLOCK, 0));

  // Get the origin from the scenario
  RealType originX = scenario.getBoundaryPos(BoundaryEdge::Left) + blockPositionX * nXNormal * cellSizeX;
  RealType originY = scenario.getBoundaryPos(BoundaryEdge::Bottom) + blockPositionY * nYNormal * cellSizeY;

  // Initialise the wave propagation block
  waveBlock->initialiseScenario(originX, originY, scenario, true);

  // Get the final simulation time from the scenario
  double endSimulationTime = scenario.getEndSimulationTime();

  // Checkpoints when output files are written
  double* checkPoints = new double[numberOfCheckPoints + 1];

  // Compute the checkpoints in time
  for (int cp = 0; cp <= numberOfCheckPoints; cp++) {
    checkPoints[cp] = cp * (endSimulationTime / numberOfCheckPoints);
  }

  /*
   * Connect blocks at boundaries
   */
  // Left and right boundaries
  Tools::Logger::logger.printString("Connecting SWE blocks at left boundaries.");
  int nleftOutflowOffset = waveBlock->getWaterHeight().getRows();
  int nleftInflowOffset = nleftOutflowOffset;
  gaspi_offset_t *leftOutflowOffset = calculateOffsets(h, Left,Outflow);
  gaspi_offset_t *leftInflowOffset = calculateOffsets(h, Left, Inflow);
  if (blockPositionX == 0) {
    waveBlock->setBoundaryType(BoundaryEdge::Left, BoundaryType::Outflow);
  }

  Tools::Logger::logger.printString("Connecting SWE blocks at right boundaries.");
  int nrightOutflowOffset = nleftOutflowOffset;
  int nrightInflowOffset = nleftOutflowOffset;
  gaspi_offset_t *rightOutflowOffset = calculateOffsets(h, Right,Outflow);
  gaspi_offset_t *rightInflowOffset = calculateOffsets(h, Right, Inflow);
  if (blockPositionX == numberOfBlocksX - 1) {
    waveBlock->setBoundaryType(BoundaryEdge::Right, BoundaryType::Outflow);
  }

  // Bottom and top boundaries
  Tools::Logger::logger.printString("Connecting SWE blocks at bottom boundaries.");
  int nbottomOutflowOffset = waveBlock->getWaterHeight().getCols();
  int nbottomInflowOffset = nbottomOutflowOffset;
  gaspi_offset_t *bottomOutflowOffset = calculateOffsets(h, Bottom,Outflow);
  gaspi_offset_t *bottomInflowOffset = calculateOffsets(h, Bottom, Inflow);
  if (blockPositionY == 0) {
    waveBlock->setBoundaryType(BoundaryEdge::Bottom, BoundaryType::Outflow);
  }

  Tools::Logger::logger.printString("Connecting SWE blocks at top boundaries.");
  int ntopOutflowOffset = nbottomOutflowOffset;
  int ntopInflowOffset = nbottomOutflowOffset;
  gaspi_offset_t *topOutflowOffset = calculateOffsets(h, Top,Outflow);
  gaspi_offset_t *topInflowOffset = calculateOffsets(h, Top, Inflow);
  if (blockPositionY == numberOfBlocksY - 1) {
    waveBlock->setBoundaryType(BoundaryEdge::Top, BoundaryType::Outflow);
  }

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
   *  -> The stride for a row is ny+2, because we have to jump over a whole column
   *     for every row-element. This holds only in the CPU-version, in CUDA a buffer is implemented.
   *     See Blocks/CUDA/CUDABlock.hpp/.cu for details.
   *  -> The stride for a column is 1, because we can access the elements linear in memory.
   */

  // Compute GPI ranks of the neighbour processes
  int leftNeighborRank   = (blockPositionX > 0) ? gpiRank - numberOfBlocksY : -1;
  int rightNeighborRank  = (blockPositionX < numberOfBlocksX - 1) ? gpiRank + numberOfBlocksY : -1;
  int bottomNeighborRank = (blockPositionY > 0) ? gpiRank - 1 : -1;
  int topNeighborRank    = (blockPositionY < numberOfBlocksY - 1) ? gpiRank + 1 : -1;

  // Print the GPI grid
  Tools::Logger::logger.getDefaultOutputStream()
    << "Neighbors: " << leftNeighborRank << " (left), " << rightNeighborRank << " (right), " << bottomNeighborRank
    << " (bottom), " << topNeighborRank << " (top)" << std::endl;

  // Intially exchange ghost and copy layers
  gaspi_segment_id_t *segmentIdLr = NULL;
  gaspi_size_t *sizeLr = NULL;
  segmentIdLr = (gaspi_segment_id_t *) calloc (nleftOutflowOffset, sizeof(gaspi_segment_id_t));
  sizeLr = (gaspi_size_t *) calloc (nleftOutflowOffset, sizeof(gaspi_size_t));
  exchangeLeftRightGhostLayers(
    leftNeighborRank,
    nleftInflowOffset, 
    leftInflowOffset,
    nleftOutflowOffset,
    leftOutflowOffset,
    rightNeighborRank,
    nrightInflowOffset,
    rightInflowOffset,
    nrightOutflowOffset,
    rightOutflowOffset,
    segmentIdLr,
    sizeLr,
    gpiRank
  );

  gaspi_segment_id_t *segmentIdBt = NULL;
  gaspi_size_t *sizeBt = NULL;
  segmentIdBt = (gaspi_segment_id_t *) calloc (nbottomOutflowOffset, sizeof(gaspi_segment_id_t));
  sizeBt = (gaspi_size_t *) calloc (nbottomOutflowOffset, sizeof(gaspi_size_t));
  exchangeBottomTopGhostLayers(
    bottomNeighborRank, 
    nbottomInflowOffset,
    bottomInflowOffset,
    nbottomOutflowOffset,
    bottomOutflowOffset,
    topNeighborRank,
    ntopInflowOffset,
    topInflowOffset,
    ntopOutflowOffset,
    topOutflowOffset,
    segmentIdBt,
    sizeBt,
    gpiRank
  );
  ASSERT(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));

  Tools::ProgressBar progressBar(endSimulationTime, gpiRank);

  Tools::Logger::logger.printOutputTime(0.0);
  progressBar.update(0.0);

  // Boundary size of the ghost layers
  Writers::BoundarySize boundarySize = {{1, 1, 1, 1}};

  std::string fileName = Writers::generateBaseFileName(baseName, blockPositionX, blockPositionY);
  auto        writer   = Writers::Writer::createWriterInstance(
    fileName,
    waveBlock->getBathymetry(),
    boundarySize,
    nXLocal,
    nYLocal,
    cellSizeX,
    cellSizeY,
    blockPositionX * nXLocal,
    blockPositionY * nYLocal,
    originX,
    originY,
    0
  );

  // Write zero time step
  writer->writeTimeStep(waveBlock->getWaterHeight(), waveBlock->getDischargeHu(), waveBlock->getDischargeHv(), 0.0);

  // Print the start message and reset the wall clock time
  progressBar.clear();
  Tools::Logger::logger.printStartMessage();
  Tools::Logger::logger.initWallClockTime(time(NULL)); // GPI_Wtime()

  double simulationTime = 0.0;
  progressBar.update(simulationTime);

  unsigned int iterations = 0;

  // Loop over checkpoints
  for (int cp = 1; cp <= numberOfCheckPoints; cp++) {
    // Do time steps until next checkpoint is reached
    while (simulationTime < checkPoints[cp]) {
      // Reset CPU-Communication clock
      Tools::Logger::logger.resetClockToCurrentTime("CPU-Communication");

      // Exchange ghost and copy layers
      exchangeLeftRightGhostLayers(
        leftNeighborRank,
        nleftInflowOffset, 
        leftInflowOffset,
        nleftOutflowOffset,
        leftOutflowOffset,
        rightNeighborRank,
        nrightInflowOffset,
        rightInflowOffset,
        nrightOutflowOffset,
        rightOutflowOffset,
        segmentIdLr,
        sizeLr,
        gpiRank
      );
      exchangeBottomTopGhostLayers(
        bottomNeighborRank, 
        nbottomInflowOffset,
        bottomInflowOffset,
        nbottomOutflowOffset,
        bottomOutflowOffset,
        topNeighborRank,
        ntopInflowOffset,
        topInflowOffset,
        ntopOutflowOffset,
        topOutflowOffset,
        segmentIdBt,
        sizeBt,
        gpiRank
      );
      ASSERT(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));

      // Reset the cpu clock
      Tools::Logger::logger.resetClockToCurrentTime("CPU");

      // Set values in ghost cells
      waveBlock->setGhostLayer();

      // Compute numerical flux on each edge
      waveBlock->computeNumericalFluxes();

      // Approximate the maximum time step
      // waveBlock->computeMaxTimeStep();

      RealType maxTimeStepWidth = waveBlock->getMaxTimeStep();

      //! Maximum allowed time steps of all blocks
      RealType maxTimeStepWidthGlobal = RealType(0.0);

      // Determine smallest time step of all blocks
      ASSERT(gaspi_allreduce(&maxTimeStepWidth, &maxTimeStepWidthGlobal, 1, GASPI_OP_MIN, GASPI_TYPE_DOUBLE, GASPI_GROUP_ALL, GASPI_BLOCK));

      // Update the cell values
      waveBlock->updateUnknowns(maxTimeStepWidthGlobal);

      // Update the cpu time in the logger
      Tools::Logger::logger.updateTime("CPU");
      Tools::Logger::logger.updateTime("CPU-Communication");

      // Print the current simulation time
      progressBar.clear();
      Tools::Logger::logger.printSimulationTime(
        simulationTime,
        "[" + std::to_string(iterations) + "]: Simulation with max. global dt " + std::to_string(maxTimeStepWidthGlobal)
          + " at time"
      );

      // Update simulation time with time step width
      simulationTime += maxTimeStepWidthGlobal;
      iterations++;
      progressBar.update(simulationTime);
    }

    // Print current simulation time of the output
    progressBar.clear();
    Tools::Logger::logger.printOutputTime(simulationTime);
    progressBar.update(simulationTime);

    // Write output
    writer->writeTimeStep(
      waveBlock->getWaterHeight(), waveBlock->getDischargeHu(), waveBlock->getDischargeHv(), simulationTime
    );
  }

  progressBar.clear();
  Tools::Logger::logger.printStatisticsMessage();
  Tools::Logger::logger.printTime("CPU", "CPU Time");
  Tools::Logger::logger.printTime("CPU-Communication", "CPU + Communication Time");
  Tools::Logger::logger.printWallClockTime(time(NULL));
  Tools::Logger::logger.printIterationsDone(iterations);

  Tools::Logger::logger.printFinishMessage();

  delete waveBlock;
  delete[] checkPoints;

  ASSERT(gaspi_proc_term(GASPI_BLOCK));

  return EXIT_SUCCESS;
}

#ifdef _MSC_VER
void fpExceptionHandler(const int signal, const int nSubCode) {
  (void)signal;
  //_fpreset();
  _clearfp();
  switch (nSubCode) {
  case _FPE_INVALID:
    throw std::logic_error("Invalid number (NaN) encountered");
    break;
  case _FPE_DENORMAL:
    throw std::logic_error("Denormal");
    break;
  case _FPE_ZERODIVIDE:
    throw std::logic_error("Division by zero");
    break;
  case _FPE_OVERFLOW:
    throw std::logic_error("Overflow error encountered");
    break;
  case _FPE_UNDERFLOW:
    throw std::logic_error("Underflow error encountered");
    break;
  case _FPE_INEXACT:
    throw std::logic_error("Inexact RealTypeing point operation encountered");
    break;
  default:
    std::stringstream ss;
    ss << "RealTypeing point error with error code " << nSubCode;
    throw std::logic_error(ss.str());
    break;
  }
}
#endif

int computeNumberOfBlockRows(int numberOfProcesses) {
  int numberOfRows = static_cast<int>(std::sqrt(numberOfProcesses));
  while (numberOfProcesses % numberOfRows != 0)
    numberOfRows--;
  return numberOfRows;
}

void exchangeLeftRightGhostLayers( 
  const int           leftNeighborRank,
  const int           nleftInflowOffset, 
  gaspi_offset_t*     leftInflowOffset,
  const int           nleftOutflowOffset,
  gaspi_offset_t*     leftOutflowOffset,
  const int           rightNeighborRank,
  const int           nrightInflowOffset,
  gaspi_offset_t*     rightInflowOffset,
  int                 nrightOutflowOffset,
  gaspi_offset_t*     rightOutflowOffset,
  gaspi_segment_id_t* segment_id,
  gaspi_size_t*       size,
  gaspi_rank_t        gpiRank
) {
  gaspi_notification_id_t id = gpiRank, fid;
  gaspi_notification_t val = 1;
  // send to left, receive from the right:
  if (leftNeighborRank >= 0)
  {
    for(int j = 0; j < 3; j++)
    {
      for (int i = 0; i < nleftOutflowOffset; i++)
      {
        segment_id[i] = j;
        size[i] = sizeof(RealType);
      }
      val = 1;
      ASSERT(gaspi_write_list_notify(
        nleftOutflowOffset,
        segment_id,
        leftOutflowOffset,
        leftNeighborRank,
        segment_id,
        rightInflowOffset,
        size,
        j,
        id,
        val,
        0,
        GASPI_BLOCK)
      );
      ASSERT(gaspi_wait(0, GASPI_BLOCK));
    }
  }
  
  // send to right, receive from the left:
  if (rightNeighborRank >= 0)
  {
    for(int j = 0; j < 3; j++)
    {
      for (int i = 0; i < nrightOutflowOffset; i++)
      {
        segment_id[i] = j;
        size[i] = sizeof(RealType);
      }
      val = 1;
      ASSERT(gaspi_write_list_notify(
        nrightOutflowOffset,
        segment_id,
        rightOutflowOffset,
        rightNeighborRank,
        segment_id,
        leftInflowOffset,
        size,
        j,
        id,
        val,
        0,
        GASPI_BLOCK)
      );
      ASSERT(gaspi_wait(0, GASPI_BLOCK));
    }
  }
  if(rightNeighborRank >= 0)
  {
    /*Wait for a notification from right neighbor for its write on my right Inflow*/
    for(int j = 0; j < 3; j++)
    {
      ASSERT(gaspi_notify_waitsome(j, rightNeighborRank, 1, &fid, GASPI_BLOCK));
      val = 0;
      ASSERT (gaspi_notify_reset(j, fid, &val));
    }
  }
  if(leftNeighborRank >= 0)
  {
    /*Wait for a notification from right neighbor for its write on my right Inflow*/
    for(int j = 0; j < 3; j++)
    {
      ASSERT(gaspi_notify_waitsome(j, leftNeighborRank, 1, &fid, GASPI_BLOCK));
      val = 0;
      ASSERT (gaspi_notify_reset(j, fid, &val));
    }
  }
}

void exchangeBottomTopGhostLayers( 
  const int           bottomNeighborRank,
  const int           nbottomInflowOffset, 
  gaspi_offset_t*     bottomInflowOffset,
  const int           nbottomOutflowOffset,
  gaspi_offset_t*     bottomOutflowOffset,
  const int           topNeighborRank,
  const int           ntopInflowOffset,
  gaspi_offset_t*     topInflowOffset,
  int                 ntopOutflowOffset,
  gaspi_offset_t*     topOutflowOffset,
  gaspi_segment_id_t* segment_id,
  gaspi_size_t*       size,
  gaspi_rank_t        gpiRank
) {
  // send to bottom, receive from the top:
  gaspi_notification_id_t id = gpiRank, fid;
  gaspi_notification_t val = 1;
  if (bottomNeighborRank >= 0)
  {
    for(int j = 0; j < 3; j++)
    {
      for (int i = 0; i < nbottomOutflowOffset; i++)
      {
        segment_id[i] = j;
        size[i] = sizeof(RealType);
      }
      val = 1;
      ASSERT(gaspi_write_list_notify(
        nbottomOutflowOffset,
        segment_id,
        bottomOutflowOffset,
        bottomNeighborRank,
        segment_id,
        topInflowOffset,
        size,
        j,
        id,
        val,
        0,
        GASPI_BLOCK)
      );
      ASSERT(gaspi_wait(0, GASPI_BLOCK));
    }
  }
  
  // send to top, receive from the bottom:
  if (topNeighborRank >= 0)
  {
    for(int j = 0; j < 3; j++)
    {
      for (int i = 0; i < nbottomOutflowOffset; i++)
      {
        segment_id[i] = j;
        size[i] = sizeof(RealType);
      }
      val = 1;
      ASSERT(gaspi_write_list_notify(
        ntopOutflowOffset,
        segment_id,
        topOutflowOffset,
        topNeighborRank,
        segment_id,
        bottomInflowOffset,
        size,
        j,
        id,
        val,
        0,
        GASPI_BLOCK)
      );
      ASSERT(gaspi_wait(0, GASPI_BLOCK));
    }
  }
  if(topNeighborRank >= 0)
  {
    for(int j = 0; j < 3; j++)
    {
      ASSERT(gaspi_notify_waitsome(j, topNeighborRank, 1, &fid, GASPI_BLOCK));
      val = 0;
      ASSERT (gaspi_notify_reset(j, fid, &val));
    }
  }
  if(bottomNeighborRank >= 0)
  {
    for(int j = 0; j < 3; j++)
    {
      ASSERT(gaspi_notify_waitsome(j, bottomNeighborRank, 1, &fid, GASPI_BLOCK));
      val = 0;
      ASSERT (gaspi_notify_reset(j, fid, &val));
    }
  }
}

gaspi_offset_t *calculateOffsets(Tools::Float2D<RealType>& grid, BoundaryEdge edge, BoundaryType type)
{
  int rows = grid.getRows(); /*nY*/
  int cols = grid.getCols(); /*nX*/
  RealType *base = grid.getData();
  gaspi_offset_t *offset = NULL;
  
  switch(edge)
  {
    case Left:
      if (type ==Outflow)
      {
        RealType *first_ele = grid.getColProxy(1).getData();
        gaspi_offset_t start = first_ele - base;
        offset = (gaspi_offset_t *) calloc (rows, sizeof(gaspi_offset_t));
        for(int i = 0; i < rows; i++)
        {
          /*second row of grid*/
          offset[i] = (start + i) * sizeof(RealType);
        }
      }
      else if (type == Inflow)
      {
        RealType *first_ele = grid.getColProxy(0).getData();
        gaspi_offset_t start = first_ele - base;
        offset = (gaspi_offset_t *) calloc (rows, sizeof(gaspi_offset_t));
        for(int i = 0; i < rows; i++)
        {
          /*first row of grid*/
          offset[i] = (start + i) * sizeof(RealType);
          //std::cout << offset[i] << "\t"; 
        }
      }
      break;
    case Right:
      if (type ==Outflow)
      {
        RealType *first_ele = grid.getColProxy(cols-2).getData();
        gaspi_offset_t start = first_ele - base;
        offset = (gaspi_offset_t *) calloc (rows, sizeof(gaspi_offset_t));
        for(int i = 0; i < rows; i++)
        {
          /*second last row of grid*/
          offset[i] = (start + i) * sizeof(RealType); 
        }
      }
      else if (type == Inflow)
      {
        RealType *first_ele = grid.getColProxy(cols-1).getData();
        gaspi_offset_t start = first_ele - base;
        offset = (gaspi_offset_t *) calloc (rows, sizeof(gaspi_offset_t));
        for(int i = 0; i < rows; i++)
        {
          /*Last row of grid*/
          offset[i] = (start + i) * sizeof(RealType); 
        }
      }
      break;
    case Bottom:
      if (type ==Outflow)
      {
        RealType *first_ele = grid.getRowProxy(1).getData();
        gaspi_offset_t start = first_ele - base;
        offset = (gaspi_offset_t *) calloc (cols, sizeof(gaspi_offset_t));
        for(int i = 0; i < cols; i++)
        {
          /*second col of grid*/
          offset[i] = (start + (rows)*i) * sizeof(RealType);
        }
      }
      else if (type == Inflow)
      {
        RealType *first_ele = grid.getRowProxy(0).getData();
        gaspi_offset_t start = first_ele - base;
        offset = (gaspi_offset_t *) calloc (cols, sizeof(gaspi_offset_t));
        for(int i = 0; i < cols; i++)
        {
          /*first col of grid*/
          offset[i] = (start + (rows)*i) * sizeof(RealType);
        }
      }
      break;
    case Top:
      if (type ==Outflow)
      {
        RealType *first_ele = grid.getRowProxy(rows-2).getData();
        gaspi_offset_t start = first_ele - base;
        offset = (gaspi_offset_t *) calloc (cols, sizeof(gaspi_offset_t));
        for(int i = 0; i < cols; i++)
        {
          /*second last col of grid*/
          offset[i] = (start + (rows)*i) * sizeof(RealType); 
        }
      }
      else if (type == Inflow)
      {
        RealType *first_ele = grid.getRowProxy(rows-1).getData();
        gaspi_offset_t start = first_ele - base;
        offset = (gaspi_offset_t *) calloc (cols, sizeof(gaspi_offset_t));
        for(int i = 0; i < cols; i++)
        {
          /*Last col of grid*/
          offset[i] = (start + (rows)*i) * sizeof(RealType); 
        }
      }
      break;  
  }
  return offset;
}