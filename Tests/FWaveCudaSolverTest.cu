#if defined(ENABLE_CUDA)
#if defined(ENABLE_VECTORIZATION) || defined(ENABLE_VECTORIZATION_WITH_Cuda)
#error "CUDA should not be enabled with Vectorization. Please edit the ccmake .."
#endif

// FIXME[epic=SWE] cuda tests is disabled fix

#include <cmath>
#include <gtest/gtest.h>

#include "Blocks/Block.hpp"
#include "Blocks/WavePropagationBlock.cuh"
#include "Blocks/WavePropagationBlock.hpp"
#include "Scenarios/BathymetryDamBreakScenario.hpp"
#include "Scenarios/RadialDamBreakScenario.hpp"
#include "Scenarios/SeaAtRestScenario.hpp"
#include "Scenarios/SplashingConeScenario.hpp"
#include "Scenarios/SplashingPoolScenario.hpp"
#include "Solvers/FWaveCudaSolver.cuh"
#include "Solvers/FWaveSolver.hpp"

TEST(Test, FWaveCudaSolverTest)
{
    int mpiRank            = 0;
    int numberOfGridCellsX = 16;
    int numberOfGridCellsY = 16;

    int numberOfBlocksY = 1;
    int numberOfBlocksX = 1;

    // Determine local block coordinates of each block
    int blockPositionX = mpiRank / numberOfBlocksY;
    int blockPositionY = mpiRank % numberOfBlocksY;
    int nXLocal        = (blockPositionX < numberOfBlocksX - 1)
                             ? numberOfGridCellsX / numberOfBlocksX
                             : numberOfGridCellsX - (numberOfBlocksX - 1) * (numberOfGridCellsX / numberOfBlocksX);
    int nYLocal        = (blockPositionY < numberOfBlocksY - 1)
                             ? numberOfGridCellsY / numberOfBlocksY
                             : numberOfGridCellsY - (numberOfBlocksY - 1) * (numberOfGridCellsY / numberOfBlocksY);
    int nXNormal       = numberOfGridCellsX / numberOfBlocksX;
    int nYNormal       = numberOfGridCellsY / numberOfBlocksY;

    // Create a simple artificial scenario
#ifdef BATHYMETRYDAMBREAKSCENARIO
    Scenarios::BathymetryDamBreakScenario scenario;
#elif RADIALDAMBREAKSCENARIO
    Scenarios::RadialDamBreakScenario scenario;
#elif SEAATRESTSCENARIO
    Scenarios::SeaAtRestScenario scenario;
#elif SPLASHINGCONESCENARIO
    Scenarios::SplashingConeScenario scenario;
#elif SPLASHINGPOOLSCENARIO
    Scenarios::SplashingPoolScenario scenario;
#else
    Scenarios::RadialDamBreakScenario scenario;
#endif

    // Compute the size of a single cell
    RealType cellSizeX = (scenario.getBoundaryPos(BoundaryEdge::Right) - scenario.getBoundaryPos(BoundaryEdge::Left))
                         / numberOfGridCellsX;
    RealType cellSizeY = (scenario.getBoundaryPos(BoundaryEdge::Top) - scenario.getBoundaryPos(BoundaryEdge::Bottom))
                         / numberOfGridCellsY;

    Blocks::Block* waveBlockFWaveSolver = new Blocks::WavePropagationBlock(
        nXLocal, nYLocal, cellSizeX, cellSizeY, std::make_shared<Solvers::FWaveSolver<RealType>>()
    );
    Blocks::Block* waveBlockFWaveCudaSolver = new cuda::WavePropagationBlock(nXLocal, nYLocal, cellSizeX, cellSizeY);

    // Get the origin from the scenario
    RealType originX = scenario.getBoundaryPos(BoundaryEdge::Left) + blockPositionX * nXNormal * cellSizeX;
    RealType originY = scenario.getBoundaryPos(BoundaryEdge::Bottom) + blockPositionY * nYNormal * cellSizeY;

    // Initialise the wave propagation block
    waveBlockFWaveCudaSolver->initialiseScenario(originX, originY, scenario, true);
    waveBlockFWaveCudaSolver->setBoundaryType(BoundaryEdge::Left, BoundaryType::Outflow);
    waveBlockFWaveCudaSolver->setBoundaryType(BoundaryEdge::Right, BoundaryType::Outflow);
    waveBlockFWaveCudaSolver->setBoundaryType(BoundaryEdge::Bottom, BoundaryType::Outflow);
    waveBlockFWaveCudaSolver->setBoundaryType(BoundaryEdge::Top, BoundaryType::Outflow);

    // Initialise the wave propagation block
    waveBlockFWaveSolver->initialiseScenario(originX, originY, scenario, true);
    waveBlockFWaveSolver->setBoundaryType(BoundaryEdge::Left, BoundaryType::Outflow);
    waveBlockFWaveSolver->setBoundaryType(BoundaryEdge::Right, BoundaryType::Outflow);
    waveBlockFWaveSolver->setBoundaryType(BoundaryEdge::Bottom, BoundaryType::Outflow);
    waveBlockFWaveSolver->setBoundaryType(BoundaryEdge::Top, BoundaryType::Outflow);

    // Get the final simulation time from the scenario
    double endSimulationTime = scenario.getEndSimulationTime();
    double simulationTime    = 0.0;

    // Loop over checkpoints
    // Do time steps until next checkpoint is reached
    while (simulationTime < endSimulationTime)
    {
        // Set values in ghost cells
        waveBlockFWaveCudaSolver->setGhostLayer();

        // Compute numerical flux on each edge
        waveBlockFWaveCudaSolver->computeNumericalFluxes();

        RealType maxTimeStepWidth = waveBlockFWaveCudaSolver->getMaxTimeStep();

        // Update the cell values
        waveBlockFWaveCudaSolver->updateUnknowns(maxTimeStepWidth);

        LOG_DBG("[FWaveCudaSolver] Max TimeStep Width= " << maxTimeStepWidth);

        // Update simulation time with time step width
        simulationTime += maxTimeStepWidth;
    }

    simulationTime = 0.0;
    while (simulationTime < endSimulationTime)
    {
        // Set values in ghost cells
        waveBlockFWaveSolver->setGhostLayer();

        // Compute numerical flux on each edge
        waveBlockFWaveSolver->computeNumericalFluxes();

        RealType maxTimeStepWidth = waveBlockFWaveSolver->getMaxTimeStep();

        // Update the cell values
        waveBlockFWaveSolver->updateUnknowns(maxTimeStepWidth);

        LOG_DBG("[FWaveSolver] Max TimeStep Width= " << maxTimeStepWidth);

        // Update simulation time with time step width
        simulationTime += maxTimeStepWidth;
    }

    auto FWaveCudaSolver_WaterHeight = waveBlockFWaveCudaSolver->getWaterHeight();
    auto FWaveCudaSolver_DischargeHu = waveBlockFWaveCudaSolver->getDischargeHu();
    auto FWaveCudaSolver_DischargeHv = waveBlockFWaveCudaSolver->getDischargeHv();
    auto FWaveCudaSolver_Bathymetry  = waveBlockFWaveCudaSolver->getBathymetry();

    auto FWaveSolver_WaterHeight = waveBlockFWaveSolver->getWaterHeight();
    auto FWaveSolver_DischargeHu = waveBlockFWaveSolver->getDischargeHu();
    auto FWaveSolver_DischargeHv = waveBlockFWaveSolver->getDischargeHv();
    auto FWaveSolver_Bathymetry  = waveBlockFWaveSolver->getBathymetry();

    constexpr auto kFaultTolerance{1e-10};
    for (int i = 0; i < FWaveSolver_WaterHeight.getRows() * FWaveSolver_WaterHeight.getCols(); i++)
    {
        EXPECT_LE(FWaveSolver_WaterHeight.getData()[i] - FWaveCudaSolver_WaterHeight.getData()[i], kFaultTolerance);
        // if(FWaveSolver_WaterHeight.getData()[i] != FWaveCudaSolver_WaterHeight.getData()[i]) {
        // std::cout << std::setprecision(std::numeric_limits<double>::max_digits10) <<
        // FWaveSolver_WaterHeight.getData()[i]
        // << ", " << FWaveCudaSolver_WaterHeight.getData()[i] << std::endl;
        // }
    }
    LOG_INF << "FWaveSolver matches FWaveCudaSolver for: \t h";

    for (int i = 0; i < FWaveSolver_DischargeHu.getRows() * FWaveSolver_DischargeHu.getCols(); i++)
    {
        EXPECT_LE(FWaveSolver_DischargeHu.getData()[i] - FWaveCudaSolver_DischargeHu.getData()[i], kFaultTolerance);
    }
    LOG_INF << "FWaveSolver matches FWaveCudaSolver for: \t hu";

    for (int i = 0; i < FWaveSolver_DischargeHv.getRows() * FWaveSolver_DischargeHv.getCols(); i++)
    {
        EXPECT_LE(FWaveSolver_DischargeHv.getData()[i] - FWaveCudaSolver_DischargeHv.getData()[i], kFaultTolerance);
    }
    LOG_INF << "FWaveSolver matches FWaveCudaSolver for: \t hv";

    for (int i = 0; i < FWaveSolver_Bathymetry.getRows() * FWaveSolver_Bathymetry.getCols(); i++)
    {
        EXPECT_LE(FWaveSolver_Bathymetry.getData()[i] - FWaveCudaSolver_Bathymetry.getData()[i], kFaultTolerance);
    }
    LOG_INF << "FWaveSolver matches FWaveCudaSolver for: \t b";
}

#endif
