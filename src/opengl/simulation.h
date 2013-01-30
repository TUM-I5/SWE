#ifndef SIMULATION_H
#define SIMULATION_H
// =====================================================================
// This file is part of SWE_CUDA (see file SWE_Block.cu for details).
// 
// Copyright (C) 2010,2011 Michael Bader, Kaveh Rahnema, Tobias Schnabel
// Copyright (C) 2012      Sebastian Rettenberger
// 
// SWE_CUDA is free software: you can redristribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// SWE_CUDA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with SWE_CUDA.  If not, see <http://www.gnu.org/licenses/>.
// =====================================================================
#include <math.h>
#include <cuda_runtime.h>

#include "blocks/cuda/SWE_BlockCUDA.hh"
#include "scenarios/SWE_simple_scenarios.hh"
#include "scenarios/SWE_VisInfo.hh"

void checkCUDAError(const char *msg);

class Simulation {

  public:
    // Constructor + Destructor 
    Simulation ();
    ~Simulation();

    // Restart simulation
    void restart();
    // Load new scenario after initialization
    void loadNewScenario(SWE_Scenario* scene);
    // Set a different resolution
    void resize(float factor);
    // Return the bathymetry data
    void setBathBuffer(float* output);
    // Simulate single timestep on graphics card
    void runCuda(struct cudaGraphicsResource **vbo_resource, struct cudaGraphicsResource **vbo_normals);

    int getNx() { return nx; }
    int getNy() { return ny; }

    const Float2D& getBathymetry() { return block->getBathymetry(); }

    void getScalingApproximation(float &bScale, float &bOffset, float &wScale);

  private:
    // Default scenario (used when no other scenario is specified)
    SWE_SplashingPoolScenario defaultScenario;

    // Instance of SWE_BlockCUDA 
    SWE_BlockCUDA* block;
    // Current scenario
    SWE_Scenario* scenario;

    // Cells in x dimension
    int nx;
    // Cells in y dimension
    int ny;

    // Current simulation time
    float curTime;
    // Is this our first simulation step?
    int isFirstStep;

    // Compute new water surface
    void calculateWaterSurface(float3* destBuffer);
    // Compute normals of the water surface for shading
    void calculateNormals(float3* vertexBuffer, float3* destBuffer);

    void updateVisBuffer(float3* _visBuffer);
    void debugVisBuffer(float3* _visBuffer);

    static void calculateNormal(float fVert1[], float fVert2[],
				float fVert3[], float fNormal[]);

};

#endif
