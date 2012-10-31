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

#include "../SWE_BlockCUDA.hh"
#include "../scenarios/SWE_Scenario.h"
#include "../scenarios/SWE_VisInfo.hh"

void checkCUDAError(const char *msg);

class Simulation {

  public:
    // Constructor + Destructor 
    Simulation (int nx, int ny, float dx, float dy, 
                SWE_Scenario* scene, SWE_BlockCUDA* init_splash);
    ~Simulation();

    // Restart simulation
    void restart();
    // Load new scenario after initialization
    void loadNewScenario(SWE_Scenario* scene);
    // Save simulation state to file
    void saveToFile();
    // Return the bathymetry data
    void setBathBuffer(float* output);
    // Simulate single timestep on graphics card
    void runCuda(struct cudaGraphicsResource **vbo_resource, struct cudaGraphicsResource **vbo_normals);

    // Debugging
    void writeDebugOutput(float3* destBuffer = NULL);

  protected:
  public:

    // Instance of SWE_BlockCUDA 
    SWE_BlockCUDA* splash;
    // Current scenario
    SWE_Scenario* myScenario;

    // Current simulation time
    float curTime;
    // Is this our first simulation step?
    int isFirstStep;
    // Maximum of cell sizes
    float maxCellSize;
    // Initialize boundaries defined by the scene
    void initBoundaries(SWE_Scenario* scene);
    // Initialize boundaries defined by an input file

    int fileNumber;
    // Use file input as boundary data
    bool useFileInput;
    // Holds maximum dimension
    int maxDim;
    // Compute new water surface
    void calculateWaterSurface(float3* destBuffer);
    // Compute normals of the water surface for shading
    void calculateNormals(float3* vertexBuffer, float3* destBuffer);

    void getScalingApproximation(float &bScale, float &bOffset, float &wScale);

    void updateVisBuffer(float3* _visBuffer);
    void debugVisBuffer(float3* _visBuffer);

    static void calculateNormal(float fVert1[], float fVert2[],
				float fVert3[], float fNormal[]);

};

#endif
