// =====================================================================
// This file is part of SWE_CUDA (see file SWE_Block.cpp for details).
// 
// Copyright (C) 2010,2011 Tobias Schnabel
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

// Project files
#include "SDL.h"
#include "../opengl/simulation.h"
#include "../opengl/visualization.h"
#include "../opengl/controller.h"
#include "../scenarios/SWE_Scenario.h"
#include "../scenarios/SWE_simple_scenarios.h"
#include "../scenarios/SWE_VtkScenarioVisInfo.h"
#include "../SWE_BlockCUDA.hh"
// #include "../SWE_RusanovBlockCUDA.hh"
#include "../SWE_WavePropagationBlockCuda.hh"

// For SDL compatibility
#undef main

// Display settings
#define SCREEN_WIDTH 800
#define SCREEN_HEIGHT 600
// Number of nodes (not cells) of grid
#define GRID_XSIZE 401
#define GRID_YSIZE 401
#define WINDOW_TITLE "Shallow Water Equations v1.2"

/**
	Main routine.
	Takes 1 command line parameter, which specifies optional input file
	
	Control keys:
	- [ESC]: quit
	- [SPACE]: pause/resume simulation
	- [RIGHT]: advance single step (only when paused)
	- [s]: save current simulation
	- [r]: reset simulation	
    - [w]: toggle rendering mode

	Mouse:
	- Left button: rotating
	- Right button: panning
	- Scroll wheel: zooming

	The class concept relies on the Model-View-Controller (MVC) 
	architecture. 
	The controller class handles all user input and updates the 
	visualization and simulation correspondingly.
	The simulation (model) advances stepwise, the visualization (view) 
	displays then the updated data.
*/
int main(int argc, char *argv[])
{  
	int done = 0;
	int gridx = GRID_XSIZE;
	int gridy = GRID_YSIZE;
	int nx = gridx-1;
	int ny = gridy-1;

	printf("Starting up...\n");

	// Initialize visualization
	Visualization visualization(SCREEN_WIDTH, SCREEN_HEIGHT, WINDOW_TITLE, gridx, gridy);
	printf("Initialized OpenGL window...\n\n");

	// Initialize scenario:
	SWE_Scenario* scene = NULL;
	SWE_VisInfo* visInfo = NULL;
	SWE_BlockCUDA* splash = NULL;
	
	// If input file specified, then read from VTK file:
	if (argc > 1) {
	   SWE_VtkScenarioVisInfo* newScene = SWE_VtkScenarioVisInfo::readVtkFile(argv[1]);
	   // NOTE: Simulation uses a fixed resolution (independent of VTK file)
           scene = newScene;
           visInfo = newScene;
           printf("Scenario read from input file %s\n\n", argv[1]);
	};
	
	if (scene == NULL) { 
	   // ... if VTK file not specified (or was not read successfully)
	   // use splashing pool scenario ...
	   SWE_SplashingPoolScenario* newScene = new SWE_SplashingPoolScenario();
	   scene = newScene;
	};
	
	// define grid size and initial time step
	float dx = (scene->getBoundaryPos(BND_RIGHT) - scene->getBoundaryPos(BND_LEFT) )/nx;
	float dy = (scene->getBoundaryPos(BND_TOP) - scene->getBoundaryPos(BND_BOTTOM) )/ny;

	SWE_Block::initGridData(nx,ny,dx,dy);

    // splash = new SWE_RusanovBlockCUDA();
    splash = new SWE_WavePropagationBlockCuda();
	// define boundary type at all four domain boundaries:
	splash->setWallBoundaries(); // walls at all boundaries

	// Initialize simulation
	printf("Init simulation\n\n");
	Simulation sim(nx,ny,dx,dy, scene, visInfo, splash); 
	printf("Init visualisation\n\n");
	visualization.init(&sim);

	// Initialize controller
	Controller controller(&sim, &visualization);
	
	printf("Start simulation\n\n");
    // sim.saveToFile();
	// Enter the main loop
	while ( ! done ) 
	{
		// Handle events
		done = controller.handleEvents();
		if (controller.isPaused()) {
			// We're paused, only render new display
			visualization.renderDisplay();
		}
		else if (controller.hasFocus()) {
			// Simulate, update visualization data
			sim.runCuda(visualization.getCudaWaterSurfacePtr(), visualization.getCudaNormalsPtr());
			// Render new data
			visualization.renderDisplay();
		}
	}

	// Clean everything up
	visualization.cleanUp();
        
	// delete scene;
	// delete splash;
	
	return 0;
}
