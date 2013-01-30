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
#include "opengl/simulation.h"
#include "opengl/visualization.h"
#include "opengl/controller.h"

#include "tools/Logger.hh"

#include <SDL/SDL.h>

// For SDL compatibility
#undef main

// Display settings
#define SCREEN_WIDTH 800
#define SCREEN_HEIGHT 600
#define WINDOW_TITLE "Shallow Water Equations v1.3"

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

	tools::Logger::logger.printStartMessage();

	// Initialize visualization
	Visualization visualization(SCREEN_WIDTH, SCREEN_HEIGHT, WINDOW_TITLE);
	printf("Initialized OpenGL window...\n\n");

	// Initialize simulation
	printf("Init simulation\n\n");
	Simulation sim;
	printf("Init visualisation\n\n");
	visualization.init(sim);

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
