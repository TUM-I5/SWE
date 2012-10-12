// =====================================================================
// This file is part of SWE_CUDA (see file SWE_Block.cu for details).
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
#include "controller.h"

#include "../scenarios/SWE_simple_scenarios_vis.h"
#ifdef ASAGI
#include "../scenarios/SWE_AsagiScenario.hpp"
#endif // ASAGI

/**
    Constructor

	@param sim				instance of simulation class
	@param vis				instance of visualization class

*/	
Controller::Controller(Simulation* sim, Visualization* vis) {
	simulation = sim;
	visualization = vis;
	isActive = true;
	done = false;
	paused = false;
	allowStep = false;

	// No scenario loaded
	memset(scenarios, 0, SCENARIO_COUNT*sizeof(SWE_Scenario*));
}

Controller::~Controller()
{
	// Delete scenarios
	for (int i = 0; i < SCENARIO_COUNT; i++)
		delete scenarios;
}

/**
    Process all user events in a loop
	Returns true, when user wants to quit

*/	
bool Controller::handleEvents() {
	// Handle events
	SDL_Event event;
	allowStep = false;

	while ( SDL_PollEvent(&event) ) 
	{
		switch( event.type )
		{
		case SDL_MOUSEMOTION:
			if (event.motion.state & SDL_BUTTON(SDL_BUTTON_LEFT)) {
				// Left mouse button dragging
				visualization->camera->orient(0.5f*event.motion.yrel, 0.5f*event.motion.xrel);
			}
			else if (event.motion.state & SDL_BUTTON(SDL_BUTTON_RIGHT)) {
				// Panning happens here
				visualization->camera->panning(event.motion.x, event.motion.y);			
			}
			break;
		case SDL_MOUSEBUTTONDOWN:
			// Scroll wheel down
			if ( event.button.button == SDL_BUTTON_WHEELDOWN ) {
				// Zoom out
				visualization->camera->zoomOut(1.2f);
			} else if (event.button.button == SDL_BUTTON_RIGHT) {
				visualization->camera->startPanning(event.button.x, event.button.y);
			break;
		case SDL_MOUSEBUTTONUP:
			// Scroll wheel up
			if ( event.button.button == SDL_BUTTON_WHEELUP ) {
				// Zoom in
				visualization->camera->zoomIn(1.15f);
			}
			} else if (event.button.button == SDL_BUTTON_MIDDLE) {
				visualization->camera->reset();
			}
			break;
		case SDL_ACTIVEEVENT:
			// Don't draw anything if we're getting minimized
			if ( event.active.state & SDL_APPACTIVE ) {
				isActive = (event.active.gain != 0);	
			}
			break;			    
		case SDL_VIDEORESIZE:
			// Our window was resized
			visualization->resizeWindow(event.resize.w, event.resize.h);
			break;
		case SDL_KEYDOWN:
			// User has pressed a key
			done = handleKeyPress( &event.key.keysym);
			break;
		case SDL_QUIT:
			// Window close request
			done = true;
			break;
		default:
			break;
		}
	}
	return done;
}

/**
	Returns true, when window has focus 

*/	
bool Controller::hasFocus() {
	return isActive && !done;
}

/**
    Return whether program is currently paused

*/	
bool Controller::isPaused() {
	return paused && !allowStep;
}

/**
    Process single keyboard event
	@param *keysym			pointer to the sdl keyevent structure

	@todo Refector!!

*/	
bool Controller::handleKeyPress( SDL_keysym *keysym) {

	switch ( keysym->sym ) {
		case SDLK_ESCAPE:
			// ESC key was pressed -> exit
			return true;
			break;
		case SDLK_r:
			// Restart simulation
			allowStep = paused;
			simulation->restart();
			break;
		case SDLK_w:
			// Switch rendering mode
			visualization->toggleRenderingMode();
			break;
		case SDLK_s:
			// Save simulation data to file
			simulation->saveToFile();
			break;
		case SDLK_RIGHT:
			// Advance single timestep when paused
			allowStep = paused;
			break;
		case SDLK_SPACE:
			// Pause/Resume
			paused = !paused;
			break;
		case SDLK_1:
			// Load scenario 1
			{
			  allowStep = paused;

			  if (scenarios[0] == 0)
				  scenarios[0] = new SWE_RadialDamBreakScenario();

			  SWE_Scenario* newScene = scenarios[0];

			  // define grid size and initial time step
			  float dx = (newScene->getBoundaryPos(BND_RIGHT) - newScene->getBoundaryPos(BND_LEFT) )/SWE_Block::getNx();
			  float dy = (newScene->getBoundaryPos(BND_TOP) - newScene->getBoundaryPos(BND_BOTTOM) )/SWE_Block::getNy();
			  SWE_Block::initGridData(SWE_Block::getNx(),SWE_Block::getNy(),dx,dy);

			  simulation->loadNewScenario(newScene, NULL);
			  visualization->updateBathymetryVBO(simulation);
			}
			break;
		case SDLK_2:
			// Load scenario 2
			{
			  allowStep = paused;

			  if (scenarios[1] == 0)
				  scenarios[1] = new SWE_BathymetryDamBreakScenario();

			  SWE_Scenario* newScene = scenarios[1];

			  // define grid size and initial time step
			  float dx = (newScene->getBoundaryPos(BND_RIGHT) - newScene->getBoundaryPos(BND_LEFT) )/SWE_Block::getNx();
			  float dy = (newScene->getBoundaryPos(BND_TOP) - newScene->getBoundaryPos(BND_BOTTOM) )/SWE_Block::getNy();
			  SWE_Block::initGridData(SWE_Block::getNx(),SWE_Block::getNy(),dx,dy);

			  simulation->loadNewScenario(newScene, NULL);
			  visualization->updateBathymetryVBO(simulation);
			}
			break;
		case SDLK_3:
			// Load scenario 3
			{
			  allowStep = paused;

			  if (scenarios[2] == 0)
				  scenarios[2] = new SWE_SplashingPoolScenarioVisInfo();

			  SWE_Scenario* newScene = scenarios[2];

			  // define grid size and initial time step
			  float dx = (newScene->getBoundaryPos(BND_RIGHT) - newScene->getBoundaryPos(BND_LEFT) )/SWE_Block::getNx();
			  float dy = (newScene->getBoundaryPos(BND_TOP) - newScene->getBoundaryPos(BND_BOTTOM) )/SWE_Block::getNy();
			  SWE_Block::initGridData(SWE_Block::getNx(),SWE_Block::getNy(),dx,dy);

			  simulation->loadNewScenario(newScene, NULL);
			  visualization->updateBathymetryVBO(simulation);
			}
			break;
#ifdef ASAGI
		case SDLK_4:
			// Load scenario 3
			{
				allowStep = paused;

				if (scenarios[3] == 0) {
					//simulation area
					float simulationArea[4];
					simulationArea[0] = -450000;
					simulationArea[1] = 6450000;
					simulationArea[2] = -2450000;
					simulationArea[3] = 1450000;
					scenarios[3] = new SWE_AsagiScenario(
							ASAGI_INPUT_DIR "tohoku_gebco_ucsb3_500m_hawaii_bath.nc",
				            ASAGI_INPUT_DIR "tohoku_gebco_ucsb3_500m_hawaii_displ.nc",
				            (float) 28800., simulationArea);
				}

				SWE_Scenario* newScene = scenarios[3];

				  // define grid size and initial time step
				  float dx = (newScene->getBoundaryPos(BND_RIGHT) - newScene->getBoundaryPos(BND_LEFT) )/SWE_Block::getNx();
				  float dy = (newScene->getBoundaryPos(BND_TOP) - newScene->getBoundaryPos(BND_BOTTOM) )/SWE_Block::getNy();
				  SWE_Block::initGridData(SWE_Block::getNx(),SWE_Block::getNy(),dx,dy);

				  simulation->loadNewScenario(newScene, NULL);
				  visualization->updateBathymetryVBO(simulation);
			}
			break;
#endif // ASAGI
		default:
			break;
	}
	return false;
}
