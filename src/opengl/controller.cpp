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

#include "scenarios/SWE_simple_scenarios.hh"
#include "scenarios/SWE_simple_scenarios_vis.hh"
#ifdef ASAGI
#include "scenarios/SWE_AsagiScenario.hh"
#include "scenarios/SWE_AsagiScenario_vis.hh"
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
	memset(visInfos, 0, SCENARIO_COUNT*sizeof(SWE_VisInfo*));
}

Controller::~Controller()
{
	// Delete scenarios
	for (int i = 0; i < SCENARIO_COUNT; i++) {
		delete scenarios[i];
		delete visInfos[i];
	}
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
			}
			break;
		case SDL_MOUSEBUTTONUP:
			// Scroll wheel up
			if ( event.button.button == SDL_BUTTON_WHEELUP ) {
				// Zoom in
				visualization->camera->zoomIn(1.2f);
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
		case SDLK_RIGHT:
			// Advance single timestep when paused
			allowStep = paused;
			break;
		case SDLK_SPACE:
			// Pause/Resume
			paused = !paused;
			break;
		case SDLK_PLUS:
			// Increase water scaling
			visualization->modifyWaterScaling(1.5);
			break;
		case SDLK_MINUS:
			// Decrease water scaling
			visualization->modifyWaterScaling(1/1.5);
			break;
// Not working
// Creating a BlockCUDA with a different resolution seams to be impossible
// at the moment; results in a thrust error.
//		case SDLK_g:
//			// Decrease resolution
//			simulation->resize(0.5);
//			break;
//		case SDLK_h:
//			// Increase resolution
//			simulation->resize(2);
//			break;
		case SDLK_1:
			// Load scenario 1
			{
			  allowStep = paused;

			  if (scenarios[0] == 0)
				  scenarios[0] = new SWE_RadialDamBreakScenario;

			  simulation->loadNewScenario(scenarios[0]);
			  visualization->init(*simulation);
			}
			break;
		case SDLK_2:
			// Load scenario 2
			{
			  allowStep = paused;

			  if (scenarios[1] == 0) {
				  scenarios[1] = new SWE_BathymetryDamBreakScenario();
				  visInfos[1] = new SWE_BathymetryDamBreakVisInfo();
			  }

			  simulation->loadNewScenario(scenarios[1]);
			  visualization->init(*simulation, visInfos[1]);
			}
			break;
		case SDLK_3:
			// Load scenario 3
			{
			  allowStep = paused;

			  if (scenarios[2] == 0)
				  scenarios[2] = new SWE_SplashingPoolScenario();

			  simulation->loadNewScenario(scenarios[2]);
			  visualization->init(*simulation);
			}
			break;
#ifdef ASAGI
		case SDLK_4:
			// Load scenario 4
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

				  simulation->loadNewScenario(scenarios[3]);
				  visualization->init(*simulation);
			}
			break;
		case SDLK_5:
			// Load scenario 5
			{
				allowStep = paused;

				if (scenarios[4] == 0) {
					//simulation area
					float simulationArea[4];
					simulationArea[0] = -450000;
					simulationArea[1] = 700000;
					simulationArea[2] = -1000000;
					simulationArea[3] = 1450000;
					scenarios[4] = new SWE_AsagiScenario(
							ASAGI_INPUT_DIR "tohoku_gebco_ucsb3_500m_hawaii_bath.nc",
				            ASAGI_INPUT_DIR "tohoku_gebco_ucsb3_500m_hawaii_displ.nc",
				            (float) 28800., simulationArea);

					visInfos[4] = new SWE_AsagiJapanSmallVisInfo();
				}

				  simulation->loadNewScenario(scenarios[4]);
				  visualization->init(*simulation, visInfos[4]);
			}
			break;
		case SDLK_6:
			// Load scenario 6
			{
				allowStep = paused;

				if (scenarios[5] == 0) {
					//simulation area
					float simulationArea[4];
					simulationArea[0] = -13775000;
					simulationArea[1] = 1650000;
					simulationArea[2] = -2750000;
					simulationArea[3] = 8840000;
					scenarios[5] = new SWE_AsagiScenario(
							ASAGI_INPUT_DIR "chile_gebco_usgs_500m_bath.nc",
				            ASAGI_INPUT_DIR "chile_gebco_usgs_500m_displ.nc",
				            (float) 28800., simulationArea);
				}

				  simulation->loadNewScenario(scenarios[5]);
				  visualization->init(*simulation);
			}
			break;
		case SDLK_7:
			// Load scenario 7
			{
				allowStep = paused;

				if (scenarios[6] == 0) {
					//simulation area
					float simulationArea[4];
					simulationArea[0] = -2275000;
					simulationArea[1] = 1650000;
					simulationArea[2] = -2265000;
					simulationArea[3] = 1870000;
					scenarios[6] = new SWE_AsagiScenario(
							ASAGI_INPUT_DIR "chile_gebco_usgs_500m_bath.nc",
				            ASAGI_INPUT_DIR "chile_gebco_usgs_500m_displ.nc",
				            (float) 28800., simulationArea);
				}

				  simulation->loadNewScenario(scenarios[6]);
				  visualization->init(*simulation);
			}
			break;
#endif // ASAGI
		default:
			break;
	}
	return false;
}
