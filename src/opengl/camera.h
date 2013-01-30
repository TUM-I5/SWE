#ifndef CAMERA_H
#define CAMERA_H
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
#include <SDL/SDL.h>
#include <SDL/SDL_opengl.h>

class Camera {
public:
	Camera(const char* window_title);
	// Set modelview matrix
	void setCamera();

	// Change viewing properties
	void reset();

	void viewDistance( float viewDistance );
	void orient( float angX, float angY );
	void zoomIn( float scaleFactor );
	void zoomOut( float scaleFactor );
	void startPanning(int xPos, int yPos);
	void panning(int newX, int newY);

	// Update framebuffer
	void displayImage();

private:
	float view_distance;

	// Position of the camera
	float cameraX;
	float cameraY;
	float cameraZ;

	// Position of the object
	GLdouble objectX;
	GLdouble objectY;
	GLdouble objectZ;

	// Zoom factor
	float zoomfactor;
	float angleX, angleY;

	// Helper variables
	unsigned int frames;
	unsigned int lastTime;
	unsigned int oldMouseX, oldMouseY;
	unsigned int newMouseX, newMouseY;
	
	// Window title
	const char* win_title;
	void rotateObject();
};

#endif
