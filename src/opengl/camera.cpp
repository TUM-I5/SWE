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
#include "camera.h"
#include <sstream>
/**
    Constructor

	@param view_distance	initial view distance from the origin 
	@param window_title		title of the current window

*/
Camera::Camera(const char* window_title)
{
	// Initialize member variables

	win_title = window_title;
	
	reset();

	// Reset framebuffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}
/**
    Set the camera via gluLookAt and set the light position afterwards

*/
void Camera::setCamera() {
	// Clear our buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	// Set camera position, center of scene and up vector 
	float factor = zoomfactor * view_distance;
	gluLookAt(cameraX*factor, cameraY*factor, cameraZ*factor,
		0, 0, 0,  // center of scene
		0, 1, 0); // up vector
	
	// Translate and rotate our object
	glTranslated(objectX, objectY, objectZ);
	rotateObject();
	glScalef(1.0f, 1.0f, -1.0f);
	
	// Position the light source
	// Note: Coordinates are multiplied by current modelview-matrix
	// by OpenGL
	GLfloat LightPosition[] = { 1.0f, 1.0f, 0.0f, 0.0f };
	glLightfv(GL_LIGHT0, GL_POSITION, LightPosition);
}

void Camera::reset()
{
	cameraX = -0.3f;
	cameraY = 1.0f;
	cameraZ = 2.0f;
	objectX = 0.0f;
	objectY = 0.0f;
	objectZ = 0.0f;
	angleX = 0.0f;
	angleY = 0.0f;
	zoomfactor = 0.85f;
	frames = 0;
	lastTime = 0;
	oldMouseX = 0;
	oldMouseY = 0;
}

/**
 * Set the view distance
 */
void Camera::viewDistance( float viewDistance )
{
	view_distance = viewDistance;
}

/**
    Rotate camera 

*/
void Camera::rotateObject() {
	glRotatef(angleX, 1.0, 0.0f, 0.0f);
	glRotatef(angleY, 0.0f, 1.0f, 0.0f);
}

/**
    Increment viewing orientation of the camera

	@param angle_dX		angle relative to the x-axis
	@param angle_dY		angle relative to the rotated y-axis

*/
void Camera::orient( float angle_dX, float angle_dY ) 
{
	angleX += angle_dX;
	angleY += angle_dY;
}

/**
    Zoom in

	@param scaleFactor		factor which is used for zooming			

*/
void Camera::zoomIn( float scaleFactor ) 
{
	if (zoomfactor*view_distance > 3.0f) {
		zoomfactor /= scaleFactor; 
	}	  
}
/**
    Zoom out

	@param scaleFactor		factor which is used for zooming			

*/
void Camera::zoomOut( float scaleFactor ) 
{
	zoomfactor *= scaleFactor;		  
}

/**
    Calculates the current framerate, updates the window title and
	swaps framebuffers to display the new image

*/
void Camera::displayImage() {
    // Set up stringstream for window title
	std::stringstream ss;
	ss.setf(std::ios::fixed, std::ios::floatfield);
	ss.setf(std::ios::showpoint);
	ss.precision(1);
	
	// Calculate framerate
	float fps;
	frames++;
	if (lastTime == 0) {
		lastTime = SDL_GetTicks();
	} else if( (SDL_GetTicks() - lastTime) >= 500) {
		fps = (float)frames/(SDL_GetTicks() - lastTime)*1000.0f;
		ss << win_title << " (" << fps << " fps)";
		SDL_WM_SetCaption(ss.str().c_str(), NULL);

		frames = 0;
        lastTime = SDL_GetTicks();
    }

	// Swap frame buffer
    SDL_GL_SwapBuffers();
}

/**
    User starts dragging. 
	Remember the old mouse coordinates.	

*/
void Camera::startPanning(int xPos, int yPos) {
	oldMouseX = xPos;
	oldMouseY = yPos;
}

/**
    User drags our object. 
	Transform screen coordinates into world coordinates
	and update the objects position

*/
void Camera::panning(int newX, int newY) {
	GLint viewport[4];
	GLdouble modelview[16];
	GLdouble projection[16];
	GLfloat winX, winY, winZ;
	GLdouble posX, posY, posZ;
	GLdouble pos2X, pos2Y, pos2Z;
	GLdouble dX, dY, dZ;
	
	// Draw invisible fake plane
	setCamera();
	glColorMask( GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE );
	GLfloat dim = 100000.0f;
	glBegin(GL_QUADS);	
		glColor3f(1.0f, 1.0f, 0.3f);
		glVertex3f( dim,0.0f, dim);
		glVertex3f(-dim,0.0f, dim);		
		glVertex3f(-dim,0.0f, -dim);	
		glVertex3f( dim,0.0f, -dim);				
	glEnd();
	glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE); 	

    // Get projection matrices
	glGetDoublev( GL_MODELVIEW_MATRIX, modelview );
	glGetDoublev( GL_PROJECTION_MATRIX, projection );
	glGetIntegerv( GL_VIEWPORT, viewport );

	// Get drag positions in object coordinates
	winX = (GLfloat)newX;
	winY = (GLfloat)viewport[3] - (GLfloat)newY;
	glReadPixels( newX, GLint(winY), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ );
	gluUnProject( winX, winY,  winZ, modelview, projection, viewport, &posX, &posY, &posZ);
	
	winX = (GLfloat)oldMouseX;
	winY = (GLfloat)viewport[3] - (GLfloat)oldMouseY;
	glReadPixels( oldMouseX, GLint(winY), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ );
	gluUnProject( winX, winY,  winZ, modelview, projection, viewport, &pos2X, &pos2Y, &pos2Z);
	
	glClear(GL_DEPTH_BUFFER_BIT);
	
	// Calculate drag distances
	dX = (posX-pos2X);
	dY = (posY-pos2Y);
	dZ = (posZ-pos2Z);

	// Convert drag distances to world coordinates
	glPushMatrix();
		glLoadIdentity();
		rotateObject();
		glScalef(1.0f, 1.0f, -1.0f);
		glGetDoublev( GL_MODELVIEW_MATRIX, modelview );
		GLdouble x = modelview[0]*dX + modelview[4]*dY + modelview[8]*dZ; 
		GLdouble y = modelview[1]*dX + modelview[5]*dY + modelview[9]*dZ; 
		GLdouble z = modelview[2]*dX + modelview[6]*dY + modelview[10]*dZ; 
	glPopMatrix();
	
	// Set the new world coordinates of our object
	objectX += x;
	objectY += y;
	objectZ += z;

	// Save our current mouse position
	oldMouseX = newX;
	oldMouseY = newY;
}
