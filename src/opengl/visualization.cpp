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

#include "visualization.h"
#include "tools/Logger.hh"

#include <limits>

/**
	Constructor. All dimensions are node-based, this means a grid consisting of 
	2x2 cells would have 3x3 nodes.

	@param: window_title	title of the window created
	@param: _grid_x_size	number of nodes of the grid (in x-direction)
	@param: _grid_y_size	number of nodes of the grid (in y-direction)

*/
Visualization::Visualization(int windowWidth, int windowHeight, const char* window_title)
	: indicesOffset(0L),
	  indicesCount(0L),
	  camera(0L),
	  windowWidth(windowWidth),
	  windowHeight(windowHeight)
{
	// Initialize member variables
	renderMode = SHADED;

	cuda_vbo_watersurface = 0L;
	cuda_vbo_normals = 0L;

	// Initialize rendering
	initSDL();
	initGLDefaults();
	initCUDA();
#ifdef USESDLTTF
	text = new Text();
	text->addText("Keys:");
#ifdef ASAGI
	text->addText("  1-7: Select scenario");
#else // ASAGI
	text->addText("  1-3: Select scenario");
#endif // ASAGI
	text->addText("  Space: Pause/Resume");
	text->addText("  ->: Next frame (when paused)");
	text->addText("  r: Restart scenario");
	text->addText("  +/-: Scale wave height");
	text->addText("Mouse:");
	text->addText("  Left button: rotate");
	text->addText("  Right button: move");
	text->addText("  Middle button: reset camera");
	text->addText("  Wheel: zoom in/out");
#endif // USESDLTTF
	
	// Load camera and shaders
	camera = new Camera(window_title);
	waterShader = new Shader("vertex.glsl", "fragment.glsl");
	if (waterShader->shadersLoaded()) {
		tools::Logger::logger.printString("Water shaders successfully loaded!");
		renderMode = WATERSHADER;
	} else
		tools::Logger::logger.printString("Error while loading water shaders");

	wScaleLocation = waterShader->getUniformLocation("scale");

	// Create openGL buffers
	vboBathymetry.init();
	vboVerticesIndex.init();
	vboWaterSurface.init();
	vboNormals.init();
	vboBathColor.init();
}

/** 
	Destructor (see note below)
*/
Visualization::~Visualization() {
	delete camera;
	delete waterShader;

#ifdef USESDLTTF
	delete text;
#endif // USESDLTTF

	delete [] indicesOffset;
	delete [] indicesCount;

	SDL_Quit();
}

/**
	Allocates memory for vertices and other geometry data.

	@param sim				instance of the simulation class
*/
void Visualization::init(Simulation &sim, SWE_VisInfo *visInfo)
{
	// Update the grid size
	grid_xsize = sim.getNx()+1;
	grid_ysize = sim.getNy()+1;

	// Set the camera distance depending on the grid size
	camera->viewDistance(1.5*grid_xsize);

	// Update the GL projection matrix
	setProjection();

	// Create indices buffer
	createIndicesVBO(grid_xsize, grid_ysize);

	// Set bathymetry
	updateBathymetryVBO(sim);

	// Create buffers for water
	createVertexVBO(vboWaterSurface, cuda_vbo_watersurface, cudaGraphicsMapFlagsNone);
	createVertexVBO(vboNormals,	cuda_vbo_normals, cudaGraphicsMapFlagsWriteDiscard);

	if (visInfo == 0L) {
		sim.getScalingApproximation(bScale, bOffset, wScale);
	} else {
		bScale = visInfo->bathyVerticalScaling();
		bOffset = visInfo->bathyVerticalOffset();
		wScale = visInfo->waterVerticalScaling();
	}
}

/**
	Frees all memory we used for geometry data 
	Needs to be called before destructor gets called 
	in order to work correctly
*/
void Visualization::cleanUp() {
	deleteCudaResource(cuda_vbo_watersurface);
	deleteCudaResource(cuda_vbo_normals);

	vboBathymetry.finialize();
	vboVerticesIndex.finialize();
	vboWaterSurface.finialize();
	vboNormals.finialize();
	vboBathColor.finialize();
}


/**
    Main rendering function. Draws the scene and updates screen

*/	
void Visualization::renderDisplay() {
	// Set camera
	camera->setCamera();
	
	// Draw Scene
	//DrawBottom();
	DrawBathymetry();
	
	// Shaded pass
	DrawWaterSurface();

#ifdef USESDLTTF
	text->startTextMode();
	SDL_Rect location = {5, windowHeight-5, 0, 0};
	while (text->showNextText(location))
		location.y -= location.h;
	text->endTextMode();
#endif // USESDLTTF

	// Update framebuffer
	camera->displayImage();  

	int errCode = glGetError();
	// Check for errors
	if (errCode != GL_NO_ERROR) {
		printf("OpenGL error occured: %s\n", gluErrorString(errCode));
	}
}

/**
    Draws the water surface geometry (triangles)
 */
void Visualization::DrawWaterSurface()
{
	if (renderMode == WATERSHADER) {
		// Depth pass first
		glColorMask( GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE );
		renderMode = SHADED;
		DrawWaterSurface();
		renderMode = WATERSHADER;
		glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);

		// Now render all visible geometry
		waterShader->enableShader();
		// Set shader parameter
		waterShader->setUniform(wScaleLocation, wScale);

		glEnable(GL_BLEND);
		glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}
	else if (renderMode == SHADED) {
		// Enable lighting
		glEnable(GL_COLOR_MATERIAL);
		glEnable(GL_LIGHTING);
	} else {
		glDisable(GL_LIGHTING);
	}
	
	// Enable array rendering
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_VERTEX_ARRAY);

	// Set rendering to VBO mode
	vboNormals.bindBuffer();
	glNormalPointer(GL_FLOAT, 0, 0);
	vboWaterSurface.bindBuffer();
    glVertexPointer(3, GL_FLOAT, 0, 0);
	vboVerticesIndex.bindBuffer(GL_ELEMENT_ARRAY_BUFFER);

	// Enable VBO access and render triangles
	glPushMatrix();
		glScalef(1.0f, wScale, 1.0f);
		glTranslatef(-(grid_xsize-1)/2.0f, 0.0f, -(grid_ysize-1)/2.0f);
		glColor3f(0.3f*1.2f, 0.45f*1.2f, 0.9f*1.2f);
		if (renderMode == WIREFRAME)
			glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );

		if (glPrimitiveRestartIndexNV) {
			glDrawElements(GL_TRIANGLE_STRIP, 2*grid_xsize*(grid_ysize-1)+grid_ysize-1,
					GL_UNSIGNED_INT, 0L);
		} else
			glDrawElements(GL_TRIANGLES, 6*(grid_xsize-1)*(grid_ysize-1), GL_UNSIGNED_INT, NULL);

		glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
	glPopMatrix();
	
	// Disable array rendering
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisable(GL_LIGHTING);
	glDisable(GL_BLEND);
	waterShader->disableShader();

}

/**
    Draws the bathymetry geometry
 */

void Visualization::DrawBathymetry() {
	GLsizei stride = 6*sizeof(float);

	// Enable lighting
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_LIGHTING);
	
	// Enable array rendering
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);
	
	// Bind buffers
	vboBathymetry.bindBuffer();
	vboVerticesIndex.bindBuffer(GL_ELEMENT_ARRAY_BUFFER);

	// Set VBO pointers
	const GLvoid* normalOffset = (GLvoid*) 12;
	glNormalPointer(GL_FLOAT, stride, normalOffset);
	glVertexPointer(3, GL_FLOAT, stride, 0);
	vboBathColor.bindBuffer();
	glColorPointer(3, GL_FLOAT, 0, 0);

	// Render triangles
	glPushMatrix();
		glScalef(1.0f, bScale, 1.0f);
		glTranslatef(-(grid_xsize-1)/2.0f, bOffset, -(grid_ysize-1)/2.0f);

		if (glPrimitiveRestartIndexNV) {
			glDrawElements(GL_TRIANGLE_STRIP, 2*grid_xsize*(grid_ysize-1)+grid_ysize-1,
					GL_UNSIGNED_INT, 0L);
		} else
			glDrawElements(GL_TRIANGLES, 6*(grid_xsize-1)*(grid_ysize-1), GL_UNSIGNED_INT, NULL);
	glPopMatrix();

	// Disable array rendering
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);
	glDisable(GL_LIGHTING);
}

/**
	Returns a pointer to the cuda memory object holding the
	vertex normals
*/
cudaGraphicsResource** Visualization::getCudaNormalsPtr() {
	return &cuda_vbo_normals;
}
/**
	Returns a pointer to the cuda memory object holding the
	vertex positions
*/
cudaGraphicsResource** Visualization::getCudaWaterSurfacePtr() {
	return &cuda_vbo_watersurface;
}

/**
    Returns, whether a special extension is supported by the current 
	graphics card

    @param	szTargetExtention	string describing the extension to look for	

*/	
bool Visualization::isExtensionSupported(const char* szTargetExtension )
{
	const unsigned char *pszExtensions = NULL;
	const unsigned char *pszStart;
	unsigned char *pszWhere, *pszTerminator;

	// Extension names should not have spaces
	pszWhere = (unsigned char *) strchr( szTargetExtension, ' ' );
	if( pszWhere || *szTargetExtension == '\0' )
		return false;

	// Get Extensions String
	pszExtensions = glGetString( GL_EXTENSIONS );

	// Search The Extensions String For An Exact Copy
	pszStart = pszExtensions;
	for(;;)
	{
		pszWhere = (unsigned char *) strstr( (const char *) pszStart, szTargetExtension );
		if( !pszWhere )
			break;
		pszTerminator = pszWhere + strlen( szTargetExtension );
		if( pszWhere == pszStart || *( pszWhere - 1 ) == ' ' )
			if( *pszTerminator == ' ' || *pszTerminator == '\0' )
				return true;
		pszStart = pszTerminator;
	}
	return false;
}

/**
    Initialize the simple directmedia layer (SDL) which is a wrapper library 
	for OpenGL

    @param	windowWidth		width in pixels of the window to create
	@param  windowHeight	height in pixels

*/	
void Visualization::initSDL() {
	// Initialize SDL system
	if ( SDL_Init(SDL_INIT_VIDEO) < 0 ) {
		fprintf(stderr, "Unable to initialize SDL: %s\n", SDL_GetError());
		exit(1);
	}
	
	// Enable double buffering and center window
	SDL_GL_SetAttribute( SDL_GL_DOUBLEBUFFER, 1 );
	SDL_putenv(const_cast<char *>("SDL_VIDEO_CENTERED=center"));
	SDL_WM_SetCaption("Loading...", NULL);

	// Disable swap control
	SDL_GL_SetAttribute(SDL_GL_SWAP_CONTROL, 0);

	// Initialize window with OpenGL support
	if ( SDL_SetVideoMode( windowWidth, windowHeight, 0,
		SDL_OPENGL | SDL_RESIZABLE | SDL_ANYFORMAT ) == NULL )
	{
		fprintf(stderr, "Unable to create OpenGL screen: %s\n", SDL_GetError());
		SDL_Quit();
		exit(2);
	}
}

void Visualization::setProjection()
{
	GLfloat ratio;

	// Protect against a divide by zero
	if ( windowHeight == 0 )
		windowHeight = 1;

	ratio = ( GLfloat )windowWidth / ( GLfloat )windowHeight;

	// Setup our viewport.
	glViewport( 0, 0, ( GLint )windowWidth, ( GLint )windowHeight );

    // change to the projection matrix and set our viewing volume.
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity( );
	float Zoom = 0.00005f;
	glFrustum(-Zoom * windowWidth, Zoom * windowWidth, -Zoom * windowHeight, Zoom * windowHeight, 0.1f, 20.0f*grid_ysize);

    // Switch back to the modelview
    glMatrixMode( GL_MODELVIEW );
}

/**
    Gets called when window gets resized

    @param	newWidth		new window width in pixels 
	@param  newHeight		height in pixels

*/
int Visualization::resizeWindow(int newWidth, int newHeight)
{
	windowWidth = newWidth;
	windowHeight = newHeight;

	if ( SDL_SetVideoMode( newWidth, newHeight, 0,
		SDL_OPENGL | SDL_RESIZABLE | SDL_ANYFORMAT ) == NULL )
	{
		fprintf( stderr, "Could not get a surface after resize: %s\n", SDL_GetError( ) );
		return 1;
	}

	setProjection();

	return 0;
}

/**
    Initializes various OpenGL settings
*/
void Visualization::initGLDefaults() {
	// Set background color to white
	glClearColor(1.0f, 1.0f, 1.0f, 0.0f);

    // Smooth (interpolated) shading
    glShadeModel( GL_SMOOTH );

    // Reset the depth buffer
    glClearDepth( 1.0f );

    // Enable depth testing
    glEnable( GL_DEPTH_TEST );

    // Comparisons in depth buffer via <=
    glDepthFunc( GL_LEQUAL );

    // Some perspective corrections
    glHint( GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST );

	// Ensure that all of the following transformations will only
	// apply to our modelview matrix (objects + camera)
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// Lighting
	GLfloat LightAmbient[]=		{ 0.0f, 0.0f, 0.0f, 1.0f };
	GLfloat LightDiffuse[]=		{ 1.0f, 1.0f, 1.0f, 1.0f };
	GLfloat LightPosition[]=	{ 1.0f, 1.0f, 0.0f, 0.0f };
	glLightfv(GL_LIGHT0, GL_POSITION,LightPosition);	// Position The Light
	glLightfv(GL_LIGHT0, GL_AMBIENT, LightAmbient);		// Setup The Ambient Light
	glLightfv(GL_LIGHT0, GL_DIFFUSE, LightDiffuse);		// Setup The Diffuse Light
	glEnable(GL_LIGHT0);	
	glEnable(GL_NORMALIZE);
	glDisable(GL_LIGHTING);

	glPrimitiveRestartIndexNV = (PFNGLPRIMITIVERESTARTINDEXNVPROC) SDL_GL_GetProcAddress("glPrimitiveRestartIndexNV");
	if (glPrimitiveRestartIndexNV) {
		glPrimitiveRestartIndexNV(std::numeric_limits<GLuint>::max());
		glEnableClientState(GL_PRIMITIVE_RESTART_NV);
	} else {
		tools::Logger::logger.printString("glPrimitiveRestartIndexNV not found");
		tools::Logger::logger.printString("Can not use tringular strips");
	}
	// TODO: Add support for glPrimitiveRestartIndex (OpenGL 3.1)
}

/**
    Calculate 1-D array offset for 2-D coordinates

    @param	x			x-position
	@param	y			y-position
	@param	width		size of a line

*/	
int Visualization::coord(int x, int y, int width) {
	if (width == -1) width = grid_xsize;
	return (y*width + x);
}

/**
    Updates vertex buffer object with new bathymetry
	data from simulation

	@param sim				pointer to an instance of the simulation class

*/	
void Visualization::updateBathymetryVBO(Simulation &sim) {
	int size = grid_xsize * grid_ysize * 6;
	
	GLfloat* vBathy = new GLfloat[size];
	sim.setBathBuffer(vBathy);
	
	vboBathymetry.setBufferData(size * sizeof(GLfloat), vBathy);
	
	delete[] vBathy;

	const Float2D &bathymetry = sim.getBathymetry();
	GLfloat* color = new GLfloat[grid_xsize*grid_ysize*3];
	for (int i = 0; i < grid_xsize; i++) {
		for (int j = 0; j < grid_ysize; j++) {
			height2Color(bathymetry[i][j], &color[((j*grid_ysize)+i)*3]);
		}
	}

	vboBathColor.setBufferData(grid_xsize*grid_ysize*3*sizeof(GLfloat), color);
	delete[] color;
}
/**
    Creates a vertex buffer object in OpenGL and an associated CUDA resource

	@param vbo				Vertex buffer object
	@param size				size in bytes to allocate 
	@param vbo_res_flags   cuda flags for memory access
	@return	vbo_res			cuda structure created by this function
	
*/	
void Visualization::createVertexVBO(VBO &vbo, struct cudaGraphicsResource *&vbo_res,
	       unsigned int vbo_res_flags)
{
	deleteCudaResource(vbo_res);

	vbo.setBufferData(grid_xsize * grid_ysize * 3 * sizeof(float), 0L,
			GL_ARRAY_BUFFER, GL_DYNAMIC_DRAW);

	cudaGraphicsGLRegisterBuffer(&vbo_res, vbo.getName(), vbo_res_flags);
	checkCUDAError("Couldn't register GL buffer");
}

/**
    Create an array buffer object which holds a list of vertex indices.
	Groups of 3 consecutive vertices corresponding to the indices will be  
	rendered to a triangle.
	This array buffer object therefore describes how single grid points (nodes)
	get transformed into a triangle mesh (tesselation).

	@param	xsize			number of grid nodes (in x-direction)
	@params ysize			number of grid nodes (in y-direction)
*/	
void Visualization::createIndicesVBO(int xsize, int ysize)
{
	// Create an array describing the vertex indices to be drawn
	int noVertices;
	if (glPrimitiveRestartIndexNV)
		noVertices = (xsize)*(ysize-1)*2 + ysize - 1;
	else
		noVertices = (xsize-1)*(ysize-1)*6;

	GLuint* vIndices = new GLuint[noVertices];

	// Create tessellation of the grid
	if (glPrimitiveRestartIndexNV) {
		for (int y = 0; y < ysize - 1; y++) {
			for (int x = 0; x < xsize; x++) {
				vIndices[coord(x, y, xsize)*2 + y] = coord(x, y);
				vIndices[coord(x, y, xsize)*2 + 1 + y] = coord(x, y+1);
			}

			// Add restart flag
			vIndices[coord(0, y+1, xsize)*2 + y] = std::numeric_limits<GLuint>::max();
		}

	} else {
		for(int y = 0; y < ysize - 1; y++) {
			for(int x = 0; x < xsize - 1; x++) {
				vIndices[coord(x,y, xsize - 1)*6] = coord(x,y);
				vIndices[coord(x,y, xsize - 1)*6 + 1] = coord(x+1,y+1);
				vIndices[coord(x,y, xsize - 1)*6 + 2] = coord(x,y+1);
				vIndices[coord(x,y, xsize - 1)*6 + 3] = coord(x,y);
				vIndices[coord(x,y, xsize - 1)*6 + 4] = coord(x+1,y);
				vIndices[coord(x,y, xsize - 1)*6 + 5] = coord(x+1,y+1);
			}
		}
	}
	
	// Initialize buffer object
	vboVerticesIndex.setBufferData(noVertices * sizeof(GLuint), vIndices,
			GL_ELEMENT_ARRAY_BUFFER, GL_STATIC_DRAW);
	
	delete[] vIndices;
}

/**
    Frees memory used by a vertex buffer object and a CUDA resource

	@param	vbo_res			pointer to a CUDA resource structure		
*/	
void Visualization::deleteCudaResource(struct cudaGraphicsResource *&vbo_res)
{
    if (vbo_res != 0L) {
		cudaGraphicsUnregisterResource(vbo_res);
		vbo_res = 0L;
    }
}

/**
    Initialize CUDA device
*/	
void Visualization::initCUDA() {
	int driver, runtime;
	cudaDeviceProp  prop;
    int dev;

	// Init cuda device
	SWE_BlockCUDA::init();

	// Display CUDA Versions
	cudaDriverGetVersion(&driver);
	cudaRuntimeGetVersion(&runtime);
	printf("CUDA Driver Version: %d, Runtime Version: %d\n", driver, runtime);
	
	// Select GPU for CUDA and OpenGL
    memset(&prop, 0, sizeof(cudaDeviceProp));
    prop.major = 1;
    prop.minor = 0;
    cudaChooseDevice( &dev, &prop );
    cudaGLSetGLDevice( dev ) ;
}

void Visualization::modifyWaterScaling(float factor)
{
	wScale *= factor;
}

/**
    Sets current rendering mode

    @param mode			rendering mode			
*/	
void Visualization::setRenderingMode(RenderMode mode) {
	renderMode = mode;
}

/**
    Switches between 3 different rendering modes:
	- Shaded: Use OpenGL shading
	- Wireframe: Only render edges of each triangle
	- Watershader: Use custom GLSL shader for water surface
*/		
void Visualization::toggleRenderingMode() {
	switch( renderMode ) {
		case SHADED:
			renderMode = WIREFRAME;
			break;
		case WIREFRAME:
			// Skip watershader if shaders not loaded
			if (waterShader->shadersLoaded()) {
				renderMode = WATERSHADER;
			} else {
				renderMode = SHADED;
			}
			break;
		case WATERSHADER:
			renderMode = SHADED;
			break;
	}
}

void Visualization::height2Color(float height, GLfloat *color)
{
	// Workaround "wrong" offset in colormap
	height += 150;

	if (height < -9000.0) {
		color[0] = 0.0;
		color[1] = 0.0;
		color[2] = 0.2;
	} else if (height < -8525.07) {
		color[0] = 0.0;
		color[1] = 0.0;
		color[2] = mix(0.2, 1.0, (-9000-height)/(-9000+8525.07));
	} else if (height < 189) {
		color[0] = 0;
		color[1] = mix(0.0, 1.0, (-8525.07-height)/(-8525.07-189));
		color[2] = 1;
	} else if (height < 190) {
		float factor = (189-height)/(189-190);
		color[0] = 0.0;
		color[1] = mix(1.0, 0.4, factor);
		color[2] = mix(1.0, 0.2, factor);
	} else if (height < 1527.7) {
		float factor = (190-height)/(190-1527.7);
		color[0] = mix(0.0, 0.952941, factor);
		color[1] = mix(0.4, 0.847059, factor);
		color[2] = mix(0.2, 0.415686, factor);
	} else if (height < 4219) {
		float factor = (1527.7-height)/(1527.7-4219);
		color[0] = mix(0.952941, 0.419577, factor);
		color[1] = mix(0.847059, 0.184253, factor);
		color[2] = mix(0.415686, 0.00648508, factor);
	} else if (height < 4496.04) {
		float factor = (4219-height)/(4219-4496.04);
		color[0] = mix(0.419577, 0.983413, factor);
		color[1] = mix(0.184253, 0.9561, factor);
		color[2] = mix(0.00648508, 0.955749, factor);
	} else if (height < 6000) {
		float factor = (4496.04-height)/(4496.04-6000);
		color[0] = mix(0.983413, 1.0, factor);
		color[1] = mix(0.9561, 1.0, factor);
		color[2] = mix(0.955749, 1.0, factor);
	} else {
		color[0] = 1.0;
		color[1] = 1.0;
		color[2] = 1.0;
	}
}

PFNGLPRIMITIVERESTARTINDEXNVPROC Visualization::glPrimitiveRestartIndexNV;
