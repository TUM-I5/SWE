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

/**
	Constructor. All dimensions are node-based, this means a grid consisting of 
	2x2 cells would have 3x3 nodes.

	@param: window_title	title of the window created
	@param: _grid_x_size	number of nodes of the grid (in x-direction)
	@param: _grid_y_size	number of nodes of the grid (in y-direction)

*/
Visualization::Visualization(int windowWidth, int windowHeight, const char* window_title, 
							 int _grid_xsize, int _grid_ysize) {
	// Initialize member variables
	grid_xsize = _grid_xsize;
	grid_ysize = _grid_ysize;
	renderMode = SHADED;

	// Initialize rendering
	initSDL(windowWidth, windowHeight);
	initGLWindow(windowWidth, windowHeight);
	initGLDefaults();
	initCUDA();
	
	// Load camera and shaders
	camera = new Camera(_grid_xsize*1.5f, window_title);
	shaders = new Shader("vertex.glsl", "fragment.glsl");
	if (shaders->shadersSupported()) {
		printf("Shaders supported!\n");
	} else {
		printf("Shaders are NOT supported! Normal rendering mode\n");
	}
	if (shaders->shadersLoaded()) {
		printf("Shaders successfully loaded!\n");
		renderMode = WATERSHADER;
	} else {
		printf("Shaders error while loading shaders\n");
	}
}

/** 
	Destructor (see note below)
*/
Visualization::~Visualization() {
	delete camera;
	delete shaders;
	SDL_Quit();
}

/**
	Frees all memory we used for geometry data 
	Needs to be called before destructor gets called 
	in order to work correctly

*/
void Visualization::cleanUp() {
	deleteVBO(&vboWaterSurface, cuda_vbo_watersurface);
	deleteVBO(&vboNormals, cuda_vbo_normals);
	deleteVBO(&verticesIndex);
	deleteVBO(&vboBathymetry);
}



/**
    Main rendering function. Draws the scene and updates screen

*/	
void Visualization::renderDisplay() {
	// Set camera
	camera->setCamera();
	
	// Draw Scene
	DrawBottom();
	DrawBathymetry(vboBathymetry, verticesIndex);
	
	// Shaded pass
	DrawWaterSurface(vboWaterSurface, vboNormals, verticesIndex);

	// Update framebuffer
	camera->displayImage();  

	// Check for errors
	if (glGetError() != GL_NO_ERROR) {
		printf("OpenGL error occured..\n");
	}
}

/**
    Draws the water surface geometry (triangles)

    @param vboID			id of the vertex buffer object storing 
							the vertex positions
	@param vboID			id of the array buffer object storing
							the corresponding normals for each vertex
	@param verticesIndex	id of the array buffer object storing the 
							vertex indices in rendering sequence 
*/
void Visualization::DrawWaterSurface(GLuint vboID, GLuint vboNormals, GLuint verticesIndex) {
	if (renderMode == WATERSHADER) {
		// Depth pass first
		glColorMask( GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE );
		renderMode = SHADED;
		DrawWaterSurface(vboID, vboNormals, verticesIndex);
		renderMode = WATERSHADER;
		glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);

		// Now render all visible geometry
		shaders->enableShader();
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
	glBindBuffer(GL_ARRAY_BUFFER,vboNormals);
	glNormalPointer(GL_FLOAT, 0, 0);
	glBindBuffer(GL_ARRAY_BUFFER, vboID);
    glVertexPointer(3, GL_FLOAT, 0, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, verticesIndex);
	
	// Enable VBO access and render triangles
	glPushMatrix();
		glTranslatef(-(grid_xsize-1)/2.0f,0.0f,-(grid_ysize-1)/2.0f);
		glColor3f(0.3f*1.2f, 0.45f*1.2f, 0.9f*1.2f);
		if (renderMode == WIREFRAME)
			glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
		glDrawElements(GL_TRIANGLES, 6*(grid_xsize - 1)*(grid_ysize - 1), GL_UNSIGNED_INT, NULL);
		glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
	glPopMatrix();
	
	// Disable array rendering
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisable(GL_LIGHTING);
	glDisable(GL_BLEND);
	shaders->disableShader();

}

/**
    Draws the bathymetry geometry

    @param vboID			id of the vertex buffer object storing 
							the vertex positions
	@param verticesIndex	id of the array buffer object storing the 
							vertex indices in rendering sequence 
*/

void Visualization::DrawBathymetry(GLuint vboID, GLuint verticesIndex) {
	GLsizei stride = 6*sizeof(float);
	
	// Enable lighting
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_LIGHTING);
	
	// Enable array rendering
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_VERTEX_ARRAY);
	
	// Bind buffers
	glBindBuffer(GL_ARRAY_BUFFER, vboID);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, verticesIndex);
	
	// Set VBO pointers
	const GLvoid* normalOffset = (GLvoid*) 12;
	glNormalPointer(GL_FLOAT, stride, normalOffset);
	glVertexPointer(3, GL_FLOAT, stride, 0);

	// Render triangles
	glPushMatrix();
		glTranslatef(-(grid_xsize-1)/2.0f,0.0f,-(grid_ysize-1)/2.0f);
		glColor3f(0.4f*1.1f, 0.36f*1.1f, 0.3f*1.1f);
		glDrawElements(GL_TRIANGLES, 6*(grid_xsize - 1)*(grid_ysize - 1), GL_UNSIGNED_INT, NULL);
	glPopMatrix();

	// Disable array rendering
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisable(GL_LIGHTING);
}

/**
    Draws the bottom plane in our scene

*/
void Visualization::DrawBottom()
{
	glBegin(GL_QUADS);	
		// Draw A Quad
		glNormal3f(0.0f,1.0f,0.0f);
		glColor3f(0.3f,0.3f,0.3f);
		glVertex3f( (grid_xsize - 1)/2.0f,0.0f, (grid_ysize - 1)/2.0f);					// Top Right Of The Quad (Bottom)
		glVertex3f(-(grid_xsize - 1)/2.0f,0.0f, (grid_ysize - 1)/2.0f);					// Top Left Of The Quad (Bottom)
		glVertex3f(-(grid_xsize - 1)/2.0f,0.0f,-(grid_ysize - 1)/2.0f);					// Bottom Left Of The Quad (Bottom)
		glVertex3f( (grid_xsize - 1)/2.0f,0.0f,-(grid_ysize - 1)/2.0f);					// Bottom Right Of The Quad (Bottom)
	glEnd();				
}



/**
	Allocates memory for vertices and other geometry data.
	
	@param sim				instance of the simulation class

*/
void Visualization::init(Simulation* sim) {
	// Create VBOs
	createBathymetryVBO(&vboBathymetry, grid_xsize * grid_ysize * 6 * sizeof(float), sim);
	createIndicesVBO(&verticesIndex, grid_xsize, grid_ysize);
	createVertexVBO(&vboWaterSurface, grid_xsize * grid_ysize * 3 * sizeof(float), 
					&cuda_vbo_watersurface, cudaGraphicsMapFlagsNone);	
	createVertexVBO(&vboNormals ,grid_xsize * grid_ysize * 3 * sizeof(float), 
					&cuda_vbo_normals, cudaGraphicsMapFlagsWriteDiscard);
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
bool Visualization::IsExtensionSupported(const char* szTargetExtension )
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
void Visualization::initSDL(int windowWidth, int windowHeight) {
	// Initialize SDL system
	if ( SDL_Init(SDL_INIT_VIDEO) < 0 ) {
		fprintf(stderr, "Unable to initialize SDL: %s\n", SDL_GetError());
		exit(1);
	}
	
	// Enable double buffering and center window
	SDL_GL_SetAttribute( SDL_GL_DOUBLEBUFFER, 1 );
	SDL_putenv(const_cast<char *>("SDL_VIDEO_CENTERED=center"));
	SDL_WM_SetCaption("Loading...", NULL);

	// Initialize window with OpenGL support
	if ( SDL_SetVideoMode(windowWidth, windowHeight, 0, SDL_OPENGL | SDL_RESIZABLE | SDL_ANYFORMAT) == NULL ) {
		fprintf(stderr, "Unable to create OpenGL screen: %s\n", SDL_GetError());
		SDL_Quit();
		exit(2);
	}

	// Check OpenGL extension(s)
	if (!IsExtensionSupported( "GL_ARB_vertex_buffer_object" )) {
		printf("Vertex Buffer Objects Extension not supported! Exit..\n");
		SDL_Quit();
		exit(1);
	}
	// Load Vertex Buffer Extension
	glGenBuffers = (PFNGLGENBUFFERSARBPROC) SDL_GL_GetProcAddress("glGenBuffersARB");
	glBindBuffer = (PFNGLBINDBUFFERARBPROC) SDL_GL_GetProcAddress("glBindBufferARB");
	glBufferData = (PFNGLBUFFERDATAARBPROC) SDL_GL_GetProcAddress("glBufferDataARB");
	glDeleteBuffers = (PFNGLDELETEBUFFERSARBPROC) SDL_GL_GetProcAddress("glDeleteBuffersARB");
	
}

/**
    Initializes OpenGL projection and viewport settings

    @param	width		width in pixels of the window to create
	@param  height		height in pixels

*/	
void Visualization::initGLWindow(int width, int height) {	
	GLfloat ratio;

	// Protect against a divide by zero
	if ( height == 0 )
		height = 1;

	ratio = ( GLfloat )width / ( GLfloat )height;
	// Setup our viewport.
	glViewport( 0, 0, ( GLint )width, ( GLint )height );

    // change to the projection matrix and set our viewing volume.
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity( );
	float Zoom = 0.00005f;
	glFrustum(-Zoom * width, Zoom * width, -Zoom * height, Zoom * height, 0.1f, 10.0f*grid_ysize);

    // Switch back to the modelview
    glMatrixMode( GL_MODELVIEW );
}

/**
    Gets called when window gets resized

    @param	newWidth		new window width in pixels 
	@param  newHeight		height in pixels

*/
int Visualization::resizeWindow(int newWidth, int newHeight) {
	if ( SDL_SetVideoMode( newWidth, newHeight, 0, 
		SDL_OPENGL | SDL_RESIZABLE | SDL_ANYFORMAT ) == NULL ) 
	{
		fprintf( stderr, "Could not get a surface after resize: %s\n", SDL_GetError( ) );
		return 1;
	} else {
		initGLWindow ( newWidth, newHeight );
		return 0;
	}
	initGLWindow(newWidth, newHeight);
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
    Creates a vertex buffer object in OpenGL 

	@param size				size in bytes to allocate
	@return	vboID			pointer to an integer depicting the id of the 
							new vertex buffer object 

*/	
void Visualization::createVertexVBO(GLuint* vboID, int size)
{
	// Create 1 buffer object
	glGenBuffers(1, vboID);
	glBindBuffer(GL_ARRAY_BUFFER, *vboID);
	
	// Initialize buffer object
	glBufferData(GL_ARRAY_BUFFER, size, NULL, GL_DYNAMIC_DRAW);
	
	// Switch to default buffer
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}

/**
    Creates a vertex buffer object in OpenGL and loads it with the bathymetry
	data from simulation


	@param size				size in bytes to allocate
	@param sim				pointer to an instance of the simulation class
	@return	vboID			pointer to an integer depicting the id of the 
							new vertex buffer object 

*/	
void Visualization::createBathymetryVBO(GLuint* vboID, int size, Simulation* sim) {
	GLfloat* vBathy = new GLfloat[size/sizeof(GLfloat)];

	// Create buffer object for vertex indices
	glGenBuffers(1, vboID);
	glBindBuffer(GL_ARRAY_BUFFER, *vboID);
	
	// Initialize buffer object
	sim->setBathBuffer(vBathy);
	glBufferData(GL_ARRAY_BUFFER, size, vBathy, GL_STATIC_DRAW);
	
	// Switch to default buffer
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	delete[] vBathy;
}

/**
    Updates vertex buffer object with new bathymetry
	data from simulation

	@param sim				pointer to an instance of the simulation class

*/	
void Visualization::updateBathymetryVBO(Simulation* sim) {
	int size = grid_xsize * grid_ysize * 6 * sizeof(float);
	GLfloat* vBathy = new GLfloat[size/sizeof(GLfloat)];
	
	// Select buffer
	glBindBuffer(GL_ARRAY_BUFFER, vboBathymetry);
	
	// Update buffer object
	sim->setBathBuffer(vBathy);
	glBufferData(GL_ARRAY_BUFFER, size, vBathy, GL_STATIC_DRAW);
	
	// Switch to default buffer
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	delete[] vBathy;
}
/**
    Creates a vertex buffer object in OpenGL and an associated CUDA resource

	@param size				size in bytes to allocate 
	@param	vbo_res_flags   cuda flags for memory access
	@return	vboID			pointer to an integer depicting the id of the 
							new vertex buffer object 
	@return	vbo_res			cuda structure created by this function
	
*/	
void Visualization::createVertexVBO(GLuint* vboID, int size, struct cudaGraphicsResource **vbo_res, 
	       unsigned int vbo_res_flags)
{
	createVertexVBO(vboID, size);
	cudaGraphicsGLRegisterBuffer(vbo_res, *vboID, vbo_res_flags);
	checkCUDAError("Couldn't register GL buffer");
}

/**
    Create an array buffer object which holds a list of vertex indices.
	Groups of 3 consecutive vertices corresponding to the indices will be  
	rendered to a triangle.
	This array buffer object therefore describes how single grid points (nodes)
	get transformed into a triangle mesh (tesselation).

    @param	vbo				pointer to an integer depicting the id of a 
							vertex buffer object 
	@param	xsize			number of grid nodes (in x-direction)
	@params ysize			number of grid nodes (in y-direction)
    @return void			
*/	
void Visualization::createIndicesVBO(GLuint* vboID, int xsize, int ysize)
{
	// Create an array describing the vertex indices to be drawn

	int noVertices = (xsize-1)*(ysize-1)*6;
	if ((xsize < 1) || (ysize < 1)) {
		noVertices = 0;
	}
	GLuint* vIndices = new GLuint[noVertices];

	for(int y=0; y < ysize - 1; y++)
	{
		for(int x=0; x < xsize - 1; x++)
		{
			// Create tessellation of the grid
			vIndices[coord(x,y, xsize - 1)*6] = coord(x,y);
			vIndices[coord(x,y, xsize - 1)*6 + 1] = coord(x+1,y+1);
			vIndices[coord(x,y, xsize - 1)*6 + 2] = coord(x,y+1);
			vIndices[coord(x,y, xsize - 1)*6 + 3] = coord(x,y);
			vIndices[coord(x,y, xsize - 1)*6 + 4] = coord(x+1,y);
			vIndices[coord(x,y, xsize - 1)*6 + 5] = coord(x+1,y+1);
			
		}
	}

	// Create buffer object for vertex indices
	glGenBuffers(1, vboID);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, *vboID);
	
	// Initialize buffer object
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*noVertices, vIndices, GL_STATIC_DRAW);
	
	// Switch to default buffer
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	delete[] vIndices;
}

/**
    Frees memory used by a vertex buffer object

    @param	vbo				pointer to an integer depicting the id of a 
							vertex buffer object 
*/	
void Visualization::deleteVBO(GLuint* vbo)
{
    if (vbo) {
		glBindBuffer(1, *vbo);
		glDeleteBuffers(1, vbo);
		*vbo = 0;
    }
}
/**
    Frees memory used by a vertex buffer object and a CUDA resource

    @param	vbo				pointer to an integer depicting the id of a 
							vertex buffer object 
	@param	vbo_res			pointer to a CUDA resource structure		
*/	
void Visualization::deleteVBO(GLuint* vbo, struct cudaGraphicsResource *vbo_res)
{
    if (vbo) {
		cudaGraphicsUnregisterResource(vbo_res);
		deleteVBO(vbo);
    }
}

/**
    Initialize CUDA device
*/	
void Visualization::initCUDA() {
	int driver, runtime;
	cudaDeviceProp  prop;
    int dev;

	// Display CUDA Versions
	cudaDriverGetVersion(&driver);
	cudaRuntimeGetVersion(&runtime);
	printf("CUDA Driver Version: %d, Runtime Version: %d\n", driver, runtime);
	
	// Select GPU for CUDA and OpenGL
    memset( &prop, 0, sizeof( cudaDeviceProp ) );
    prop.major = 1;
    prop.minor = 0;
    cudaChooseDevice( &dev, &prop );
    cudaGLSetGLDevice( dev ) ;
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
			if (shaders->shadersLoaded()) {
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

