#ifndef VISUALIZATION_H
#define VISUALIZATION_H
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
#include <SDL.h>
#include <SDL_opengl.h>
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>
#include "camera.h"
#include "simulation.h"
#include "shader.h"
#ifdef USESDLTTF
#include "text.h"
#endif // USESDLTTF

void checkCUDAError(const char *msg);
typedef enum RenderMode {
   SHADED, WIREFRAME, WATERSHADER
} RenderMode;

class Visualization {
public:
	// Constructor and Destructor
	Visualization(int windowWidth, int windowHeight, const char* window_title, int _grid_xsize, int _grid_ysize);
	~Visualization();
	void init(Simulation* sim);
	void cleanUp();
	Camera* camera;

	// Access to CUDA VBO pointers
	cudaGraphicsResource** getCudaNormalsPtr();
	cudaGraphicsResource** getCudaWaterSurfacePtr();

	// Main rendering function
	void renderDisplay();
	
	// Helper functions
	void setRenderingMode(RenderMode mode);
	void updateBathymetryVBO(Simulation* sim);
	void toggleRenderingMode();
	int resizeWindow(int newWidth, int newHeight);
private:	
	// Init helper functions
	void initSDL(int windowWidth, int windowHeight);
	void initGLWindow(int width, int height);
	void initGLDefaults();
	void initCUDA();
	bool IsExtensionSupported(const char* szTargetExtension );

	// Drawing functions
	void DrawWaterSurface(GLuint vboID, GLuint vboNormals, GLuint verticesIndex);
	void DrawBathymetry(GLuint vboID, GLuint verticesIndex);
	void DrawBottom();
	int grid_xsize;
	int grid_ysize;

	// Vertex Buffer objects
	GLuint vboBathymetry;
	GLuint verticesIndex;
	GLuint vboWaterSurface;
	GLuint vboNormals;

	struct cudaGraphicsResource* cuda_vbo_watersurface;
	struct cudaGraphicsResource* cuda_vbo_normals;

	// VBO management functions
	void createIndicesVBO(GLuint* vboID, int xsize, int ysize);
	void createVertexVBO(GLuint* vboID, int size);
	void createBathymetryVBO(GLuint* vboID, int size, Simulation* sim);
	void createVertexVBO(GLuint* vboID, int size, struct cudaGraphicsResource **vbo_res, 
		unsigned int vbo_res_flags);
	void deleteVBO(GLuint* vbo);
	void deleteVBO(GLuint* vbo, struct cudaGraphicsResource *vbo_res);

	// Rendering mode
	RenderMode renderMode;

	// Shader helper class
	Shader* shaders;

#ifdef USESDLTTF
	// Text helper class
	Text* text;
	int windowHeight;
#endif // USESDLTTF

	// Helper function
	int coord(int x, int y, int width = -1);

	// VBO Extension Function Pointers
	PFNGLGENBUFFERSARBPROC glGenBuffers;					// VBO Name Generation Procedure
	PFNGLBINDBUFFERARBPROC glBindBuffer;					// VBO Bind Procedure
	PFNGLBUFFERDATAARBPROC glBufferData;					// VBO Data Loading Procedure
	PFNGLDELETEBUFFERSARBPROC glDeleteBuffers;			// VBO Deletion Procedure
};
#endif
