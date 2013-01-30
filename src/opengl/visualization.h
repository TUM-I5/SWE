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
#include <SDL/SDL.h>
#include <SDL/SDL_opengl.h>
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>
#include "camera.h"
#include "simulation.h"
#include "shader.h"
#include "vbo.h"
#ifdef USESDLTTF
#include "text.h"
#endif // USESDLTTF

#include "scenarios/SWE_VisInfo.hh"

void checkCUDAError(const char *msg);
typedef enum RenderMode {
   SHADED, WIREFRAME, WATERSHADER
} RenderMode;

class Visualization {
public:
	// Constructor and Destructor
	Visualization(int windowWidth, int windowHeight, const char* window_title);
	~Visualization();

	void init(Simulation &sim, SWE_VisInfo *visInfo = 0L);
	void cleanUp();
	Camera* camera;

	// Access to CUDA VBO pointers
	cudaGraphicsResource** getCudaNormalsPtr();
	cudaGraphicsResource** getCudaWaterSurfacePtr();

	// Main rendering function
	void renderDisplay();
	
	// Rescale water
	void modifyWaterScaling(float factor);

	// Helper functions
	void setRenderingMode(RenderMode mode);
	void toggleRenderingMode();
	int resizeWindow(int newWidth, int newHeight);

	static bool isExtensionSupported(const char* szTargetExtension );
private:
	// Init helper functions
	void initSDL();
	void initGLDefaults();
	void initCUDA();

	void setProjection();

	void updateBathymetryVBO(Simulation &sim);

	// Drawing functions
	void DrawWaterSurface();
	void DrawBathymetry();
	int grid_xsize;
	int grid_ysize;

	// Vertex Buffer objects
	VBO vboBathymetry;
	VBO vboVerticesIndex;
	VBO vboWaterSurface;
	VBO vboNormals;
	// Bathymetry color
	VBO vboBathColor;

	/**
	 * When using trinagle strip mode, this array hold the start
	 * index of each strip
	 */
	GLvoid** indicesOffset;
	/**
	 * When using triangle strip mode, this array holds the length
	 * of each stripe. All our strips have equal length.
	 */
	GLsizei* indicesCount;

	struct cudaGraphicsResource* cuda_vbo_watersurface;
	struct cudaGraphicsResource* cuda_vbo_normals;

	// VBO management functions
	void createIndicesVBO(int xsize, int ysize);
	void createVertexVBO(VBO &vbo, struct cudaGraphicsResource *&vbo_res,
		unsigned int vbo_res_flags);

	void deleteCudaResource(struct cudaGraphicsResource *&vbo_res);

	// Rendering mode
	RenderMode renderMode;

	// Water/Bathymetry scaling/offset
	float wScale, bScale, bOffset;

	/** Location of the water scaling in the shader */
	GLint wScaleLocation;

	// Shaders
	Shader* waterShader;

#ifdef USESDLTTF
	// Text helper class
	Text* text;
#endif // USESDLTTF

	int windowWidth;
	int windowHeight;

	// Helper function
	int coord(int x, int y, int width = -1);

	static void height2Color(float height, GLfloat *color);
	static GLfloat mix(GLfloat a, GLfloat b, float factor)
	{
		return a * (1-factor) + b * factor;
	}

	static PFNGLPRIMITIVERESTARTINDEXNVPROC glPrimitiveRestartIndexNV;
};
#endif
