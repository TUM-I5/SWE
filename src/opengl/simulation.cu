// =====================================================================
// This file is part of SWE_CUDA (see file SWE_Block.cu for details).
// 
// Copyright (C) 2010,2011 Michael Bader, Kaveh Rahnema, Tobias Schnabel
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
#include "simulation.h"

#include "blocks/cuda/SWE_WavePropagationBlockCuda.hh"

#include <cstring>

// Taken form FWaveCuda.h
// TODO: Put it in a common header file
const float dryTol = 100.;

#define DEFAULT_NX 560
#define DEFAULT_NY 560

/**
    Constructor. 
	Initializes SWE_BlockCUDA and creates a new instance of it.

*/
Simulation::Simulation ()
	: block(0L), scenario(0L),
	  nx(DEFAULT_NX), ny(DEFAULT_NY)
{
	loadNewScenario(&defaultScenario);
}

/** 
	Destructor.
*/
Simulation::~Simulation () {
	delete block;
}

void Simulation::loadNewScenario(SWE_Scenario* scene)
{
	// Load new scene
	scenario = scene;

	restart();
}

void Simulation::resize(float factor)
{
	this->nx *= factor;
	this->ny *= factor;

	restart();
}

/**
    This is the main simulation procedure. Simulates one timestep and computes
	normals afterwards.
	
    @param vbo_resource			cuda resource holding the vertex positions
	@param vbo_normals			cuda resource holding the normals
*/
void Simulation::runCuda(struct cudaGraphicsResource **vbo_resource, struct cudaGraphicsResource **vbo_normals)
{
    // map OpenGL buffer object for writing from CUDA
    float3 *dptr, *dptr2;
	std::size_t num_bytes, num_bytes2;
    cudaGraphicsMapResources(1, vbo_resource, 0);
    
    cudaGraphicsResourceGetMappedPointer((void **)&dptr, &num_bytes,  
						       *vbo_resource);

	// run simulation and fill our vbo
	calculateWaterSurface(dptr);

    // unmap buffer object
    cudaGraphicsUnmapResources(1, vbo_resource, 0);
	checkCUDAError("Fehler bei Simulations-VBO\n");

	//Calculate normals
    cudaGraphicsMapResources(1, vbo_resource, 0);
    cudaGraphicsResourceGetMappedPointer((void **)&dptr, &num_bytes,  
						       *vbo_resource);
	cudaGraphicsMapResources(1, vbo_normals, 0);
    cudaGraphicsResourceGetMappedPointer((void **)&dptr2, &num_bytes2,  
						       *vbo_normals);

	// calculate normals of our new water surface
	calculateNormals(dptr, dptr2);

    // unmap buffer objects
	cudaGraphicsUnmapResources(1, vbo_resource, 0);
	checkCUDAError("Fehler bei Simulations-VBO\n");
    cudaGraphicsUnmapResources(1, vbo_normals, 0);
	checkCUDAError("Fehler bei Normalen-VBO\n");
}


/**
    Restarts the simulation. Restores the initial bondaries.

*/
void Simulation::restart()
{
	curTime = 0.0f;
	isFirstStep = 1;

	// define grid size
	float dx = (scenario->getBoundaryPos(BND_RIGHT) - scenario->getBoundaryPos(BND_LEFT) )/nx;
	float dy = (scenario->getBoundaryPos(BND_TOP) - scenario->getBoundaryPos(BND_BOTTOM) )/ny;

	// Create the wavepropagation block
	delete block;
	block = new SWE_WavePropagationBlockCuda(nx, ny, dx, dy);

	// Initialize the scenario
	float l_originX, l_originY;

	// get the origin from the scenario
	l_originX = scenario->getBoundaryPos(BND_LEFT);
	l_originY = scenario->getBoundaryPos(BND_BOTTOM);

	block->initScenario(l_originX, l_originY, *scenario);

	block->setGhostLayer();
}


/**
    Sets the bathymetry buffer. Buffer contains vertex position and 
    vertex normal in sequence.
    @param bath				float array in which computed values will be stored
*/
void Simulation::setBathBuffer(float* bath) {
        const Float2D& b = block->getBathymetry();
        int nx = b.getRows()-2;
        int ny = b.getCols()-2;
        
        // Set vertex coordinates
	for (int j=0; j<ny+1;j++) {
		for (int i=0;i<nx+1;i++) {
				bath[(j*(ny+1) + i)*6] = (float) i;
				bath[(j*(ny+1) + i)*6 + 1]= 0.25f * (b[i][j]+b[i+1][j]+b[i][j+1]+b[i+1][j+1]);
				bath[(j*(ny+1) + i)*6 + 2] = (float) j;
				bath[(j*(ny+1) + i)*6 + 3] = 0.0f;
				bath[(j*(ny+1) + i)*6 + 4] = 0.0f;
				bath[(j*(ny+1) + i)*6 + 5] = 0.0f;
		}
	}
	// Calculate normals
	for(int j=0; j < ny; j++)
	{
		for(int i=0; i < nx; i++)
		{
			// Calculate normal vectors for each triangle
			float normal1[3];
			float normal2[3];

			calculateNormal(&bath[(j*(ny+1) + i)*6], 
				&bath[((j+1)*(ny+1) + i + 1)*6],
				&bath[((j+1)*(ny+1) + i)*6],
				normal1);

			calculateNormal(&bath[(j*(ny+1) + i)*6],
				&bath[(j*(ny+1) + i + 1)*6],
				&bath[((j+1)*(ny+1) + i + 1)*6],
				normal2);
			// Copy normals to array
			for (int k=0; k < 3; k++) {
				bath[(j*(ny+1) + i)*6 + 3 + k] = (normal1[k]+normal2[k])*0.5f;
			}
		}
	}

	// Fill boundary regions
	for(int x=0; x < ny; x++) {
		for (int i=0; i < 3; i++) {
			bath[(ny*(ny+1) + x)*6 + 3 + i] = bath[((ny-1)*(ny+1) + x)*6 + 3 + i];
		}
	}
	for(int y=0; y < nx; y++) {
		for (int i=0; i < 3; i++) {
			bath[(y*(ny+1) + nx)*6 + 3 + i] = bath[(y*(ny+1) + nx - 1)*6 + 3 + i];
		}
	}
	for (int i=0; i < 3; i++) {
			bath[(ny*(ny+1) + nx)*6 + 3 + i] = bath[((ny-1)*(ny+1) + nx - 1)*6 + 3 + i];
	}
}

/**
    Advance single step in simulation on graphics card 
	and copy results to visualization buffer

    @param destBuffer			vertex buffer object holding all water surface data
*/
void Simulation::calculateWaterSurface(float3* destBuffer) {
	if (isFirstStep == 1) {
		isFirstStep = 0;
	} else {
		block->setGhostLayer();
		block->computeNumericalFluxes();
		float dt = block->getMaxTimestep();
		block->updateUnknowns(dt);
		curTime += dt;
	}
	// splash->updateVisBuffer(destBuffer, wAverage, wScale, wOffset);
	updateVisBuffer(destBuffer);
}

/**
    Computes a first approximation of the scaling values needed
	for visualization.
	Gets called before simulation starts and determines the average,
	mininimum and maximum values of the bathymetry and water surface data.
	Uses latter values to estimate the scaling factors.
*/	
void Simulation::getScalingApproximation(float &bScale, float &bOffset, float &wScale)
{
	const Float2D &h = block->getWaterHeight();
	const Float2D &b = block->getBathymetry();

	// Minimum values
	float minB, minH;
	// Maximum values
	float maxB, maxH;
        
        int nx = h.getRows()-2;
        int ny = h.getCols()-2;
	int maxDim = (nx > ny) ? nx : ny;

	minB = b[1][1];
	minH = h[1][1];
	maxB = b[1][1];
	maxH = h[1][1];

	for(int i=1; i<=nx; i++) {
		for(int j=1; j<=ny; j++) {
			// Update minima
			if ((h[i][j] + b[i][j]) < minH)
				minH = (h[i][j] + b[i][j]);
			if (b[i][j] < minB)
				minB = b[i][j];
			// Update maxima
			if ((h[i][j] + b[i][j]) > maxH)
				maxH = (h[i][j] + b[i][j]);
			if (b[i][j] > maxB)
				maxB = b[i][j];			
		}
	}
	bOffset = 0;	// This should be !=0 only in some artificial scenarios
	bScale = -80/minB;
	cout << "Scaling of bathymetry: " << bScale << endl;

	if ((maxH - minH) < 0.0001f) {
		wScale = 1.0f/(maxH- minH);
	} else {
		wScale = 1.0f;
	}
	wScale = (maxDim/50.0)*wScale;
	cout << "Scaling of water level: " << wScale << endl;

}

/**
    Compute normal of a triangle
	@param fVert1-fVert3	vertices of the triangle
	@param fNormal			resulting normal
*/	
void Simulation::calculateNormal(float fVert1[], float fVert2[],
                                 float fVert3[], float fNormal[]) {
   float Qx, Qy, Qz, Px, Py, Pz;

   Qx = fVert2[0]-fVert1[0];
   Qy = fVert2[1]-fVert1[1];
   Qz = fVert2[2]-fVert1[2];
   Px = fVert3[0]-fVert1[0];
   Py = fVert3[1]-fVert1[1];
   Pz = fVert3[2]-fVert1[2];

   fNormal[0] = Py*Qz - Pz*Qy;
   fNormal[1] = Pz*Qx - Px*Qz;
   fNormal[2] = Px*Qy - Py*Qx;
}

//==================================================================
// member functions for visualisation and output
//==================================================================
/**
    Scale function used for visualization 
*/	
__device__ __host__ 
float scaleFunction(float val, float scale) {
	return  val * scale;
}

// /**
//     Return coordinates of visBuffer[x][y]
// 	@param x,y		coordinates
// 	@param ny		length in y-direction
// */	
// __device__
// int getVisBufCoord(int x, int y, int ny) {
// 	return y*(ny+1) + x;
// }

/**
	Transform cell-centered values into node-centered values
	for visualization purposes
	
	@param visBuffer		pointer to visualization buffer
	@param hd				h-data on device
	@param bd				b-data on device
	@param nx, ny			x- and y-dimensions
	@param dx, dy			cellsize in x- and y-direction
	@param average			center (average) of h-data (water surface) 
							for visualization
	@param scale			scaling factor which is applied for visualization
	@param offset			offset factor for visualization
 */
__global__
void kernelCalcVisBuffer(float3* visBuffer, const float* hd, const float* bd, 
                         int nx, int ny)
{
   int i = TILE_SIZE*blockIdx.x + threadIdx.x;
   int j = TILE_SIZE*blockIdx.y + threadIdx.y;
   if ((i <= nx) && (j <= ny)) {
	   int index = i*(nx+2) + j;
	   int index2 = (i+1)*(nx+2) + j;	
	   if (hd[index] <= dryTol
			   || hd[index+1] <= dryTol
			   || hd[index2] <= dryTol
			   || hd[index2+1] <= dryTol)
		   visBuffer[j*(ny+1)+i] = make_float3(i, 0, j);
   	   else
   		   visBuffer[j*(ny+1) + i] = make_float3(
   				   i,
   				   0.25 * (
   						   hd[index]+hd[index + 1]+ hd[index2] + hd[index2 + 1] +
   						   bd[index]+bd[index + 1]+ bd[index2] + bd[index2 + 1]),
				   j);
   } 
   //Corresponding C-Code:
   //	for (int j=0; j<ny+1;j++)
   //		for (int i=0;i<nx+1;i++)
   //			Vtk_file << i*dx<<" "<<j*dy<<" "
   //			         << 0.25*(h[i][j]+h[i+1][j]+h[i][j+1]+h[i+1][j+1]
   //				         +b[i][j]+b[i+1][j]+b[i][j+1]+b[i+1][j+1]) 
}

/**
	copy h-values into visualization buffer
	@param _visBuffer		visualization buffer

*/
void Simulation::updateVisBuffer(float3* _visBuffer) {
	
        const float* hd = block->getCUDA_waterHeight();
        const float* bd = block->getCUDA_bathymetry();
        int nx = block->getNx();
        int ny = block->getNy();
//         // Fill ghost layer corner cells
// 	kernelHdBufferEdges<<<1,1>>>(hd, nx, ny);
	// Interpolate cell centered h-values
	dim3 dimBlock(TILE_SIZE,TILE_SIZE);
	dim3 dimGrid((nx+TILE_SIZE)/TILE_SIZE,(ny+TILE_SIZE)/TILE_SIZE);
	kernelCalcVisBuffer<<<dimGrid,dimBlock>>>(_visBuffer, hd, bd, nx, ny);
}
/**
	Function used for debugging. Outputs the current visBuffer 
	@param visBuffer		current visualization buffer
	
*/
void Simulation::debugVisBuffer(float3* _visBuffer) {
        int nx = block->getNx();
        int ny = block->getNy();
	float* dbgBuf = new float[(nx+1)*(ny+1)*3];
	int size = (nx+1)*(ny+1)*sizeof(float3);
	cudaMemcpy(dbgBuf,_visBuffer, size, cudaMemcpyDeviceToHost);
    checkCUDAError("memory of visualization buffer not transferred");
	for (int i = 0; i < 3*(nx+1)*(ny+1); i = i + 3) {
		std::cout << "(" << dbgBuf[i] << ", " << dbgBuf[i+1] << ", " << dbgBuf[i+2] << ") \n";
	}
	delete[] dbgBuf;
}

/**
    Compute normal of a triangle
	@param fVert1-fVert3	vertices of the triangle
	@return	resulting normal
*/	
inline __device__
float3 calculateNormal(float3 fVert1, float3 fVert2,
                             float3 fVert3)
   {
   float Qx, Qy, Qz, Px, Py, Pz;

   Qx = fVert2.x-fVert1.x;
   Qy = fVert2.y-fVert1.y;
   Qz = fVert2.z-fVert1.z;
   Px = fVert3.x-fVert1.x;
   Py = fVert3.y-fVert1.y;
   Pz = fVert3.z-fVert1.z;

   return make_float3(Py*Qz - Pz*Qy, Pz*Qx - Px*Qz, Px*Qy - Py*Qx);
}

/**
    Compute new normals resulting from updated water surface
	@param visBuffer			vertex buffer holding water surface vertices
	@param normalBuffer			buffer which will contain new normals
	@param nx, ny				
*/	
__global__
void kernelCalcNormals(float3* visBuffer, float3* normalBuffer, int nx, int ny)
{
   int i = TILE_SIZE*blockIdx.x + threadIdx.x;
   int j = TILE_SIZE*blockIdx.y + threadIdx.y;
   if ((i < nx + 1) && (j < ny + 1)) {
	   int _i = i;
	   int _j = j;
	   // Handle boundaries
	   if (i >= nx) _i = i - 1;
	   if (j >= ny) _j = j - 1;
	   normalBuffer[j*(ny+1) + i] = calculateNormal(visBuffer[(_j*(ny+1) + _i)], 
				visBuffer[(_j+1)*(ny+1) + _i + 1],
				visBuffer[(_j+1)*(ny+1) + _i]);
   } 
}

/**
    Compute normals of the new water surface on graphics card

    @param vertexBuffer			cuda resource holding the vertex positions
	@param destBuffer			cuda resource holding the new normals
*/
void Simulation::calculateNormals(float3* vertexBuffer, float3* destBuffer) {
// 	splash->calculateNormals(vertexBuffer, destBuffer);
	int nx = block->getNx();
	int ny = block->getNy();
	dim3 dimBlock(TILE_SIZE,TILE_SIZE);
	dim3 dimGrid((nx+TILE_SIZE)/TILE_SIZE,(ny+TILE_SIZE)/TILE_SIZE);
	kernelCalcNormals<<<dimGrid,dimBlock>>>(vertexBuffer, destBuffer, nx, ny);
}

/**
	Helper function to check for errors in CUDA calls
	param msg			custom error message
	source: NVIDIA

 */
