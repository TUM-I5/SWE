/**
 * @file
 * This file is part of SWE.
 *
 * @author Michael Bader, Kaveh Rahnema, Tobias Schnabel
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
 *
 * @section LICENSE
 *
 * SWE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SWE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SWE.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * @section DESCRIPTION
 *
 * TODO
 */

#ifndef __SWE_BLOCKCUDA_HH
#define __SWE_BLOCKCUDA_HH

#include "blocks/SWE_Block.hh"

#include "tools/help.hh"

#include <iostream>
#include <fstream>

#include <cuda_runtime.h>

using namespace std;

void checkCUDAError(const char *msg);
void tryCUDA(cudaError_t err, const char *msg);

const int TILE_SIZE=16;
//const int TILE_SIZE=8;

/**
 * SWE_BlockCUDA extends the base class SWE_Block towards  
 * a base class for a CUDA implementation of the shallow water equations.
 * It adds the respective variables in GPU memory, and provides 
 * methods for data transfer between main and GPU memory.
 */
class SWE_BlockCUDA : public SWE_Block {

  public:
    // Constructor und Destructor
    SWE_BlockCUDA(int l_nx, int l_ny,
    		float l_dx, float l_dy);
    virtual ~SWE_BlockCUDA();
    
  // object methods

// ---> COULD BE IMPLEMENTED TO PROVIDE A DEFAULT IMPLEMENTATION
//     // determine maximum possible time step
//     virtual float getMaxTimestep();

    // deliver a pointer to proxy class that represents
    // the layer that is copied to an external ghost layer 
    virtual SWE_Block1D* registerCopyLayer(BoundaryEdge edge);
    // "grab" the ghost layer in order to set these values externally
    virtual SWE_Block1D* grabGhostLayer(BoundaryEdge edge);

    // access to CUDA variables
    /**
     *  @return	pointer to the array #hd (water height) in device memory 
     */
    const float* getCUDA_waterHeight() { return hd; };
    /**
     *  @return	pointer to the array #hb (bathymetry) in device memory 
     */
    const float* getCUDA_bathymetry() { return bd; };

  protected:
     
    // synchronisation Methods
    virtual void synchAfterWrite();
    virtual void synchWaterHeightAfterWrite();
    virtual void synchDischargeAfterWrite();
    virtual void synchBathymetryAfterWrite();
    virtual void synchGhostLayerAfterWrite();

    virtual void synchBeforeRead();
    virtual void synchWaterHeightBeforeRead();
    virtual void synchDischargeBeforeRead();
    virtual void synchBathymetryBeforeRead();
    virtual void synchCopyLayerBeforeRead();
    
    // set boundary conditions in ghost layers (set boundary conditions)
    virtual void setBoundaryConditions();

    // define arrays for main unknowns in CUDA global memory: 
    // hd, hud, hvd, and bd are CUDA arrays corresp. to h, hu, hv, and b
    float* hd;
    float* hud;
    float* hvd;
    float* bd;
	
  private:
     
    // separate memory to hold bottom and top ghost and copy layer 
    // in main memory allowing non-strided access
    float* bottomLayer;
    float* topLayer;
    SWE_Block1D* bottomGhostLayer;
    SWE_Block1D* bottomCopyLayer;
    SWE_Block1D* topGhostLayer;
    SWE_Block1D* topCopyLayer;
    // and resp. memory on the CUDA device:
    float* bottomLayerDevice;
    float* topLayerDevice;

    // helper arrays: store maximum height and velocities to determine time step
    float* maxhd;
    float* maxvd;

  public:
    // print information about the CUDA device
    static void printDeviceInformation();

    /**
     *  Initializes the cuda device
     *  Has to be called once at the beginning.
     */
    static void init(int i_cudaDevice = 0);
    /** Cleans up the cuda device */
    static void finalize();
};

/**
    Return index of hd[i][j] in linearised array
	@param i,j		x- and y-coordinate of grid cell
	@param ny		grid size in y-direction (without ghost layers)
*/	
inline __device__
int getCellCoord(int x, int y, int ny) {
   return x*(ny+2) + y;
}


/**
    Return index of edge-data Fhd[i][j] or Ghd[i][j] in linearised array
	@param i,j		x- and y-coordinate of grid cell
	@param ny		grid size in y-direction (without ghost layers)
*/	
inline __device__
int getEdgeCoord(int x, int y, int ny) {
   return x*(ny+1) + y;
}

/**
    Return index of a specific element in the arrays of bathymetry source terms
	@param i,j		x- and y-coordinate of grid cell
	@param ny		grid size in y-direction (without ghost layers)
*/	
inline __device__
int getBathyCoord(int x, int y, int ny) {
   return x*ny + y;
}


#endif
