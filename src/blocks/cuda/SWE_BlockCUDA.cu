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

#include "SWE_BlockCUDA.hh"
#include "SWE_BlockCUDA_kernels.hh"

#include "tools/help.hh"
#include "tools/Logger.hh"

#include <cassert>
#include <cstdlib>
#include <cmath>

//const int TILE_SIZE=16;
//const int TILE_SIZE=8;

// #include "SWE_BlockCUDA_kernels.cu"

/*
 * helper function to read CUDA error codes
 * (implementation in swe.cu */
void checkCUDAError(const char *msg)
{
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err)
    {
        fprintf(stderr, "\nCuda error (%s): %s.\n", msg, cudaGetErrorString( err) );
        exit(-1);
    }
}

/*
 * helper function to read CUDA error codes
 * (implementation in swe.cu */
void tryCUDA(cudaError_t err, const char *msg)
{
    if( cudaSuccess != err)
    {
        fprintf(stderr, "\nCuda error (%s): %s.\n", msg, cudaGetErrorString( err) );
        exit(-1);
    }
}

/**
 * Constructor: allocate variables for simulation
 *
 * unknowns h,hu,hv,b are defined on grid indices [0,..,nx+1]*[0,..,ny+1]
 * -> computational domain is [1,..,nx]*[1,..,ny]
 * -> plus ghost cell layer
 *
 * flux terms are defined for edges with indices [0,..,nx]*[1,..,ny]
 * or [1,..,nx]*[0,..,ny] (for horizontal/vertical edges)
 * Flux term with index (i,j) is located on the edge between 
 * cells with index (i,j) and (i+1,j) or (i,j+1)
 *
 * bathymetry source terms are defined for cells with indices [1,..,nx]*[1,..,ny]
 *
 *
 * @param i_cudaDevice ID of the CUDA-device, which should be used.
 */
SWE_BlockCUDA::SWE_BlockCUDA(
		int l_nx, int l_ny,
		float l_dx, float l_dy)
	: SWE_Block(l_nx, l_ny, l_dx, l_dy)
{
  if (nx % TILE_SIZE != 0) {
    cout << "WARNING: nx not a multiple of TILE_SIZE  -> will lead to crashes!" 
         << endl << flush;
  };
  if (ny % TILE_SIZE != 0) {
    cout << "WARNING: ny not a multiple of TILE_SIZE  -> will lead to crashes!" 
         << endl << flush;
  };

  int size = (nx+2)*(ny+2)*sizeof(float);
  // allocate CUDA memory for unknows h,hu,hv and bathymetry b
  cudaMalloc((void**)&hd, size); 
     checkCUDAError("allocate device memory for h");
  cudaMalloc((void**)&hud, size);
     checkCUDAError("allocate device memory for hu");
  cudaMalloc((void**)&hvd, size);
     checkCUDAError("allocate device memory for hv");
  cudaMalloc((void**)&bd, size); 
     checkCUDAError("allocate device memory for bd");

  // allocate consecutive memory for 2 columns with three unknowns each
  // (h, hu, hv, excluding b) for copy/ghost layer at bottom/top boundary
  size = nx+2;
  bottomLayer = new float[6*size];
  bottomGhostLayer = new SWE_Block1D( bottomLayer, bottomLayer+size, bottomLayer+(2*size), size );
  bottomCopyLayer  = new SWE_Block1D( bottomLayer+(3*size), bottomLayer+(4*size), bottomLayer+(5*size), size );
  
  // same for top boundary:
  topLayer = new float[6*size];
  topGhostLayer = new SWE_Block1D( topLayer, topLayer+size, topLayer+(2*size), size );
  topCopyLayer  = new SWE_Block1D( topLayer+(3*size), topLayer+(4*size), topLayer+(5*size), size );

  // allocate resp. memory on the CUDA device
  cudaMalloc((void**)&bottomLayerDevice, 6*size*sizeof(float));
     checkCUDAError("allocate device memory for bottom copy/ghost layer");
  cudaMalloc((void**)&topLayerDevice, 6*size*sizeof(float));
     checkCUDAError("allocate device memory for top copy/ghost layer");

}

/**
 * Destructor: de-allocate all variables
 */
SWE_BlockCUDA::~SWE_BlockCUDA() {
  cudaFree(hd); cudaFree(hud); cudaFree(hvd); cudaFree(bd);
//  cudaFree(maxhd); cudaFree(maxvd);
  
  // free memory for copy/ghost layers
  delete bottomLayer;
  delete topLayer;
  
  cudaFree(bottomLayerDevice);
  cudaFree(topLayerDevice);
}

//==================================================================
// methods for simulation
//==================================================================

/**
 * set the values of all ghost cells depending on the specifed 
 * boundary conditions
 */
void SWE_BlockCUDA::setBoundaryConditions() {
#ifdef DBG
 cout << "Call kernel to compute h in ghost layer corner (for visualisation only) " 
      << flush << endl;
#endif
  // Fill ghost layer corner cells
  kernelHdBufferEdges<<<1,1>>>(hd, nx, ny);

#ifdef DBG
 cout << "Call kernel to compute left/right boundaries " << flush << endl;
#endif
//   synchWaterHeightAfterWrite();
//   synchDischargeAfterWrite();

  if (boundary[BND_LEFT] == PASSIVE || boundary[BND_LEFT] == CONNECT) {
     // nothing to be done: 
     // ghost values are copied by SWE_BlockCUDA::synchGhostLayerAfterWrite(...)
  }
  else {
     dim3 dimBlock(1,TILE_SIZE);
     dim3 dimGrid(1,ny/TILE_SIZE);
     kernelLeftBoundary<<<dimGrid,dimBlock>>>(
        hd,hud,hvd,nx,ny,boundary[BND_LEFT]);
  };

  if (boundary[BND_RIGHT] == PASSIVE || boundary[BND_RIGHT] == CONNECT) {
     // nothing to be done: 
     // ghost values are copied by SWE_BlockCUDA::synchGhostLayerAfterWrite(...)
  }
  else {
     dim3 dimBlock(1,TILE_SIZE);
     dim3 dimGrid(1,ny/TILE_SIZE);
     kernelRightBoundary<<<dimGrid,dimBlock>>>(
        hd,hud,hvd,nx,ny,boundary[BND_RIGHT]);
  };

#ifdef DBG
  cout << "Call kernel to compute bottom/top boundaries " << flush << endl;
#endif
  switch (boundary[BND_BOTTOM]) {
    case CONNECT:
    {
      // set ghost layer data in auxiliary data structure for ghost layer:
      for(int i=0; i<=nx+1; i++) {
        bottomGhostLayer->h[i]  = neighbour[BND_BOTTOM]->h[i];
        bottomGhostLayer->hu[i] = neighbour[BND_BOTTOM]->hu[i];
        bottomGhostLayer->hv[i] = neighbour[BND_BOTTOM]->hv[i];
      };
    }
    case PASSIVE: /* also executed for case CONNECT */ 
    {
      // copy ghost layer data from buffer bottomLayerDevice
      // into bottom ghost layer of unknowns
      dim3 dimBlock(TILE_SIZE,1);
      dim3 dimGrid(nx/TILE_SIZE,1);
      kernelBottomGhostBoundary<<<dimGrid,dimBlock>>>(
         hd,hud,hvd,bottomLayerDevice,nx,ny);
      break;
    }
    default: 
    {
      // set simple boundary conditions (OUTFLOW, WALL) by resp. kernel:
      dim3 dimBlock(TILE_SIZE,1);
      dim3 dimGrid(nx/TILE_SIZE,1);
      kernelBottomBoundary<<<dimGrid,dimBlock>>>(
         hd,hud,hvd,nx,ny,boundary[BND_BOTTOM]);
    }
  };

  switch (boundary[BND_TOP]) {
    case CONNECT:
    {
      // set ghost layer data in auxiliary data structure for ghost layer:
      for(int i=0; i<=nx+1; i++) {
        topGhostLayer->h[i]  = neighbour[BND_TOP]->h[i];
        topGhostLayer->hu[i] = neighbour[BND_TOP]->hu[i];
        topGhostLayer->hv[i] = neighbour[BND_TOP]->hv[i];
      };
    }
    case PASSIVE: /* also executed for case CONNECT */
    {
      // copy ghost layer data from buffer bottomLayerDevice
      // into bottom ghost layer of unknowns
      dim3 dimBlock(TILE_SIZE,1);
      dim3 dimGrid(nx/TILE_SIZE,1);
      kernelTopGhostBoundary<<<dimGrid,dimBlock>>>(
         hd,hud,hvd,topLayerDevice,nx,ny);
      break;
    }
    default:
    { 
      // set simple boundary conditions (OUTFLOW, WALL) by resp. kernel: 
      dim3 dimBlock(TILE_SIZE,1);
      dim3 dimGrid(nx/TILE_SIZE,1);
      kernelTopBoundary<<<dimGrid,dimBlock>>>(
         hd,hud,hvd,nx,ny,boundary[BND_TOP]);
    }
  };
}

/**
 * synchronise the ghost layer content of h, hu, and hv in main memory 
 * with device memory and auxiliary data structures, i.e. transfer 
 * memory from main/auxiliary memory into device memory
 */
void SWE_BlockCUDA::synchGhostLayerAfterWrite() {

  if (boundary[BND_LEFT] == PASSIVE || boundary[BND_LEFT] == CONNECT) {
#ifdef DBG
  cout << "Set left passive boundary" << endl;
#endif
     // transfer h, hu, and hv from left ghost layer to resp. device memory
     cudaMemcpy(hd, h[0], (ny+2)*sizeof(float), cudaMemcpyHostToDevice);
     checkCUDAError("left ghost layer not transferred to device");
     cudaMemcpy(hud, hu[0], (ny+2)*sizeof(float), cudaMemcpyHostToDevice);
     checkCUDAError("left ghost layer not transferred to device");
     cudaMemcpy(hvd, hv[0], (ny+2)*sizeof(float), cudaMemcpyHostToDevice);
     checkCUDAError("left ghost layer not transferred to device");
  };
  if (boundary[BND_RIGHT] == PASSIVE || boundary[BND_RIGHT] == CONNECT) {
#ifdef DBG
  cout << "Set right passive boundary" << endl;
#endif
     // transfer h, hu, and hv from right ghost layer to resp. device memory
     cudaMemcpy(hd+((nx+1)*(ny+2)), h[nx+1], (ny+2)*sizeof(float), cudaMemcpyHostToDevice);
     checkCUDAError("right ghost layer not transferred to device");
     cudaMemcpy(hud+((nx+1)*(ny+2)), hu[nx+1], (ny+2)*sizeof(float), cudaMemcpyHostToDevice);
     checkCUDAError("right ghost layer not transferred to device");
     cudaMemcpy(hvd+((nx+1)*(ny+2)), hv[nx+1], (ny+2)*sizeof(float), cudaMemcpyHostToDevice);
     checkCUDAError("right ghost layer not transferred to device");
  };

  // bottom and top boundary ghost layers are replicated (for efficiency reasons)
  // in the  memory regions starting from bottomLayer and topLayer
  // -> these need to be transfered to device memory
  if (boundary[BND_BOTTOM] == PASSIVE || boundary[BND_BOTTOM] == CONNECT) {
#ifdef DBG
  cout << "Set bottom passive boundary" << endl;
  cout << " * transfer " << 3*(nx+2)*sizeof(float) << " bytes." << endl;
#endif
     // transfer bottom ghost layer 
     // (3 arrays - h, hu, hv - of nx+2 floats, consecutive in memory)
     cudaMemcpy(bottomLayerDevice, bottomLayer, 3*(nx+2)*sizeof(float), cudaMemcpyHostToDevice);
     checkCUDAError("bottom ghost layer not transferred to device");
  };

  if (boundary[BND_TOP] == PASSIVE || boundary[BND_TOP] == CONNECT) {
#ifdef DBG
  cout << "Set top passive boundary" << endl;
  cout << " * transfer " << 3*(nx+2)*sizeof(float) << " bytes." << endl;
#endif
     // transfer top ghost layer
     // (3 arrays - h, hu, hv - of nx+2 floats, consecutive in memory)
     cudaMemcpy(topLayerDevice, topLayer, 3*(nx+2)*sizeof(float), cudaMemcpyHostToDevice);
     checkCUDAError("top ghost layer not transferred to device");
  };

}

/**
 * Update (for heterogeneous computing) variables h, hu, hv in copy layers
 * before an external access to the unknowns 
 * (only required for PASSIVE and CONNECT boundaries)
 * - copy (up-to-date) content from device memory into main memory
 */
void SWE_BlockCUDA::synchCopyLayerBeforeRead() {
#ifdef DBG
  cout << "Transfer copy layers from device to main memory" << flush << endl;
#endif

  // copy values in copy layers to main memory
  // left copy layer:
  int offset = ny+2;
  cudaMemcpy(h[1], hd+offset, (ny+2)*sizeof(float), cudaMemcpyDeviceToHost);
  checkCUDAError("left ghost layer not transferred to device");
  cudaMemcpy(hu[1], hud+offset, (ny+2)*sizeof(float), cudaMemcpyDeviceToHost);
  checkCUDAError("left ghost layer not transferred to device");
  cudaMemcpy(hv[1], hvd+offset, (ny+2)*sizeof(float), cudaMemcpyDeviceToHost);
  checkCUDAError("left ghost layer not transferred to device");
  // right copy layer
  offset = nx*(ny+2);
  cudaMemcpy(h[nx], hd+offset, (ny+2)*sizeof(float), cudaMemcpyDeviceToHost);
  checkCUDAError("right ghost layer not transferred to device");
  cudaMemcpy(hu[nx], hud+offset, (ny+2)*sizeof(float), cudaMemcpyDeviceToHost);
  checkCUDAError("right ghost layer not transferred to device");
  cudaMemcpy(hv[nx], hvd+offset, (ny+2)*sizeof(float), cudaMemcpyDeviceToHost);
  checkCUDAError("right ghost layer not transferred to device");

  int size = 3*(nx+2);
  // bottom copy layer
  if (boundary[BND_BOTTOM] == PASSIVE || boundary[BND_BOTTOM] == CONNECT) {
     dim3 dimBlock(TILE_SIZE,1);
     dim3 dimGrid(nx/TILE_SIZE,1);
     kernelBottomCopyLayer<<<dimGrid,dimBlock>>>(
        hd,hud,hvd,bottomLayerDevice+size,nx,ny);

     cudaMemcpy(bottomLayer+size, bottomLayerDevice+size, size*sizeof(float), cudaMemcpyDeviceToHost);
     checkCUDAError("bottom copy layer not transferred from device");
  };

  // top copy layer
  if (boundary[BND_TOP] == PASSIVE || boundary[BND_TOP] == CONNECT) {
     dim3 dimBlock(TILE_SIZE,1);
     dim3 dimGrid(nx/TILE_SIZE,1);
     kernelTopCopyLayer<<<dimGrid,dimBlock>>>(
        hd,hud,hvd,topLayerDevice+size,nx,ny);
     
     cudaMemcpy(topLayer+size, topLayerDevice+size, size*sizeof(float), cudaMemcpyDeviceToHost);
     checkCUDAError("top copy layer not transferred from device");
  };

}

/**
 * register the row or column layer next to a boundary as a "copy layer",
 * from which values will be copied into the ghost layer or a neighbour;
 * @return	a SWE_Block1D object that contains row variables h, hu, and hv
 */
SWE_Block1D* SWE_BlockCUDA::registerCopyLayer(BoundaryEdge edge){

  // for TOP and BOTTOM layer, the implementation is identical to that in SWE_Block
  // for LEFT and RIGHT layer, separate layers are used that avoid strided copies 
  // when transferring memory between host and device memory
  switch (edge) {
    case BND_LEFT:
      return new SWE_Block1D( h.getColProxy(1), hu.getColProxy(1), hv.getColProxy(1));
    case BND_RIGHT:
      return new SWE_Block1D( h.getColProxy(nx), hu.getColProxy(nx), hv.getColProxy(nx));
    case BND_BOTTOM:
      // transfer bottom ghost and copy layer to extra SWE_Block1D
      for(int i=0; i<=nx+1; i++) {
        bottomGhostLayer->h[i]  = h[i][0];
        bottomGhostLayer->hu[i] = hu[i][0];
        bottomGhostLayer->hv[i] = hv[i][0];
        bottomCopyLayer->h[i]  = h[i][1];
        bottomCopyLayer->hu[i] = hu[i][1];
        bottomCopyLayer->hv[i] = hv[i][1];
      };
      return bottomCopyLayer;
    case BND_TOP:
      // transfer bottom ghost and copy layer to extra SWE_Block1D
      for(int i=0; i<=nx+1; i++) {
        topGhostLayer->h[i]  = h[i][ny+1];
        topGhostLayer->hu[i] = hu[i][ny+1];
        topGhostLayer->hv[i] = hv[i][ny+1];
        topCopyLayer->h[i]  = h[i][ny];
        topCopyLayer->hu[i] = hu[i][ny];
        topCopyLayer->hv[i] = hv[i][ny];
      };
      return topCopyLayer;
  };
  return NULL;
}

/**
 * "grab" the ghost layer at the specific boundary in order to set boundary values 
 * in this ghost layer externally. 
 * The boundary conditions at the respective ghost layer is set to PASSIVE, 
 * such that the grabbing program component is responsible to provide correct 
 * values in the ghost layer, for example by receiving data from a remote 
 * copy layer via MPI communication. 
 * @param	specified edge
 * @return	a SWE_Block1D object that contains row variables h, hu, and hv
 */
SWE_Block1D* SWE_BlockCUDA::grabGhostLayer(BoundaryEdge edge){

  // for TOP and BOTTOM layer, the implementation is identical to that in SWE_Block
  // for LEFT and RIGHT layer, separate layers are used that avoid strided copies 
  // when transferring memory between host and device memory
  boundary[edge] = PASSIVE;
  switch (edge) {
    case BND_LEFT:
      return new SWE_Block1D( h.getColProxy(0), hu.getColProxy(0), hv.getColProxy(0));
    case BND_RIGHT:
      return new SWE_Block1D( h.getColProxy(nx+1), hu.getColProxy(nx+1), hv.getColProxy(nx+1));
    case BND_BOTTOM:
      return bottomGhostLayer;
    case BND_TOP:
      return topGhostLayer;
  };
  return NULL;
}

/**
 * Print some available information about the CUDA devices.
 */
void SWE_BlockCUDA::printDeviceInformation()
{
	tools::Logger::logger.printString("Printing device information");

  //! id of the CUDA device.
  int l_deviceId;
  cudaGetDevice(&l_deviceId);

  //! total number of CUDA devices on this host.
  int l_deviceCount;
  cudaGetDeviceCount(&l_deviceCount);

  //! drive and runtime version
  int l_driverVersion, l_runtimeVersion;
  cudaDriverGetVersion(&l_driverVersion);
  cudaRuntimeGetVersion(&l_runtimeVersion);

  //! device properties
  cudaDeviceProp l_deviceProperty;
  cudaGetDeviceProperties(&l_deviceProperty, l_deviceId);

  // print information about the current device

  tools::Logger::logger.cout() << "Current CUDA device (relative to host): " << l_deviceId
                     << " ( " << l_deviceCount << " in total)" << std::endl;

  tools::Logger::logger.cout() << "CUDA device properties: "
                     << l_deviceProperty.name << " (name), "
                     << l_driverVersion << "/" << l_runtimeVersion << " (driver/runtime version), "
                     << l_deviceProperty.major << "." << l_deviceProperty.minor << " (compute capability)"
                     << std::endl;
}



//==================================================================
// protected member functions for memory model: 
// in case of temporary variables (especial in non-local memory, for 
// example on accelerators), the main variables h, hu, hv, and b 
// are not necessarily updated after each time step.
// The following methods are called to synchronise before or after 
// external read or write to the variables.
//==================================================================

/**
 * Update all temporary and non-local (for heterogeneous computing) variables
 * after an external update of the main variables h, hu, hv, and b.
 */
void SWE_BlockCUDA::synchAfterWrite() {
  // update h, hu, hv, b in device memory
  synchWaterHeightAfterWrite();
  synchDischargeAfterWrite();
  synchBathymetryAfterWrite();

  // update the auxiliary data structures for copy and ghost layers 
  // at bottom (and top, see below) boundaries
  // -> only required for PASSIVE and CONNECT boundaries
  if (boundary[BND_BOTTOM] == PASSIVE || boundary[BND_BOTTOM] == CONNECT) {
     // transfer bottom ghost and copy layer to extra SWE_Block1D
     for(int i=0; i<=nx+1; i++) {
       bottomGhostLayer->h[i]  = h[i][0];
       bottomGhostLayer->hu[i] = hu[i][0];
       bottomGhostLayer->hv[i] = hv[i][0];
       bottomCopyLayer->h[i]  = h[i][1];
       bottomCopyLayer->hu[i] = hu[i][1];
       bottomCopyLayer->hv[i] = hv[i][1];
     };
  };
  
  if (boundary[BND_TOP] == PASSIVE || boundary[BND_TOP] == CONNECT) {
     // transfer bottom ghost and copy layer to extra SWE_Block1D
     for(int i=0; i<=nx+1; i++) {
       topGhostLayer->h[i]  = h[i][ny+1];
       topGhostLayer->hu[i] = hu[i][ny+1];
       topGhostLayer->hv[i] = hv[i][ny+1];
       topCopyLayer->h[i]  = h[i][ny];
       topCopyLayer->hu[i] = hu[i][ny];
       topCopyLayer->hv[i] = hv[i][ny];
     };
  };

}

/**
 * Update temporary and non-local (for heterogeneous computing) variables
 * after an external update of the water height h
 */
void SWE_BlockCUDA::synchWaterHeightAfterWrite() {
#ifdef DBG
  cout << "Load water height h into device memory" << flush << endl;
#endif
  int size = (nx+2)*(ny+2)*sizeof(float);
  cudaMemcpy(hd,h.elemVector(), size, cudaMemcpyHostToDevice);
     checkCUDAError("memory of h not transferred");
}

/**
 * Update temporary and non-local (for heterogeneous computing) variables
 * after an external update of the discharge variables hu and hv
 */
void SWE_BlockCUDA::synchDischargeAfterWrite() {
#ifdef DBG
  cout << "Load discharge hu and hv into device memory" << flush << endl;
#endif
  int size = (nx+2)*(ny+2)*sizeof(float);
  cudaMemcpy(hud,hu.elemVector(), size, cudaMemcpyHostToDevice);
     checkCUDAError("memory of hu not transferred");
  cudaMemcpy(hvd,hv.elemVector(), size, cudaMemcpyHostToDevice);
     checkCUDAError("memory of hv not transferred");
}

/**
 * Update temporary and non-local (for heterogeneous computing) variables
 * after an external update of the bathymetry b
 */
void SWE_BlockCUDA::synchBathymetryAfterWrite() {
#ifdef DBG
  cout << "Load bathymetry unknowns into device memory" << flush << endl;
#endif
  int size = (nx+2)*(ny+2)*sizeof(float);
  cudaMemcpy(bd,b.elemVector(), size, cudaMemcpyHostToDevice);
     checkCUDAError("memory of b not transferred");
  
//  computeBathymetrySources();
}

/**
 * Update the main variables h, hu, hv, and b before an external read access:
 * copies current content of the respective device variables hd, hud, hvd, bd
 */
void SWE_BlockCUDA::synchBeforeRead() {
   synchWaterHeightBeforeRead();
   synchDischargeBeforeRead();
   synchBathymetryBeforeRead();

/* --- the following part is excluded:
   --- it should not be necessary to update the auxiliary data structures 
   --- for top/bottom copy/ghost layers in main memory: these need to be 
   --- kept consistent by class SWE_BlockCUDA (in particular, they cannot be 
   --- changed externally).
   --- */  
//   if (boundary[BND_BOTTOM] == PASSIVE || boundary[BND_BOTTOM] == CONNECT) {
//      // transfer bottom ghost and copy layer to extra SWE_Block1D
//      for(int i=0; i<=nx+1; i++) {
//        h[i][0] = bottomGhostLayer->h[i]; 
//        hu[i][0]= bottomGhostLayer->hu[i];
//        hv[i][0]= bottomGhostLayer->hv[i];
//        h[i][1] = bottomCopyLayer->h[i];  
//        hu[i][1]= bottomCopyLayer->hu[i]; 
//        hv[i][1]= bottomCopyLayer->hv[i]; 
//      };
//   };
//   
//   if (boundary[BND_TOP] == PASSIVE || boundary[BND_TOP] == CONNECT) {
//      // transfer bottom ghost and copy layer to extra SWE_Block1D
//      for(int i=0; i<=nx+1; i++) {
//        h[i][ny+1]  = topGhostLayer->h[i]; 
//        hu[i][ny+1] = topGhostLayer->hu[i];
//        hv[i][ny+1] = topGhostLayer->hv[i];
//        h[i][ny]  = topCopyLayer->h[i];  
//        hu[i][ny] = topCopyLayer->hu[i]; 
//        hv[i][ny] = topCopyLayer->hv[i]; 
//      };
//   };

}

/**
 * Update temporary and non-local (for heterogeneous computing) variables
 * before an external access to the water height h
 */
void SWE_BlockCUDA::synchWaterHeightBeforeRead() {
  int size = (nx+2)*(ny+2)*sizeof(float);
#ifdef DBG
  cout << "Copy water height h from device" << flush << endl;
#endif
  cudaMemcpy(h.elemVector(),hd, size, cudaMemcpyDeviceToHost);
     checkCUDAError("memory of h not transferred");

#ifdef DBG
  cout << "Set water height in ghost-layer corner cells" << flush << endl;
#endif
  // only required for visualisation: set values in corner ghost cells
  h[0][0] = h[1][1];
  h[0][ny+1] = h[1][ny];
  h[nx+1][0] = h[nx][1];
  h[nx+1][ny+1] = h[nx][ny];
}

/**
 * Update temporary and non-local (for heterogeneous computing) variables
 * before an external access to the discharge variables hu and hv
 */
void SWE_BlockCUDA::synchDischargeBeforeRead() {
  int size = (nx+2)*(ny+2)*sizeof(float);
#ifdef DBG
  cout << "Copy discharge hu and hv from device" << flush << endl;
#endif
  cudaMemcpy(hu.elemVector(),hud, size, cudaMemcpyDeviceToHost);
     checkCUDAError("memory of hu not transferred");
  cudaMemcpy(hv.elemVector(),hvd, size, cudaMemcpyDeviceToHost);
     checkCUDAError("memory of hv not transferred");

}

/**
 * Update temporary and non-local (for heterogeneous computing) variables
 * before an external access to the bathymetry b
 */
void SWE_BlockCUDA::synchBathymetryBeforeRead() {
  int size = (nx+2)*(ny+2)*sizeof(float);
#ifdef DBG
  cout << "Copy water bathymetry b from device" << flush << endl;
#endif
  cudaMemcpy(b.elemVector(),bd, size, cudaMemcpyDeviceToHost);
     checkCUDAError("memory of b not transferred");
}

void SWE_BlockCUDA::init(int i_cudaDevice)
{
	  tools::Logger::logger.setProcessRank(i_cudaDevice);

	  cudaSetDevice(i_cudaDevice);

	  // check for a valid CUDA device id
	  #ifndef NDEBUG
	  int l_deviceCount;
	  cudaGetDeviceCount(&l_deviceCount);
	  assert( (i_cudaDevice >= 0) && (i_cudaDevice < l_deviceCount) );
	  #endif

	  printDeviceInformation();

	  // Make sure the cuda device is reset at exit
	  atexit( SWE_BlockCUDA::finalize );
}

void SWE_BlockCUDA::finalize()
{
	// reset the cuda device
	tools::Logger::logger.printString("Resetting the CUDA devices");
	cudaDeviceReset();
}
