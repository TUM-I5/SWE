/**
 * @file
 * This file is part of SWE.
 *
 * @author Michael Bader, Kaveh Rahnema, Tobias Schnabel
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

#include <math.h>
#include "tools/help.hh"
#include "SWE_BlockCUDA.hh"
#include "SWE_RusanovBlockCUDA.hh"

//const int TILE_SIZE=16;
//const int TILE_SIZE=8;

// #include "SWE_BlockCUDA_kernels.cu"
#include "SWE_RusanovBlockCUDA_kernels.hh"

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
 */
SWE_RusanovBlockCUDA::SWE_RusanovBlockCUDA(float _offsetX, float _offsetY, const int i_cudaDevice)
 : SWE_BlockCUDA(_offsetX,_offsetY, i_cudaDevice)
#ifdef DBG
 , Fh(nx+1,ny+1), Fhu(nx+1,ny+1), Fhv(nx+1,ny+1),
   Gh(nx+1,ny+1), Ghu(nx+1,ny+1), Ghv(nx+1,ny+1)
#endif
{

  int size = (nx+1)*(ny+1)*sizeof(float);
  // allocate CUDA memory for flow unknows
  cudaMalloc((void**)&Fhd, size); 
     checkCUDAError("allocate device memory for Fh");
  cudaMalloc((void**)&Fhud, size);
     checkCUDAError("allocate device memory for Fhu");
  cudaMalloc((void**)&Fhvd, size);
     checkCUDAError("allocate device memory for Fhv");

  // allocate CUDA memory for flow unknows
  cudaMalloc((void**)&Ghd, size); 
     checkCUDAError("allocate device memory for Fh");
  cudaMalloc((void**)&Ghud, size);
     checkCUDAError("allocate device memory for Ghu");
  cudaMalloc((void**)&Ghvd, size);
     checkCUDAError("allocate device memory for Ghv");

  // allocate CUDA unknowns: bathymetry source terms 
  size = nx*ny*sizeof(float);
  cudaMalloc((void**)&Bxd, size); 
     checkCUDAError("allocate device memory for Bxd");
  cudaMalloc((void**)&Byd, size);
     checkCUDAError("allocate device memory for Byd");

  // allocate CUDA unknowns: maxmimum height and velocity
  size = (nx/TILE_SIZE)*(ny/TILE_SIZE)*sizeof(float);
  cudaMalloc((void**)&maxhd, size);
     checkCUDAError("allocate device memory for maxhd");
  cudaMalloc((void**)&maxvd, size);
     checkCUDAError("allocate device memory for maxvd");

}

/**
 * Destructor: de-allocate all variables
 */
SWE_RusanovBlockCUDA::~SWE_RusanovBlockCUDA() {
  cudaFree(Fhd); cudaFree(Fhud); cudaFree(Fhvd);
  cudaFree(Ghd); cudaFree(Ghud); cudaFree(Ghvd);
  cudaFree(Bxd); cudaFree(Byd);
  cudaFree(maxhd); cudaFree(maxvd);
  
}


//==================================================================
// methods for simulation
//==================================================================

/**
 * compute the largest allowed time step for the current grid block;
 * will be stored in the member variable SWE_Block::maxTimestep.
 * This specific CUDA implementation requires velocity data precomputed
 * by the kernel kernelEulerTimestep, and therefore needs to be called
 * directly after calling this kernel.
 */
void SWE_RusanovBlockCUDA::computeMaxTimestepCUDA() {

  float hmax = 0.0;
  float vmax = 0.0;
  float meshSize = (dx<dy) ? dx : dy;
   
  int numTiles = (nx/TILE_SIZE)*(ny/TILE_SIZE);
#ifdef DBG
  cout << "Number of tiles: " << numTiles << " -> ";
#endif
  // use 512 threads for maxmimum computation ...
  int threads = 512;
  // ... or the largest power-of-two smaller than numTiles:
  while (threads > numTiles) threads >>=1;
#ifdef DBG
  cout << "use " << threads << " threads per tile" << endl << flush;
#endif
  
  // loop over arrays of maxhd and maxvd to determine maximum
  // GPU maximum kernel called for 512 (or threads) elements at a time
  // CPU maintains maximum of all groups-of-512
  for(int tile=numTiles; tile>0; tile -= threads) {
     float newmax = 0;
     dim3 dimBlock(threads,1,1);
     dim3 dimGrid(1,1);

     // Compute maximum of arrays maxhd and maxvd
     // (for elements in range [tile-threads .. tile-1])
     int start = tile-threads;
     // use element range [0 .. threads-1], if tile-threads < 0
     if (start<0) start=0; 
//   cout << "Call kernel to compute time step " << start << flush << endl;
     kernelMaximum<<<dimGrid,dimBlock>>>(maxhd,maxvd,start,threads);
     checkCUDAError("compute maximum height and velocity");

     cudaMemcpy(&newmax,maxhd, sizeof(float), cudaMemcpyDeviceToHost);
  checkCUDAError("memory of hmax not transferred");
     if (newmax>hmax) hmax = newmax;
     cudaMemcpy(&newmax,maxvd, sizeof(float), cudaMemcpyDeviceToHost);
  checkCUDAError("memory of vmax not transferred");
     if (newmax>vmax) vmax = newmax;
     
#ifdef DBG
     cout << "Current maxmimum of tile " << start << ":"
          << hmax << ',' << vmax << endl << flush;
#endif
  };

#ifdef DBG
  cout << "Computed Time step (CUDA): " 
     << " from " << hmax << " and " << vmax << " -> "
     << meshSize/(sqrt(g*hmax) + vmax) << endl << flush;
#endif
  
  // sqrt(g*hmax) + vmax is the velocity of a characteristic shallow-water wave
  // such that a wave must not propagate farther than dx in a single time step
  // (uses 0.5 as a pessimistic factor to reduce time step size)
  maxTimestep = 0.5*meshSize/(sqrt(g*hmax) + vmax);
}


/**
 * Depending on the current values of h, hu, hv (incl. ghost layers)
 * update these unknowns in each grid cell
 * (ghost layers and bathymetry are not updated).
 * The Rusanov CUDA-implementation of simulateTimestep
 * subsequently calls the functions computeNumericalFluxes (to compute
 * all fluxes on grid edges), and updateUnknowns (to update the variables
 * according to flux values, typically according to an Euler time step).
 * @param dt  size of the time step
 */
void SWE_RusanovBlockCUDA::simulateTimestep(float dt) {

  // compute the numerical fluxes on all edges
  computeNumericalFluxes();

  // update the unknowns according to an Euler timestep
  updateUnknowns(dt);

}

/**
 * perform forward-Euler time steps, starting with simulation time tStart,:
 * until simulation time tEnd is reached; 
 * device-global variables hd, hud, hvd are updated;
 * unknowns h, hu, hv in main memory are not updated.
 * Ghost layers and bathymetry sources are updated between timesteps.
 * intended as main simulation loop between two checkpoints
 */
__host__
float SWE_RusanovBlockCUDA::simulate(float tStart, float tEnd) {
  float t = tStart;
  do {
     // set values in ghost cells:
     setGhostLayer();
     
     // execute Euler time step:
     simulateTimestep(maxTimestep);

     t += maxTimestep;

#ifdef DBG
     cout << "Simulation at time " << t << endl << flush;
#endif
  } while(t < tEnd);

  return t;
}

/**
 * implements interface function updateUnknowns:
 * based on the (Rusanov) fluxes computed on each edge
 * (and stored in the variables Fh, Gh, etc.);
 * compute the balance terms for each cell, and update the
 * unknowns according to an Euler time step.
 * It will force an update of the copy layer in the main
 * memory by calling synchCopyLayerBeforeRead(),
 * and provide an compute the maximum allowed time step size
 * by calling computeMaxTimestepCUDA().
 * @param dt  size of the time step.
 */
__host__
void SWE_RusanovBlockCUDA::updateUnknowns(float dt) {

  dim3 dimBlock(TILE_SIZE,TILE_SIZE);
  dim3 dimGrid(nx/TILE_SIZE,ny/TILE_SIZE);
#ifdef DBG
cout << "Call kernel for Euler timestep " << flush << endl;
#endif
  kernelEulerTimestep<<<dimGrid,dimBlock>>>(hd,hud,hvd,
                                            Fhd,Fhud,Fhvd,Ghd,Ghud,Ghvd,
              Bxd,Byd,
              maxhd,maxvd,
              nx,ny,dt,1.0f/dx,1.0f/dy);


  computeMaxTimestepCUDA();

  synchCopyLayerBeforeRead();
}

/**
 * compute the flux terms on all edges
 */
void SWE_RusanovBlockCUDA::computeNumericalFluxes() {

  // time step required to compute Lax Friedrichs constant dx/dt
  float dt = 2*maxTimestep;

  // compute bathymetry source terms (depend on h)
  computeBathymetrySources();

  // compute F-fluxes, i.e. all fluxes in x-direction (on left/right edges)
  dim3 dimBlock(TILE_SIZE,TILE_SIZE);
  dim3 dimGrid(nx/TILE_SIZE,ny/TILE_SIZE);
#ifdef DBG
cout << "Call kernel to compute F-fluxes " << flush << endl;
#endif
  kernelComputeFluxesF<<<dimGrid,dimBlock>>>(hd,hud,hvd, Fhd,Fhud,Fhvd,
            ny,g,dx/dt,1);

  dim3 dimLeftBlock(1,TILE_SIZE);
  dim3 dimLeftGrid(1,ny/TILE_SIZE);
  kernelComputeFluxesF<<<dimLeftGrid,dimLeftBlock>>>(hd,hud,hvd, Fhd,Fhud,Fhvd,
                 ny,g,dx/dt,0);
  
  // compute G-fluxes, i.e. all fluxes in y-direction (on top/bottom edges)
#ifdef DBG
cout << "Call kernel to compute G-fluxes " << flush << endl;
#endif
  kernelComputeFluxesG<<<dimGrid,dimBlock>>>(hd,hud,hvd, Ghd,Ghud,Ghvd,
            ny,g,dy/dt,1);

  dim3 dimTopBlock(TILE_SIZE,1);
  dim3 dimTopGrid(nx/TILE_SIZE,1);
  kernelComputeFluxesG<<<dimTopGrid,dimTopBlock>>>(hd,hud,hvd, Ghd,Ghud,Ghvd,
               ny,g,dy/dt,0);

// cout << (*this) << flush;

}

/**
 * compute the bathymetry source terms in all cells
 */
void SWE_RusanovBlockCUDA::computeBathymetrySources() {

  dim3 dimBlock(TILE_SIZE,TILE_SIZE);
  dim3 dimGrid(nx/TILE_SIZE,ny/TILE_SIZE);
#ifdef DBG
cout << "Call kernel to compute bathymetry sources" << flush << endl;
#endif
  kernelComputeBathymetrySources<<<dimGrid,dimBlock>>>(hd,bd,Bxd,Byd,ny,g);

// cout << (*this) << flush;

}


//==================================================================

/**
 * overload operator<< such that data can be written via cout <<
 * -> needs to be declared as friend to be allowed to access private data
 */
ostream& operator<<(ostream& os, const SWE_RusanovBlockCUDA& swe) {

  os << "Grid dimensions: " << swe.nx << "x" << swe.ny << endl;

  cout << "Water height:" << endl;
  for(int j=swe.ny+1; j>=0; j--) {
    for(int i=0; i<=swe.nx+1; i++) {
      os << swe.h[i][j] << "  ";
    };
    os << endl;
  };

  cout << "Momentum in x-direction:" << endl;
  for(int j=swe.ny+1; j>=0; j--) {
    for(int i=0; i<=swe.nx+1; i++) {
      os << swe.hu[i][j] << "  ";
    };
    os << endl;
  };

  cout << "Momentum in y-direction:" << endl;
  for(int j=swe.ny+1; j>=0; j--) {
    for(int i=0; i<=swe.nx+1; i++) {
      os << swe.hv[i][j] << "  ";
    };
    os << endl;
  };

#ifdef DBG
  cout << "Ghost/Copy Layer bottom:" << endl;
     for(int i=0; i<=swe.nx+1; i++) {
       os << swe.bottomGhostLayer->h[i]  << "  "; 
       os << swe.bottomGhostLayer->hu[i] << "  ";
       os << swe.bottomGhostLayer->hv[i] << "  ";
       os << swe.bottomCopyLayer->h[i]  << "  ";  
       os << swe.bottomCopyLayer->hu[i] << "  "; 
       os << swe.bottomCopyLayer->hv[i] << "  "; 
       cout << endl;
     };
  
  cout << "Ghost/Copy Layer top:" << endl;
     for(int i=0; i<=swe.nx+1; i++) {
       os << swe.topGhostLayer->h[i]  << "  "; 
       os << swe.topGhostLayer->hu[i] << "  ";
       os << swe.topGhostLayer->hv[i] << "  ";
       os << swe.topCopyLayer->h[i]  << "  ";  
       os << swe.topCopyLayer->hu[i] << "  "; 
       os << swe.topCopyLayer->hv[i] << "  "; 
       cout << endl;
     };
#endif


// #ifdef DBG
//   cout << "Fluss - Wellenhoehe:" << endl;
//   for(int i=0; i<=swe.nx; i++) {
//     for(int j=1; j<=swe.ny; j++) {
//       os << swe.Fh[i][j] << "  ";
//     };
//     os << endl;
//   };
// 
//   cout << "Fluss - Durchfluss in x-Richtung:" << endl;
//   for(int i=0; i<=swe.nx; i++) {
//     for(int j=1; j<=swe.ny; j++) {
//       os << swe.Fhu[i][j] << "  ";
//     };
//     os << endl;
//   };
// 
//   cout << "Fluss - Durchfluss in y-Richtung:" << endl;
//   for(int i=0; i<=swe.nx; i++) {
//     for(int j=1; j<=swe.ny; j++) {
//       os << swe.Fhv[i][j] << "  ";
//     };
//     os << endl;
//   };
// 
//   cout << "Fluss - Wellenhoehe:" << endl;
//   for(int i=1; i<=swe.nx; i++) {
//     for(int j=0; j<=swe.ny; j++) {
//       os << swe.Gh[i][j] << "  ";
//     };
//     os << endl;
//   };
// 
//   cout << "Fluss - Durchfluss in x-Richtung:" << endl;
//   for(int i=1; i<=swe.nx; i++) {
//     for(int j=0; j<=swe.ny; j++) {
//       os << swe.Ghu[i][j] << "  ";
//     };
//     os << endl;
//   };
// 
//   cout << "Fluss - Durchfluss in y-Richtung:" << endl;
//   for(int i=1; i<=swe.nx; i++) {
//     for(int j=0; j<=swe.ny; j++) {
//       os << swe.Ghv[i][j] << "  ";
//     };
//     os << endl;
//   };
// #endif
  
  os << flush;

  return os;
}

