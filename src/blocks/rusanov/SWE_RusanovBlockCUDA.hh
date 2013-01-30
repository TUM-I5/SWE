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

#ifndef __SWE_RUSANOVBLOCKCUDA_HH
#define __SWE_RUSANOVBLOCKCUDA_HH

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <cuda_runtime.h>
#include "tools/help.hh"
#include "SWE_Block.hh"
#include "SWE_BlockCUDA.hh"

using namespace std;

/**
 * SWE_RusanovBlockCUDA extends the base class SWE_BlockCUDA, 
 * and provides a concrete CUDA implementation of a simple 
 * shallow water model based on Rusanov Flux computation on the 
 * edges and explicit time stepping.
 */
class SWE_RusanovBlockCUDA : public SWE_BlockCUDA {

  public:
    // Constructor und Destructor
    SWE_RusanovBlockCUDA(float _offsetX = 0, float _offsetY = 0, const int i_cudaDevice = 0);
    virtual ~SWE_RusanovBlockCUDA();
    
  // object methods

    virtual void computeNumericalFluxes();
    // simulate for specified time range
    // execute Euler time step
    virtual void updateUnknowns(float dt);
    /// execute a single time step of the simulation
    virtual void simulateTimestep(float dt);
    // compute flux terms on edges
    virtual float simulate(float tStart, float tEnd);
    
  private:
     
    // compute bathymetry source terms
    void computeBathymetrySources();

    // determine maximum possible time step
    void computeMaxTimestepCUDA();

    // arrays to hold the values of the flux terms at cell edges
    float* Fhd;
    float* Fhud;
    float* Fhvd;
    float* Ghd;
    float* Ghud;
    float* Ghvd;

    // arrays to hold the bathymetry source terms for the hu and hv equations
    float* Bxd;
    float* Byd;
    
    // helper arrays: store maximum height and velocities to determine time step
    float* maxhd;
    float* maxvd;

    // overload operator<< such that data can be written via cout <<
    // -> needs to be declared as friend to be allowed to access private data
    friend ostream& operator<< (ostream& os, const SWE_RusanovBlockCUDA& swe);

#ifdef DBG
    // --- only required for debugging purposes ---
    // arrays for fluxes for h,hu,hv in main memory
    Float2D Fh; 
    Float2D Fhu;
    Float2D Fhv;
    Float2D Gh; 
    Float2D Ghu;
    Float2D Ghv;
    // dump fluxes for h,hu,hv from CUDA device memory into main memory
    void cudaDumpFlux();
#endif
    
};

ostream& operator<< (ostream& os, const SWE_RusanovBlockCUDA& swe);

#endif
