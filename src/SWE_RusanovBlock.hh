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

#ifndef __SWE_RUSANOVBLOCK_HH
#define __SWE_RUSANOVBLOCK_HH

#include <iostream>
#include <stdio.h>
#include <fstream>

#include "tools/help.hh"
#include "SWE_Block.hh"


using namespace std;

/**
 * SWE_RusanovBlock is an implementation of the SWE_Block abstract class. 
 * It uses a simple Rusanov flux (aka local Lax-Friedrich) in the model,
 * with some simple modifications to obtain a well-balanced scheme. 
 */
class SWE_RusanovBlock : public SWE_Block {

  public:
    // Constructor und Destructor
    SWE_RusanovBlock(float _offsetX = 0, float _offsetY = 0);
    virtual ~SWE_RusanovBlock();
    
  // object methods
    /// execute a single time step of the simulation
    virtual void simulateTimestep(float dt);
    /// compute simulate from specified start to end time
    virtual float simulate(float tStart, float tEnd);
    
    /// compute flux terms on edges
    virtual void computeNumericalFluxes();
    /// update unknowns according to fluxes (Euler time step)
    virtual void updateUnknowns(float dt);

  protected:
     
    /// compute source terms
    virtual void computeBathymetrySources();

    static float computeFlux(float fLoc, float fNeigh, float xiLoc, float xiNeigh, float llf);
    float computeLocalSV(int i, int j, char dir);

    // compute the largest allowed time step for the current grid block
    virtual void computeMaxTimestep() {
       SWE_Block::computeMaxTimestep();
       // more pessimistic choice of the time step
       maxTimestep *= 0.5; 
    };

    // define additional arrays for temporary unknowns: 
    // - arrays to hold the values of the flux terms at cell edges
    Float2D Fh;
    Float2D Fhu;
    Float2D Fhv;
    Float2D Gh;
    Float2D Ghu;
    Float2D Ghv;
    // - arrays to hold the bathymetry source terms for the hu and hv equations
    Float2D Bx;
    Float2D By;
    
    // overload operator<< such that data can be written via cout <<
    // -> needs to be declared as friend to be allowed to access private data
    friend ostream& operator<< (ostream& os, const SWE_RusanovBlock& swe);
  
};

ostream& operator<< (ostream& os, const SWE_RusanovBlock& swe);



#endif
