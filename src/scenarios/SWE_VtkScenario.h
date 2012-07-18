#ifndef __SWE_VTK_SCENARIO_H
#define __SWE_VTK_SCENARIO_H
// =====================================================================
// This file is part of SWE_CUDA (see file SWE_Block.cu for details).
// 
// Copyright (C) 2010,2011 Michael Bader, Kaveh Rahnema, Tobias Schnabel
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

#include <iostream>
#include <stdio.h>
#include <fstream>

#include <string>
#include <vector>
#include "../tools/help.hh"
#include "SWE_Scenario.h"

using namespace std;

/**
 * Scenario "VTK File":
 * the object will read the values stored in a VTK file 
 * (i.e., use such a visualisation file as checkpoint for restarting a simulation)
 *
 * Note: - Boundary Conditions are NOT stored in a VTK file; 
 *         OUTFLOW conditions are assumed at all boundaries
 *       - similar, 0.1 is always returnedby function endSimulation()
 *
 * TODO: Option to read a VTK container file
 *       (for example: read by different MPI processes, such that each MPI process
 *                     will represent one of the container's VTK files)
 */
class SWE_VtkScenario : public SWE_Scenario {

  public:

    // SWE_VtkScenario does not offer a public constructor
    // -> use this method to create a SWE_VtkScenario object
    static SWE_VtkScenario* readVtkFile(char const * filename);

    virtual float getWaterHeight(float x, float y);
    virtual float getVeloc_u(float x, float y);
    virtual float getVeloc_v(float x, float y);
    virtual float getBathymetry(float x, float y);

    virtual float waterHeightAtRest() { return heightAtRest; };
    virtual void setWaterHeightAtRest(float whar) { heightAtRest = whar; };

    virtual BoundaryType getBoundaryType(BoundaryEdge edge) { return OUTFLOW; };
    virtual float getBoundaryPos(BoundaryEdge edge) { return bndPos[edge]; };

  protected:
  
    SWE_VtkScenario(int nx, int ny);
    
    Float2D h;
    Float2D u;
    Float2D v;
    Float2D b;
    
    float heightAtRest;
    float bndPos[4];
    float dx;
    float dy;

    void computeHeightAtRest();

    static void stringExplode(string str, string separator, vector<string>* results);

    friend ostream& operator<< (ostream& os, SWE_VtkScenario& scene);

};

ostream& operator<< (ostream& os, const SWE_VtkScenario& scene);



#endif
