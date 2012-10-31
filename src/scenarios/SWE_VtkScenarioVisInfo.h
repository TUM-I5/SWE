#ifndef __SWE_VTK_SCENARIO_VISINFO_H
#define __SWE_VTK_SCENARIO_VISINFO_H
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

#include "../tools/help.hh"
#include "SWE_VtkScenario.h"
#include "SWE_VisInfo.hh"

/**
 * Scenario "VTK File" incl. info for visualisation:
 * extends Scenario SWE_VtkScenario
 */
class SWE_VtkScenarioVisInfo : public SWE_VtkScenario, public SWE_VisInfo {

  public:

    // SWE_VtkScenario does not offer a public constructor
    // -> use this method to create a SWE_VtkScenario object
    static SWE_VtkScenarioVisInfo* readVtkFile(char const * filename);

    // Scaling factors for custom visualization
    virtual float waterHeightAtRest() { return hAverage; };
    virtual float waterDistanceFromGround() { return hOffset; };
    virtual float waterVerticalScaling() { return hScale; };
    virtual float bathyVerticalCenter() { return bAverage; };
    virtual float bathyDistanceFromGround() { return bOffset; };
    virtual float bathyVerticalScaling() { return bScale; };

  protected:
  
    SWE_VtkScenarioVisInfo(int nx, int ny);
    
    void computeScalingApproximation();
    
    float bAverage;
    float hAverage;
    float bScale;
    float hScale;
    float bOffset;
    float hOffset;
};



#endif
