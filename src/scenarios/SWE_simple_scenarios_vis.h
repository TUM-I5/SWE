#ifndef __SWE_SIMPLE_SCENARIOS_VIS_H
#define __SWE_SIMPLE_SCENARIOS_VIS_H
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

#include <math.h>

#include "SWE_simple_scenarios.h"
#include "SWE_VisInfo.h"

/**
 * Scenario "Radial Dam Break":
 * elevated water in the center of the domain
 */
class SWE_RadialDamBreakScenarioVisInfo 
: public SWE_RadialDamBreakScenario,
  public SWE_VisInfo {

  virtual float waterVerticalScaling() { return 5.0f; };
  /* use inherited default implementations */
};

/**
 * Scenario "Splashing Pool with custom scaling":
 * intial water surface has a fixed slope (diagonal to x,y)
 * shows how to use custom scaling to enhance visualization 
 * results
 */
class SWE_SplashingPoolScenarioVisInfo 
 : public SWE_SplashingPoolScenario, 
   public SWE_VisInfo {

  public:

    float waterDistanceFromGround() { return 9.0f; };
    float waterVerticalScaling() { return 5.0f; };
    float bathyVerticalCenter() { return 0.0f; };
    float bathyDistanceFromGround() { return 0.0f; };
    float bathyVerticalScaling() { return 0.0f; };

    virtual BoundaryType getBoundaryType(BoundaryEdge edge) { return OUTFLOW; };
};



#endif
