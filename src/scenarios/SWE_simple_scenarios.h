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

#ifndef __SWE_SIMPLE_SCENARIOS_H
#define __SWE_SIMPLE_SCENARIOS_H

#include <math.h>

#include "SWE_Scenario.h"

/**
 * Scenario "Radial Dam Break":
 * elevated water in the center of the domain
 */
class SWE_RadialDamBreakScenario : public SWE_Scenario {

  public:

    float getWaterHeight(float x, float y) { 
       return ( sqrt( (x-0.5f)*(x-0.5f) + (y-0.5f)*(y-0.5f) ) < 0.1f ) ? 12.0f: 10.0f;
    };

    virtual BoundaryType getBoundaryType(BoundaryEdge edge) { return OUTFLOW; };
};

/**
 * Scenario "Bathymetry Dam Break":
 * uniform water depth, but elevated bathymetry in the centre of the domain
 */
class SWE_BathymetryDamBreakScenario : public SWE_Scenario {

  public:

    float getBathymetry(float x, float y) { 
       // return ( sqrt( (x-0.3f)*(x-0.3f) + (y-0.8f)*(y-0.8f) ) < 0.1f ) ? 0.1f: 0.0f;
       return ( sqrt( (x-0.5f)*(x-0.5f) + (y-0.5f)*(y-0.5f) ) < 0.1f ) ? 0.1f: 0.0f;
    };
    virtual float endSimulation() { return 0.2f; };

    virtual BoundaryType getBoundaryType(BoundaryEdge edge) { return OUTFLOW; };
};

/**
 * Scenario "Sea at Rest":
 * flat water surface ("sea at rest"), 
 * but non-uniform bathymetry (id. to "Bathymetry Dam Break")
 * test scenario for "sea at rest"-solution 
 */
class SWE_SeaAtRestScenario : public SWE_Scenario {

  public:

    float getWaterHeight(float x, float y) { 
       return ( sqrt( (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) ) < 0.1f ) ? 9.9f: 10.0f;
    };
    float getBathymetry(float x, float y) { 
       return ( sqrt( (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) ) < 0.1f ) ? 0.1f: 0.0f;
    };

};

/**
 * Scenario "Splashing Pool":
 * intial water surface has a fixed slope (diagonal to x,y)
 */
class SWE_SplashingPoolScenario : public SWE_Scenario {

  public:

    float getWaterHeight(float x, float y) { return 10.0f+(1.0f-(x+y));};
    float endSimulation() { return 1.0f; };

};

/**
 * Scenario "Splashing Cone":
 * bathymetry forms a circular cone
 * intial water surface designed to form "sea at rest"
 * but: elevated water region in the centre (similar to radial dam break)
 */
class SWE_SplashingConeScenario : public SWE_Scenario {

  public:

    float getWaterHeight(float x, float y) { 
       float r = sqrt( (x-0.5f)*(x-0.5f) + (y-0.5f)*(y-0.5f) );
       float h = 4.0f-4.5f*(r/0.5f);

       if (r<0.1f) h = h+1.0f;

       return (h>0.0f) ? h : 0.0f;
    };

    float getBathymetry(float x, float y) { 
       float r = sqrt( (x-0.5f)*(x-0.5f) + (y-0.5f)*(y-0.5f) );
       return 1.0f+9.0f*( (r < 0.5f) ? r : 0.5f);
    };
    
    float waterHeightAtRest() { return 4.0f; };
    float endSimulation() { return 0.5f; };

    virtual BoundaryType getBoundaryType(BoundaryEdge edge) { return OUTFLOW; };
};

#endif
