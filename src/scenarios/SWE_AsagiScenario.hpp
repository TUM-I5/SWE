/**
 * @file
 * This file is part of SWE.
 *
 * @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
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
 * Access to bathymetry and displacement files with ASAGI.
 */

#ifndef ASAGI_HPP_
#define ASAGI_HPP_

#include <cassert>
#include <string>
#include <asagi.h>
#include "SWE_Scenario.h"

namespace scenarios {
  class Asagi;
}

class scenarios::Asagi: public SWE_Scenario {
//private:
    //! pointer to the Asagi bathymetry grid
    asagi::Grid* bathymetryGrid;

    //! pointer the the Asagi displacement grid
    asagi::Grid* displacementGrid;

    //! flag whether the displacement is dynamic or static
    const bool dynamicDisplacement;

    //! duration of the simulation
    const float duration;

#ifndef NDEBUG
    //! ranges of the bathymetry file, 0: min(x), 1: max(x), 2: min(y), 3: max(y)
    float bathymetryRange[4];
#endif

    //! ranges of the displacement file, 0: min(x), 1: max(x), 2: min(y), 3: max(y)
    float displacementRange[4];

    //! ranges of the time is a dynamic displacement is used 0: min(t), 1: max(t)
    float dynamicDisplacementTimeRange[2];

    //! ranges of the simulation area: min(x), 1: max(x), 2: min(y), 3: max(y)
   float simulationArea[4];

  public:
    /**
     * Constructor of an Asagi Scenario, which initializes the corresponding Asagi grids.
     *
     * @param i_originX origin of the simulation area (x-direction)
     * @param i_originY origin of the simulation area (y-direction)
     * @param i_bathymetryFile path to the netCDF-bathymetry file
     * @param i_displacementFile path to the netCDF-displacement file
     * @param i_duration time the simulation runs (in seconds)
     */
    Asagi ( const std::string i_bathymetryFile,
            const std::string i_displacementFile,
            const float i_duration,
            const float i_simulationArea[4],
            const bool i_dynamicDisplacement = false ):
              dynamicDisplacement(i_dynamicDisplacement),
              duration(i_duration) {
      //create the bathymetry grid
      bathymetryGrid = asagi::Grid::create( asagi::Grid::FLOAT );
      //create the displacement grid
      displacementGrid = asagi::Grid::create( asagi::Grid::FLOAT );

      int l_asagiOpen = bathymetryGrid->open(i_bathymetryFile.c_str());
      //open the bathymetry grid
      if( l_asagiOpen != 0 ) {
        std::cout << "Could not open bathymetry file: " << i_bathymetryFile << std::endl;
        std::cout << "Error code: " << l_asagiOpen << std::endl;
        assert(false);
      }

      l_asagiOpen = displacementGrid->open(i_displacementFile.c_str());
      //open the displacement grid
      if( l_asagiOpen != 0 ) {
        std::cout << "Could not open displacement file: " << i_displacementFile << std::endl;
        std::cout << "Error code: " << l_asagiOpen << std::endl;
        assert(false);
      }

#ifndef NDEBUG
      //read grid information
      bathymetryRange[0] = bathymetryGrid->getXMin();
      bathymetryRange[1] = bathymetryGrid->getXMax();
      bathymetryRange[2] = bathymetryGrid->getYMin();
      bathymetryRange[3] = bathymetryGrid->getYMax();
#endif

      displacementRange[0] = displacementGrid->getXMin();
      displacementRange[1] = displacementGrid->getXMax();
      displacementRange[2] = displacementGrid->getYMin();
      displacementRange[3] = displacementGrid->getYMax();
      if(dynamicDisplacement == false) {
        dynamicDisplacementTimeRange[0] = dynamicDisplacementTimeRange[1] = 0;
      }
      else {
        dynamicDisplacementTimeRange[0] = displacementGrid->getZMin();
        dynamicDisplacementTimeRange[1] = displacementGrid->getZMax();
      }

      simulationArea[0] = i_simulationArea[0];
      simulationArea[1] = i_simulationArea[1];
      simulationArea[2] = i_simulationArea[2];
      simulationArea[3] = i_simulationArea[3];

#ifndef NDEBUG
#ifdef PRINT_ASAGI_INFORMATION
      //print information
      std::cout << "  *** scenarios::Asagi created" << std::endl
                << "    i_bathymetryFile=" << i_bathymetryFile << std::endl
                << "    i_displacementFile=" << i_displacementFile << std::endl
                << "    duration= " << duration << std::endl
                << "    bathymetryRange[0]=" << bathymetryRange[0] << std::endl
                << "    bathymetryRange[1]=" << bathymetryRange[1] << std::endl
                << "    bathymetryRange[2]=" << bathymetryRange[2] << std::endl
                << "    bathymetryRange[3]=" << bathymetryRange[3] << std::endl
                << "    displacementRange[0]=" << displacementRange[0] << std::endl
                << "    displacementRange[1]=" << displacementRange[1] << std::endl
                << "    displacementRange[2]=" << displacementRange[2] << std::endl
                << "    displacementRange[3]=" << displacementRange[3] << std::endl
                << "    dynamicDisplacementTimeRange[0]=" << dynamicDisplacementTimeRange[0] << std::endl
                << "    dynamicDisplacementTimeRange[1]=" << dynamicDisplacementTimeRange[1] << std::endl
                << "    simulationArea[0]=" << simulationArea[0] << std::endl
                << "    simulationArea[1]=" << simulationArea[1] << std::endl
                << "    simulationArea[2]=" << simulationArea[2] << std::endl
                << "    simulationArea[3]=" << simulationArea[3] << std::endl;
#endif
#endif

    }

    virtual ~Asagi() {
    }

    void deleteGrids() {
      delete bathymetryGrid;
      delete displacementGrid;
    }

    //methods from SWE_SCENARIO

    /**
     * Get the water height at a specific location (before the initial displacement).
     *
     * @param i_positionX position relative to the origin of the bathymetry grid in x-direction
     * @param i_positionY position relative to the origin of the bathymetry grid in y-direction
     * @return water height (before the initial displacement)
     */
    float getWaterHeight( float i_positionX,
                          float i_positionY ) {
      assert(i_positionX > bathymetryRange[0]);
      assert(i_positionX < bathymetryRange[1]);
      assert(i_positionY > bathymetryRange[2]);
      assert(i_positionY < bathymetryRange[3]);

      float bathymetryValue = bathymetryGrid->getFloat2D(i_positionX, i_positionY);

      if( bathymetryValue > (float)0. ) {
        return 0.;
      }
      else {
        return -bathymetryValue;
      }
    }

    /**
     * Get the bathymetry including static displacement at a specific location
     *
     * @param i_positionX position relative to the origin of the displacement grid in x-direction
     * @param i_positionY position relative to the origin of the displacement grid in y-direction
     * @return bathymetry (after the initial displacement (static displacement)
     */
    float getBathymetry( const float i_positionX,
                         const float i_positionY ) {
      //assert that the 2D wrapper is not used for 3D displacements
      //assert(dynamicDisplacement == false);
      // no assertation for compability

      return getBathymetryAndDynamicDisplacement(i_positionX, i_positionY, 0);
    }

    /**
     * Get the bathymetry including dynamic displacement at a specific location
     *
     * @param i_positionX position relative to the origin of the displacement grid in x-direction
     * @param i_positionY position relative to the origin of the displacement grid in y-direction
     * @param i_time time relative to the origin of the dynamic displacement
     * @return bathymetry (after the initial displacement (static displacement), after the specified amount of time (dynamic displacement))
     */
    float getBathymetryAndDynamicDisplacement( const float i_positionX,
                                               const float i_positionY,
                                               const float i_time ) {
      assert(i_positionX > bathymetryRange[0]);
      assert(i_positionX < bathymetryRange[1]);
      assert(i_positionY > bathymetryRange[2]);
      assert(i_positionY < bathymetryRange[3]);

      float bathymetryValue = bathymetryGrid->getFloat2D(i_positionX, i_positionY);

      //bathymetryValue = (float) 0.; //TODO: remove: old file format

      float displacementValue = (float) 0.;

      if ( i_positionX > displacementRange[0] &&
           i_positionX < displacementRange[1] &&
           i_positionY > displacementRange[2] &&
           i_positionY < displacementRange[3] ) {
        if(dynamicDisplacement == false)
          displacementValue = displacementGrid->getFloat2D(i_positionX, i_positionY);
        else
          displacementValue = displacementGrid->getFloat3D(i_positionX, i_positionY, i_time);
      }

      return bathymetryValue + displacementValue;
    }

    /**
     * Check if there is an dynamic displacement is available for the corresponding time.
     * @param i_time current simulation time
     * @return true if there is data available, false else
     */
    bool dynamicDisplacementAvailable(const float i_time) {
      if( i_time > dynamicDisplacementTimeRange[0] &&
          i_time < dynamicDisplacementTimeRange[1] )
        return true;
      else
        return false;
    }

    /**
     * Get the number of seconds, the simulation should run.
     * @return number of seconds, the simulation should run
     */
    float endSimulation() {
      return duration;
    };

    /**
     * Get the boundary types of the simulation
     * @param edge specific edge
     * @return type of the edge
     */
    BoundaryType getBoundaryType( BoundaryEdge i_edge ) {
      //nothing other than outflow/transparent boundary conditions makes sense in a real simulation
      return OUTFLOW;
    }

    /** Get the boundary positions
     *
     * @param i_edge which edge
     * @return value in the corresponding dimension
     */
    float getBoundaryPos(BoundaryEdge i_edge) {
       if ( i_edge == BND_LEFT )
         return simulationArea[0];
       else if ( i_edge == BND_RIGHT)
         return simulationArea[1];
       else if ( i_edge == BND_BOTTOM )
         return simulationArea[2];
       else
         return simulationArea[3];
    };
};

#endif /* ASAGI_HPP_ */
