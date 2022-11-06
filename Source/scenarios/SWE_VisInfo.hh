/**
 * @file
 * This file is part of SWE.
 *
 * @author Michael Bader
 * @author Kaveh Rahnema
 * @author Tobias Schnabel
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

#ifndef __SWE_VISINFO_H
#define __SWE_VISINFO_H

#include "SWE_Scenario.hh"

/**
 * SWE_VisInfo defines an interface that can be used for online 
 * visualization of a shallow water simulation.
 * In particular, it provides information required for proper 
 * scaling of the involved variables.
 *
 * For water height:
 * displayedWaterHeight = waterVerticalScaling() * simulatedWaterHeight
 *
 * For bathymetry:
 * displayedBatyhmetry = bathyVerticalScaling() * realBathymetry
 *  + bathyVerticalOffset()
 *
 * The default water height should be 0. In this case a bathymetry value
 * smaller than 0 means water and a value greater than 0 is land. Therefore
 * bathyVerticalOffset should 0 for all real scenarios.
 *
 * If you do not not provide an SWE_VisInfo for scenario,
 * (water|bathy)VerticalScaling will be guessed form the value initial
 * values. bathyVerticalOffset is always 0 in this case.
 */
class SWE_VisInfo {

 public:
	/**
	 * Empty virtual destructor
	 */
	virtual ~SWE_VisInfo() {};

    /**
     * @return The vertical scaling factor of the water
     */
    virtual float waterVerticalScaling() { return 10.0f; };

    /**
     * @return The vertical offset for the bathymetry. Should be
     *  0 for "real" scenarios (scenarios with dry areas)
     */
    virtual float bathyVerticalOffset() { return 0.0f; };

    /**
     * @return The vertical scaling factor for the bathymetry
     */
    virtual float bathyVerticalScaling() { return 10.0f; };

};

#endif
