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

#pragma once

#include "Tools/RealType.hpp"
#include "Types/BoundaryEdge.hpp"
#include "Types/BoundaryType.hpp"

namespace Scenarios {

  /**
   * Scenarios::Scenario defines an interface to initialise the unknowns of a
   * shallow water simulation - i.e. to initialise water height, velocities,
   * and bathymatry according to certain scenarios.
   * Scenarios::Scenario can act as stand-alone scenario class, providing a very
   * basic scenario (all functions are constant); however, the idea is
   * to provide derived classes that implement the Scenarios::Scenario interface
   * for more interesting scenarios.
   */
  class Scenario {
  public:
    virtual ~Scenario() = default;

    virtual RealType getWaterHeight(RealType x, RealType y) const;
    virtual RealType getVelocityU(RealType x, RealType y) const;
    virtual RealType getVelocityV(RealType x, RealType y) const;
    virtual RealType getBathymetry(RealType x, RealType y) const;

    virtual RealType getWaterHeightAtRest() const;
    virtual double   getEndSimulationTime() const;

    virtual BoundaryType getBoundaryType(BoundaryEdge edge) const;
    virtual RealType     getBoundaryPos(BoundaryEdge edge) const;
  };

} // namespace Scenarios
