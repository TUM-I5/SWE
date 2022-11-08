/**
 * @file
 * This file is part of SWE.
 *
 * @author Michael Bader, Kaveh Rahnema, Tobias Schnabel
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de,
 * http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
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

#include "Scenario.hpp"

namespace Scenarios {

  /**
   * Scenario "Splashing Cone":
   * bathymetry forms a circular cone
   * intial water surface designed to form "sea at rest"
   * but: elevated water region in the centre (similar to radial dam break)
   */
  class SplashingConeScenario: public Scenario {
  public:
    ~SplashingConeScenario() override = default;

    RealType getWaterHeightAtRest() const;
    RealType getWaterHeight(RealType x, RealType y) const override;
    RealType getBathymetry(RealType x, RealType y) const override;

    double getEndSimulationTime() const override;

    BoundaryType getBoundaryType(BoundaryEdge edge) const override;
  };

} // namespace Scenarios
