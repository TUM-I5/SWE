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

#include "SplashingConeScenario.hpp"

#include <cmath>

RealType Scenarios::SplashingConeScenario::getWaterHeightAtRest() const { return utilities::smart_cast<RealType>(4.0); }

RealType Scenarios::SplashingConeScenario::getWaterHeight(RealType x, RealType y) const {
  RealType r = utilities::smart_cast<RealType>(sqrt((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)));
  RealType h = utilities::smart_cast<RealType>(4.0 - 4.5 * (r / 0.5));

  if (r < 0.1)
    h = h + 1.0;

  return (h > 0.0) ? h : 0.0;
}

RealType Scenarios::SplashingConeScenario::getBathymetry([[maybe_unused]] RealType x, [[maybe_unused]] RealType y)
  const {
  RealType r = utilities::smart_cast<RealType>(sqrt((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)));
  return 1.0 + 9.0 * ((r < 0.5) ? r : 0.5);
}

double Scenarios::SplashingConeScenario::getEndSimulationTime() const { return 0.5; }

BoundaryType Scenarios::SplashingConeScenario::getBoundaryType([[maybe_unused]] BoundaryEdge edge) const {
  return BoundaryType::Outflow;
}
