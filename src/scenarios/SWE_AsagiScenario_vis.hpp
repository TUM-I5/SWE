/**
 * @file
 * This file is part of SWE.
 *
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
 * Rescale water height in small Japan scenario
 */

#ifndef SWEASAGISCENARIO_VIS_HPP_
#define SWEASAGISCENARIO_VIS_HPP_

#include "SWE_VisInfo.hh"

class SWE_AsagiJapanSmallVisInfo : public SWE_VisInfo
{
public:
    virtual float waterVerticalScaling() { return 4.0f; };
    virtual float bathyVerticalScaling() { return 0.010313f; };
};

#endif // SWEASAGISCENARIO_VIS_HPP_
