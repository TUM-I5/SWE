#ifndef __SWE_SIMPLE_SCENARIOS_VIS_H
#define __SWE_SIMPLE_SCENARIOS_VIS_H
// =====================================================================
// This file is part of SWE_CUDA (see file SWE_Block.cu for details).
// 
// Copyright (C) 2010,2011 Michael Bader, Kaveh Rahnema, Tobias Schnabel
// Copyright (C) 2012      Sebastian Rettenberger
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

#include "SWE_VisInfo.hh"

/**
 * VisInfo "Bathymetry Dam Break":
 * uniform water depth, but elevated bathymetry in the center of the domain
 * Set bathymetry offset hence it is visible in the screen
 */
class SWE_BathymetryDamBreakVisInfo : public SWE_VisInfo {
public:
	float bathyVerticalOffset() { return 245.0f; };
};

#endif
