#ifndef __SWE_VISINFO_H
#define __SWE_VISINFO_H
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

#include "SWE_Scenario.h"

/**
 * SWE_VisInfo defines an interface that can be used for online 
 * visualisation of a shallow water simulation. 
 * In particular, it provides information required for proper 
 * scaling of the involved variables. 
 */
class SWE_VisInfo {

 public:
	/**
	 * Empty virtual destructor
	 */
	virtual ~SWE_VisInfo() {};

    // Scaling factors for custom visualization
    virtual float waterVerticalScaling() { return 10.0f; };
    virtual float bathyVerticalOffset() { return 0.0f; };
    virtual float bathyVerticalScaling() { return 10.0f; };

};

#endif
