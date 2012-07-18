#ifndef __SWE_VISINFO_H
#define __SWE_VISINFO_H
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

/**
 * SWE_VisInfo defines an interface that can be used for online 
 * visualisation of a shallow water simulation. 
 * In particular, it provides information required for proper 
 * scaling of the involved variables. 
 */
class SWE_VisInfo {

 public:

    // Scaling factors for custom visualization
//     virtual bool useCustomScaling() { return false; };
    virtual float waterHeightAtRest() { return 10.0f; };
    virtual float waterDistanceFromGround() { return 10.0f; };
    virtual float waterVerticalScaling() { return 10.0f; };
    virtual float bathyVerticalCenter() { return 0.0f; };
    virtual float bathyDistanceFromGround() { return 0.0f; };
    virtual float bathyVerticalScaling() { return 10.0f; };

};

#endif
