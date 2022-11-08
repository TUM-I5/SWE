/**
 * @file
 * This file is part of SWE.
 *
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
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "Tools/Float2D.hpp"
#include "Tools/RealType.hpp"

namespace Writers {

  std::string generateBaseFileName(const std::string& baseName, int blockPositionX, int blockPositionY);

  /**
   * This struct is used so we can initialize this array
   * in the constructor.
   */
  struct BoundarySize {
    /**
     * boundarySize[0] == Left
     * boundarySize[1] == Right
     * boundarySize[2] == Bottom
     * boundarySize[3] == Top
     */
    int boundarySize[4];

    int& operator[](unsigned int i);
    int  operator[](unsigned int i) const;
  };

  class Writer {
  protected:
    const std::string fileName_;

    const int nX_;
    const int nY_;

    const Tools::Float2D<RealType>& bathymetry_;
    const BoundarySize              boundarySize_;

    int timeStep_;

  public:
    static std::shared_ptr<Writer> createWriterInstance(
      const std::string&              fileName,
      const Tools::Float2D<RealType>& bathymetry,
      const BoundarySize&             boundarySize,
      int                             nX,
      int                             nY,
      RealType                        dX,
      RealType                        dY,
      RealType                        offsetX,
      RealType                        offsetY,
      RealType                        originX,
      RealType                        originY,
      int                             flush
    );

    Writer(
      const std::string&              fileName,
      const Tools::Float2D<RealType>& bathymetry,
      const BoundarySize&             boundarySize,
      int                             nX,
      int                             nY
    );

    virtual ~Writer() = default;

    /**
     * Writes one time step
     *
     * @param h water heights at a given time step.
     * @param hu momentums in x-direction at a given time step.
     * @param hv momentums in y-direction at a given time step.
     * @param time simulation time of the time step.
     */
    virtual void writeTimeStep(
      const Tools::Float2D<RealType>& h,
      const Tools::Float2D<RealType>& hu,
      const Tools::Float2D<RealType>& hv,
      double                          time
    ) = 0;
  };

} // namespace Writers