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

#include "Writer.hpp"

namespace Writers {

  class VTKWriter: public Writer {
    RealType dX_, dY_;

    RealType offsetX_, offsetY_;

    std::string generateFileName() const;

  public:
    /**
     * Creates a vtk file for each time step.
     * Any existing file will be replaced.
     *
     * @param baseName base name of the netCDF-file to which the data will be written to.
     * @param nX number of cells in the horizontal direction.
     * @param nY number of cells in the vertical direction.
     * @param dX cell size in x-direction.
     * @param dY cell size in y-direction.
     * @param offsetX x-offset of the block
     * @param offsetY y-offset of the block
     * @param bathymetry
     *
     * @TODO: This version can only handle a boundary layer of size 1
     */
    VTKWriter(
      const std::string&              fileName,
      const Tools::Float2D<RealType>& bathymetry,
      const BoundarySize&             boundarySize,
      int                             nX,
      int                             nY,
      RealType                        dX,
      RealType                        dY,
      RealType                        offsetX = 0,
      RealType                        offsetY = 0
    );

    ~VTKWriter() override = default;

    void writeTimeStep(
      const Tools::Float2D<RealType>& h,
      const Tools::Float2D<RealType>& hu,
      const Tools::Float2D<RealType>& hv,
      double                          time
    ) override;
  };

} // namespace Writers