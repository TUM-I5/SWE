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

#include "Writer.hpp"
#if defined(ENABLE_WRITERS)
#if defined(ENABLE_NETCDF)
#include "NetCDFWriter.hpp"
#else
#include "VTKWriter.hpp"
#endif
#endif

#include <sstream>

std::string Writers::generateBaseFileName(const std::string& baseName, int blockPositionX, int blockPositionY) {
  std::ostringstream fileName;

  fileName << baseName << "_" << blockPositionX << blockPositionY;
  return fileName.str();
}

int& Writers::BoundarySize::operator[](unsigned int i) { return boundarySize[i]; }

int Writers::BoundarySize::operator[](unsigned int i) const { return boundarySize[i]; }

#if defined(ENABLE_WRITERS)
std::shared_ptr<Writers::Writer> Writers::Writer::createWriterInstance(
  [[maybe_unused]] const std::string&              fileName,
  [[maybe_unused]] const Tools::Float2D<RealType>& bathymetry,
  [[maybe_unused]] const BoundarySize&             boundarySize,
  [[maybe_unused]] int                             nX,
  [[maybe_unused]] int                             nY,
  [[maybe_unused]] RealType                        dX,
  [[maybe_unused]] RealType                        dY,
  [[maybe_unused]] RealType                        offsetX,
  [[maybe_unused]] RealType                        offsetY,
  [[maybe_unused]] RealType                        originX,
  [[maybe_unused]] RealType                        originY,
  [[maybe_unused]] int                             flush
) {

#if defined(ENABLE_NETCDF)
  auto writer = std::make_shared<NetCDFWriter>(
    fileName, bathymetry, boundarySize, nX, nY, dX, dY, originX, originY, flush
  );
#else
  auto writer = std::make_shared<VTKWriter>(fileName, bathymetry, boundarySize, nX, nY, dX, dY, offsetX, offsetY);
#endif

  return writer;
}
#endif

Writers::Writer::Writer(
  const std::string&              fileName,
  const Tools::Float2D<RealType>& bathymetry,
  const BoundarySize&             boundarySize,
  int                             nX,
  int                             nY
):
  fileName_(fileName),
  nX_(nX),
  nY_(nY),
  bathymetry_(bathymetry),
  boundarySize_(boundarySize),
  timeStep_(0) {}
