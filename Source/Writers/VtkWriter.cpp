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

#include "VTKWriter.hpp"

#include <cassert>
#include <fstream>
#include <sstream>

std::string Writers::VTKWriter::generateFileName() const {
  std::ostringstream name;

  name << fileName_ << '.' << timeStep_ << ".vts";
  return name.str();
}

Writers::VTKWriter::VTKWriter(
  const std::string&              baseName,
  const Tools::Float2D<RealType>& bathymetry,
  const BoundarySize&             boundarySize,
  int                             nX,
  int                             nY,
  RealType                        dX,
  RealType                        dY,
  RealType                        offsetX,
  RealType                        offsetY
):
  Writer(baseName, bathymetry, boundarySize, nX, nY),
  dX_(dX),
  dY_(dY),
  offsetX_(offsetX),
  offsetY_(offsetY) {}

void Writers::VTKWriter::writeTimeStep(
  const Tools::Float2D<RealType>& h, const Tools::Float2D<RealType>& hu, const Tools::Float2D<RealType>& hv, double time
) {
  std::ofstream vtkFile(generateFileName().c_str());
  assert(vtkFile.good());
  (void)time; // TODO: Write time as additional information

  // VTK Header
  vtkFile
    << "<?xml version=\"1.0\"?>" << std::endl
    << "<VTKFile type=\"StructuredGrid\">" << std::endl
    << "<StructuredGrid WholeExtent=\"" << offsetX_ << " " << offsetX_ + nX_ << " " << offsetY_ << " " << offsetY_ + nY_
    << " 0 0\">" << std::endl
    << "<Piece Extent=\"" << offsetX_ << " " << offsetX_ + nX_ << " " << offsetY_ << " " << offsetY_ + nY_ << " 0 0\">"
    << std::endl;

  vtkFile
    << "<Points>" << std::endl
    << "<DataArray NumberOfComponents=\"3\" type=\"Float32\" format=\"ASCII\">" << std::endl;

  // Grid points
  for (int j = 0; j < nY_ + 1; j++) {
    for (int i = 0; i < nX_ + 1; i++) {
      vtkFile << (offsetX_ + i) * dX_ << " " << (offsetY_ + j) * dY_ << " 0" << std::endl;
    }
  }

  vtkFile << "</DataArray>" << std::endl << "</Points>" << std::endl;
  vtkFile << "<CellData>" << std::endl;

  // Water surface height h
  vtkFile << "<DataArray Name=\"h\" type=\"Float32\" format=\"ASCII\">" << std::endl;
  for (int j = 1; j < nY_ + 1; j++) {
    for (int i = 1; i < nX_ + 1; i++) {
      vtkFile << h[i][j] << std::endl;
    }
  }
  vtkFile << "</DataArray>" << std::endl;

  // Momentums
  vtkFile << "<DataArray Name=\"hu\" type=\"Float32\" format=\"ASCII\">" << std::endl;
  for (int j = 1; j < nY_ + 1; j++) {
    for (int i = 1; i < nX_ + 1; i++) {
      vtkFile << hu[i][j] << std::endl;
    }
  }
  vtkFile << "</DataArray>" << std::endl;

  vtkFile << "<DataArray Name=\"hv\" type=\"Float32\" format=\"ASCII\">" << std::endl;
  for (int j = 1; j < nY_ + 1; j++) {
    for (int i = 1; i < nX_ + 1; i++) {
      vtkFile << hv[i][j] << std::endl;
    }
  }
  vtkFile << "</DataArray>" << std::endl;

  // Bathymetry
  vtkFile << "<DataArray Name=\"b\" type=\"Float32\" format=\"ASCII\">" << std::endl;
  for (int j = 1; j < nY_ + 1; j++) {
    for (int i = 1; i < nX_ + 1; i++) {
      vtkFile << bathymetry_[i][j] << std::endl;
    }
  }
  vtkFile << "</DataArray>" << std::endl;

  vtkFile << "</CellData>" << std::endl << "</Piece>" << std::endl;

  vtkFile << "</StructuredGrid>" << std::endl << "</VTKFile>" << std::endl;

  timeStep_++;
}
