/**
 * @file
 * This file is part of SWE.
 *
 * @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
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
 * A writer for the netCDF-format: http://www.unidata.ucar.edu/software/netcdf/
 */

#include "NetCDFWriter.hpp"

#if defined(ENABLE_WRITERS) && defined(ENABLE_NETCDF)

#include <cassert>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>
#ifdef ENABLE_MPI
#include <mpi.h>
#ifndef MPI_INCLUDED
#define MPI_INCLUDED
#define MPI_INCLUDED_NETCDF
#endif
#endif
#include <netcdf.h>
#ifdef MPI_INCLUDED_NETCDF
#undef MPI_INCLUDED
#undef MPI_INCLUDED_NETCDF
#endif

void Writers::NetCDFWriter::ncPutAttText(int varid, const char* name, const char* value) {
  nc_put_att_text(dataFile_, varid, name, strlen(value), value);
}

Writers::NetCDFWriter::NetCDFWriter(
  const std::string&              baseName,
  const Tools::Float2D<RealType>& bathymetry,
  const BoundarySize&             boundarySize,
  int                             nX,
  int                             nY,
  RealType                        dX,
  RealType                        dY,
  RealType                        originX,
  RealType                        originY,
  unsigned int                    flush
):
  Writer(baseName + ".nc", bathymetry, boundarySize, nX, nY),
  flush_(flush) {
  int status = -1;

  // Create a netCDF-file, an existing file will be replaced
  status = nc_create(fileName_.c_str(), NC_NETCDF4, &dataFile_);

  // Check if the netCDF-file creation constructor succeeded.
  if (status != NC_NOERR) {
    assert(false);
    return;
  }

#ifndef NDEBUG
  std::cout << "   *** Writers::NetCDFWriter::NetCDFWriter" << std::endl;
  std::cout << "     created/replaced: " << fileName_ << std::endl;
  std::cout << "     dimensions(nx, ny): " << nX << ", " << nY << std::endl;
  std::cout << "     cell width(dx, dy): " << dX << ", " << dY << std::endl;
  std::cout << "     origin(x, y): " << originX << ", " << originY << std::endl;
#endif

  // Dimensions
  int l_timeDim, l_xDim, l_yDim;
  nc_def_dim(dataFile_, "time", NC_UNLIMITED, &l_timeDim);
  nc_def_dim(dataFile_, "x", nX, &l_xDim);
  nc_def_dim(dataFile_, "y", nY, &l_yDim);

  // Variables (TODO: add rest of CF-1.5)
  int l_xVar, l_yVar;

  nc_def_var(dataFile_, "time", NC_FLOAT, 1, &l_timeDim, &timeVar_);
  ncPutAttText(timeVar_, "long_name", "Time");
  ncPutAttText(
    timeVar_, "units", "seconds since simulation start"
  ); // The word "since" is important for the ParaView reader

  nc_def_var(dataFile_, "x", NC_FLOAT, 1, &l_xDim, &l_xVar);
  nc_def_var(dataFile_, "y", NC_FLOAT, 1, &l_yDim, &l_yVar);

  // Variables, fastest changing index is on the right (C syntax), will be mirrored by the library
  int dims[] = {l_timeDim, l_yDim, l_xDim};
  nc_def_var(dataFile_, "h", NC_FLOAT, 3, dims, &hVar_);
  nc_def_var(dataFile_, "hu", NC_FLOAT, 3, dims, &huVar_);
  nc_def_var(dataFile_, "hv", NC_FLOAT, 3, dims, &hvVar_);
  nc_def_var(dataFile_, "b", NC_FLOAT, 2, &dims[1], &bVar_);

  // Set attributes to match CF-1.5 convention
  ncPutAttText(NC_GLOBAL, "Conventions", "CF-1.5");
  ncPutAttText(NC_GLOBAL, "title", "Computed tsunami solution");
  ncPutAttText(NC_GLOBAL, "history", "SWE");
  ncPutAttText(
    NC_GLOBAL,
    "institution",
    "Technische Universitaet Muenchen, Department of Informatics, Chair of Scientific Computing"
  );
  ncPutAttText(NC_GLOBAL, "source", "Bathymetry and displacement data.");
  ncPutAttText(NC_GLOBAL, "references", "http://www5.in.tum.de/SWE");
  ncPutAttText(
    NC_GLOBAL,
    "comment",
    "SWE is free software and licensed under the GNU General Public License. Remark: In general this does not hold for "
    "the used input data."
  );

  // Setup grid size
  float gridPosition = originX + (float).5 * dX;
  for (size_t i = 0; i < nX; i++) {
    nc_put_var1_float(dataFile_, l_xVar, &i, &gridPosition);

    gridPosition += dX;
  }

  gridPosition = originY + (float).5 * dY;
  for (size_t j = 0; j < nY; j++) {
    nc_put_var1_float(dataFile_, l_yVar, &j, &gridPosition);

    gridPosition += dY;
  }
}

Writers::NetCDFWriter::~NetCDFWriter() { nc_close(dataFile_); }

void Writers::NetCDFWriter::writeVarTimeDependent(const Tools::Float2D<RealType>& matrix, int ncVariable) {
  // Write column wise, necessary to get rid of the boundary
  // Storage in Float2D is column wise
  // Read carefully, the dimensions are confusing
  std::size_t start[] = {static_cast<std::size_t>(timeStep_), 0, 0};
  std::size_t count[] = {1, static_cast<std::size_t>(nY_), 1};
  for (unsigned int col = 0; col < nX_; col++) {
    start[2] = col; // Select column (dim "x")
    nc_put_vara_double(
      dataFile_,
      ncVariable,
      start,
      count,
      &matrix[col + boundarySize_[0]][boundarySize_[2]]
    ); // Write column
  }
}

void Writers::NetCDFWriter::writeVarTimeIndependent(const Tools::Float2D<RealType>& matrix, int ncVariable) {
  // Write column wise, necessary to get rid of the boundary
  // Storage in Float2D is column wise
  // Read carefully, the dimensions are confusing
  std::size_t start[] = {0, 0};
  std::size_t count[] = {static_cast<std::size_t>(nY_), 1};
  for (unsigned int col = 0; col < nX_; col++) {
    start[1] = col; // Select column (dim "x")
    nc_put_vara_double(
      dataFile_,
      ncVariable,
      start,
      count,
      &matrix[col + boundarySize_[0]][boundarySize_[2]]
    ); // Write column
  }
}

void Writers::NetCDFWriter::writeTimeStep(
  const Tools::Float2D<RealType>& h, const Tools::Float2D<RealType>& hu, const Tools::Float2D<RealType>& hv, double time
) {
  if (timeStep_ == 0) {
    // Write bathymetry
    writeVarTimeIndependent(bathymetry_, bVar_);
  }

  // Write time
  std::size_t timeStep = timeStep_;
  nc_put_var1_double(dataFile_, timeVar_, &timeStep, &time);

  // Write water height
  writeVarTimeDependent(h, hVar_);

  // Write momentum in x-direction
  writeVarTimeDependent(hu, huVar_);

  // Write momentum in y-direction
  writeVarTimeDependent(hv, hvVar_);

  // Increment timeStep for next call
  timeStep_++;

  if (flush_ > 0 && timeStep_ % flush_ == 0) {
    nc_sync(dataFile_);
  }
}

#endif
