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

#pragma once

#include "Writer.hpp"

namespace Writers {

#if defined(ENABLE_WRITERS) && defined(ENABLE_NETCDF)

  class NetCDFWriter: public Writer {
  private:
    int dataFile_;

    int timeVar_, hVar_, huVar_, hvVar_, bVar_;

    /** Flush after every x write operation? */
    unsigned int flush_;

    /**
     * Writes time dependent data to a netCDF-file (-> constructor) with respect to the boundary sizes.
     *
     * boundarySize[0] == left
     * boundarySize[1] == right
     * boundarySize[2] == bottom
     * boundarySize[3] == top
     *
     * @param matrix array which contains time dependent data.
     * @param boundarySize size of the boundaries.
     * @param ncVariable time dependent netCDF-variable to which the output is written to.
     */
    void writeVarTimeDependent(const Tools::Float2D<RealType>& matrix, int ncVariable);

    /**
     * Write time independent data to a netCDF-file (-> constructor) with respect to the boundary sizes.
     * Variable is time-independent
     * boundarySize[0] == left
     * boundarySize[1] == right
     * boundarySize[2] == bottom
     * boundarySize[3] == top
     *
     * @param matrix array which contains time independent data.
     * @param boundarySize size of the boundaries.
     * @param ncVariable time independent netCDF-variable to which the output is written to.
     */
    void writeVarTimeIndependent(const Tools::Float2D<RealType>& matrix, int ncVariable);

    /**
     * This is a small wrapper for `nc_put_att_text` which automatically sets the length.
     */
    void ncPutAttText(int varid, const char* name, const char* value);

  public:
    NetCDFWriter(
      const std::string&              fileName,
      const Tools::Float2D<RealType>& bathymetry,
      const BoundarySize&             boundarySize,
      int                             nX,
      int                             nY,
      RealType                        dX,
      RealType                        dY,
      RealType                        originX = 0.,
      RealType                        originY = 0.,
      unsigned int                    flush   = 0
    );

    ~NetCDFWriter() override;

    /**
     * Writes the unknwons to a netCDF-file (-> constructor) with respect to the boundary sizes.
     *
     * boundarySize[0] == left
     * boundarySize[1] == right
     * boundarySize[2] == bottom
     * boundarySize[3] == top
     *
     * @param h water heights at a given time step.
     * @param hu momentums in x-direction at a given time step.
     * @param hv momentums in y-direction at a given time step.
     * @param boundarySize size of the boundaries.
     * @param time simulation time of the time step.
     */
    void writeTimeStep(
      const Tools::Float2D<RealType>& h,
      const Tools::Float2D<RealType>& hu,
      const Tools::Float2D<RealType>& hv,
      double                          time
    ) override;
  };

#endif

} // namespace Writers
