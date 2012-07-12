/**
 * @file
 * This file is part of SWE.
 *
 * @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
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

#ifndef NETCDFWRITER_HH_
#define NETCDFWRITER_HH_

#include <string>
#include <vector>
#include <netcdfcpp.h>
#include "help.hh"

namespace io {
  class NetCdfWriter;
}

class io::NetCdfWriter {
  //private:
    //! dimensions of the grid in x- and y-direction.
    const int nX, nY;

    //! current time step of the netCDF-file
    int timeStep;

    //! file name of the netCDF-file
    const std::string fileName;

    // writer time dependent variables.
    void writeVarTimeDependent( const Float2D &i_matrix,
                                const int i_boundarySize[4],
                                NcVar* o_ncVariable);

    // writes time independent variables.
    void writeVarTimeIndependent( const Float2D &i_matrix,
                                  const int i_boundarySize[4],
                                  NcVar* o_ncVariable);


  public:
    NetCdfWriter(const std::string &i_fileName,
                 const int &i_nX,
                 const int &i_nY);
    virtual ~NetCdfWriter();

    // creates a netCDF-file with the specified properties.
    void createNetCdfFile( const float &i_dX,
                           const float &i_dY,
                           const float &i_originX = 0.,
                           const float &i_originY = 0.);

    // writes the unknowns at a given time step to the netCDF-file.
    void writeUnknowns( const Float2D &i_h,
                        const Float2D &i_hu,
                        const Float2D &i_hv,
                        const int i_boundarySize[4],
                        const float &i_time);

    // writes the time independent bathymetry (static displacement) to the netCDF-file.
    void writeBathymetry( const Float2D &i_b,
                          const int i_boundarySize[4] );

    // writes one dimensional unknowns to the netCDF-file (not used in SWE).
    void writeUnknownsOneDimensional( const std::vector<double> &i_h,
                                      const std::vector<double> &i_hu,
                                      const int &i_leftBoundarySize,
                                      const int &i_rightBoundarySize,
                                      const double &i_time);

    // writes one dimensional bathymetry to the netCDF-file (not used in SWE).
    void writeBathymetryOneDimensional( const std::vector<double> &i_b,
                                        const int &i_leftBoundarySize,
                                        const int &i_rightBoundarySize );

};

#endif /* NETCDFWRITER_HH_ */
