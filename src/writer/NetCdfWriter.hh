/**
 * @file
 * This file is part of SWE.
 *
 * @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
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

#include <cstring>
#include <string>
#include <vector>
#include <netcdf.h>

#include "writer/Writer.hh"

namespace io {
  class NetCdfWriter;
}

class io::NetCdfWriter : public io::Writer {
private:
    /** netCDF file id*/
    int dataFile;

    /** Variable ids */
    int timeVar, hVar, huVar, hvVar, bVar;

    /** Flush after every x write operation? */
    unsigned int flush;

    // writer time dependent variables.
    void writeVarTimeDependent( const Float2D &i_matrix,
                                int i_ncVariable);

    // writes time independent variables.
    void writeVarTimeIndependent( const Float2D &i_matrix,
                                  int i_ncVariable);


  public:
    NetCdfWriter(const std::string &i_fileName,
    			 const Float2D &i_b,
                 const BoundarySize &i_boundarySize,
                 int i_nX, int i_nY,
                 float i_dX, float i_dY,
                 float i_originX = 0., float i_originY = 0.,
                 unsigned int i_flush = 0);
    virtual ~NetCdfWriter();

    // writes the unknowns at a given time step to the netCDF-file.
    void writeTimeStep( const Float2D &i_h,
                        const Float2D &i_hu,
                        const Float2D &i_hv,
                        float i_time);

  private:
    /**
     * This is a small wrapper for `nc_put_att_text` which automatically sets the length.
     */
    void ncPutAttText(int varid, const char* name, const char *value)
    {
    	nc_put_att_text(dataFile, varid, name, strlen(value), value);
    }

};

#endif /* NETCDFWRITER_HH_ */
