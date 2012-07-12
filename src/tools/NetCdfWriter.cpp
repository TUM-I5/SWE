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

#include "NetCdfWriter.hh"
#include <string>
#include <vector>
#include <iostream>
#include <cassert>
#include <netcdfcpp.h>

/**
 * Constructor of a netCDF-writer.
 *
 * @param i_filename file of the netCDF-file to which the data will be written to.
 * @param i_nX number of cells in the horizontal direction.
 * @param i_nY number of cells in the vertical direction.
 */
io::NetCdfWriter::NetCdfWriter( const std::string &i_filename,
                                const int &i_nX,
                                const int &i_nY ):
  nX(i_nX), nY(i_nY), timeStep(0), fileName(i_filename) {
}

/**
 * Destructor of a netCDF-writer.
 */
io::NetCdfWriter::~NetCdfWriter() {
}

/**
 * Create a netCdf-file (-> constructor)
 * Any existing file will be replaced.
 *
 * @param i_dX cell size in x-direction.
 * @param i_dY cell size in y-direction.
 * @param i_originX
 * @param i_originY
 * @param i_dynamicBathymetry
 */
void io::NetCdfWriter::createNetCdfFile( const float &i_dX,
                                         const float &i_dY,
                                         const float &i_originX,
                                         const float &i_originY ) {//,
                                         //const bool  &i_dynamicBathymetry) { //!TODO
  //create a netCDF-file, an existing file will be replaced
  NcFile l_dataFile = NcFile(fileName.c_str(), NcFile::Replace);

  //check if the netCDF-file creation constructor succeeded.
  if (!l_dataFile.is_valid()) {
    return;
  }

#ifdef PRINT_NETCDFWRITER_INFORMATION
  std::cout << "   *** io::NetCdfWriter::createNetCdfFile" << std::endl;
  std::cout << "     created/replaced: " << fileName << std::endl;
  std::cout << "     dimensions(nx, ny): " << nX << ", " << nY << std::endl;
  std::cout << "     cell width(dx,dy): " << i_dX << ", " << i_dY << std::endl;
  std::cout << "     origin(x,y): " << i_originX << ", " << i_originY << std::endl;
#endif

  //dimensions
  NcDim* l_timeDim = l_dataFile.add_dim("time", NC_UNLIMITED);
  NcDim* l_xDim = l_dataFile.add_dim("x", nX);
  NcDim* l_yDim = l_dataFile.add_dim("y", nY);

  //variables (TODO: add rest of CF-1.5)
  NcVar* l_ncVariable;

  l_ncVariable = l_dataFile.add_var("time", ncFloat, l_timeDim);
  l_ncVariable->add_att("long_name", "Time");
  l_ncVariable->add_att("units", "seconds since simulation start"); // the word "since" is important for the paraview reader

  l_dataFile.add_var("x", ncFloat, l_xDim);
  l_dataFile.add_var("y", ncFloat, l_yDim);

  //variables, fastest changing index is on the right (C syntax), will be mirrored by the library
  l_dataFile.add_var("h",  ncFloat, l_timeDim, l_xDim, l_yDim);
  l_dataFile.add_var("hu", ncFloat, l_timeDim, l_xDim, l_yDim);
  l_dataFile.add_var("hv", ncFloat, l_timeDim, l_xDim, l_yDim);
  l_dataFile.add_var("b",  ncFloat, l_xDim, l_yDim);

  //set attributes to match CF-1.5 convention
  l_dataFile.add_att("Conventions", "CF-1.5");
  l_dataFile.add_att("title", "Computed tsunami solution");
  l_dataFile.add_att("history", "SWE");
  l_dataFile.add_att("institution", "Technische Universitaet Muenchen, Department of Informatics, Chair of Scientific Computing");
  l_dataFile.add_att("source", "Bathymetry and displacement data.");
  l_dataFile.add_att("references", "http://www5.in.tum.de/SWE");
  l_dataFile.add_att("comment", "SWE is free software and licensed under the GNU General Public License. Remark: In general this does not hold for the used input data.");


  //setup grid size
  float gridPosition = i_originX + (float).5 * i_dX;
  for(int i = 0; i < nX; i++) {
    l_ncVariable = l_dataFile.get_var("x");
    l_ncVariable->set_cur(i);
    l_ncVariable->put(&gridPosition, 1);

    gridPosition += i_dX;
  }

  gridPosition = i_originY + (float).5 * i_dY;
  for(int j = 0; j < nY; j++) {
    l_ncVariable = l_dataFile.get_var("y");
    l_ncVariable->set_cur(j);
    l_ncVariable->put(&gridPosition, 1);

    gridPosition += i_dY;
  }

  //l_dataFile closed automatically (out of scope)
}

/**
 * Writes time dependent data to a netCDF-file (-> constructor) with respect to the boundary sizes.
 *
 * boundarySize[0] == left
 * boundarySize[1] == right
 * boundarySize[2] == bottom
 * boundarySize[3] == top
 *
 * @param i_matrix array which contains time dependent data.
 * @param i_boundarySize size of the boundaries.
 * @param o_ncVariable time dependent netCDF-variable to which the output is written to.
 */
void io::NetCdfWriter::writeVarTimeDependent( const Float2D &i_matrix,
                                              const int i_boundarySize[4],
                                              NcVar* o_ncVariable ) {
  //write col wise, necessary to get rid of the boundary
  //storage in Float2D is col wise
  //read carefully, the dimensions are confusing
  for(int col = 0; col < nX; col++) {
    o_ncVariable->set_cur(timeStep, col, 0); //select col (dim "x")
    o_ncVariable->put(&i_matrix[col+i_boundarySize[0]][i_boundarySize[2]], 1, 1, nY); //write col
  }
}

/**
 * Write time independent data to a netCDF-file (-> constructor) with respect to the boundary sizes.
 * Variable is time-independent
 * boundarySize[0] == left
 * boundarySize[1] == right
 * boundarySize[2] == bottom
 * boundarySize[3] == top
 *
 * @param i_matrix array which contains time independent data.
 * @param i_boundarySize size of the boundaries.
 * @param o_ncVariable time independent netCDF-variable to which the output is written to.
 */
void io::NetCdfWriter::writeVarTimeIndependent( const Float2D &i_matrix,
                                                const int i_boundarySize[4],
                                                NcVar* o_ncVariable ) {
  //write col wise, necessary to get rid of the boundary
  //storage in Float2D is col wise
  //read carefully, the dimensions are confusing
  for(int col = 0; col < nX; col++) {
    o_ncVariable->set_cur(col, 0); //select col (dim "x")
    o_ncVariable->put(&i_matrix[col+i_boundarySize[0]][i_boundarySize[2]], 1, nY); //write col
  }
}

/**
 * Writes the unknwons to a netCDF-file (-> constructor) with respect to the boundary sizes.
 *
 * boundarySize[0] == left
 * boundarySize[1] == right
 * boundarySize[2] == bottom
 * boundarySize[3] == top
 *
 * @param i_h water heights at a given time step.
 * @param i_hu momentums in x-direction at a given time step.
 * @param i_hv momentums in y-direction at a given time step.
 * @param i_boundarySize size of the boundaries.
 * @param i_time simulation time of the time step.
 */
void io::NetCdfWriter::writeUnknowns( const Float2D &i_h,
                                      const Float2D &i_hu,
                                      const Float2D &i_hv,
                                      const int i_boundarySize[4],
                                      const float &i_time) {
  //open dataFile
  NcFile l_dataFile = NcFile(fileName.c_str(), NcFile::Write);

  //check if the netCDF file open constructor succeeded.
  if (!l_dataFile.is_valid()) {
    std::cerr << "could not open " << fileName << std::endl;
    assert(false);
  }
  else {
    //set time
    timeStep = l_dataFile.get_dim("time")->size();
  }

  //write i_time
  NcVar* l_ncVariable = l_dataFile.get_var("time");
  l_ncVariable->set_cur(timeStep);
  l_ncVariable->put(&i_time, 1);

  //write water height
  l_ncVariable = l_dataFile.get_var("h");
  writeVarTimeDependent(i_h, i_boundarySize, l_ncVariable);

  //write momentum in x-direction
  l_ncVariable = l_dataFile.get_var("hu");
  writeVarTimeDependent(i_hu, i_boundarySize, l_ncVariable);

  //write momentum in y-direction
  l_ncVariable = l_dataFile.get_var("hv");
  writeVarTimeDependent(i_hv, i_boundarySize, l_ncVariable);
}


/**
 * Write (static) bathymetry data to a netCDF-file (-> constructor) with respect to the boundary sizes.
 *
 * @param i_b bathymetry data.
 * @param i_boundarySize size of the boundaries.
 */
//dynamic possible -> writeUnknowns
void io::NetCdfWriter::writeBathymetry( const Float2D &i_b,
                                        const int i_boundarySize[4] ) {
  //open dataFile
  NcFile l_dataFile = NcFile(fileName.c_str(), NcFile::Write);

  //check if the netCDF file open constructor succeeded.
  if (!l_dataFile.is_valid()) {
    std::cerr << "could not open " << fileName << std::endl;
    assert(false);
  }
  NcVar* l_ncVariable = l_dataFile.get_var("b");
  writeVarTimeIndependent(i_b, i_boundarySize, l_ncVariable);
}

/**
 * Writes one dimensional unknowns to a netCDF-file (-> constructor) with respect to the boundary sizes.
 *
 * @param i_h water heights.
 * @param i_hu momentums in x-direction.
 * @param i_leftBoundarySize size of the left boundary.
 * @param i_rightBoundarySize size of the right boundary.
 * @param i_time simulation time of the time step.
 */
void io::NetCdfWriter::writeUnknownsOneDimensional( const std::vector<double> &i_h,
                                                    const std::vector<double> &i_hu,
                                                    const int &i_leftBoundarySize,
                                                    const int &i_rightBoundarySize,
                                                    const double &i_time ) {
  Float2D l_hTwoDimensional( i_h.size(), 1);
  Float2D l_huTwoDimensional( i_hu.size(), 1);
  Float2D l_hvTwoDimensional( i_hu.size(), 1);

  for(unsigned int i = 0; i < i_h.size(); i++) {
    l_hTwoDimensional[i][0] = i_h[i];
    l_huTwoDimensional[i][0] = i_hu[i];
    l_hvTwoDimensional[i][0] = 0.;
  }

  int l_boundarySizeTwoDimensional[4];
  l_boundarySizeTwoDimensional[0] = i_leftBoundarySize;
  l_boundarySizeTwoDimensional[1] = i_rightBoundarySize;
  l_boundarySizeTwoDimensional[2] = l_boundarySizeTwoDimensional[3] = 0;

  writeUnknowns(l_hTwoDimensional, l_huTwoDimensional, l_hvTwoDimensional, l_boundarySizeTwoDimensional, i_time);
}

/**
 * Writes one dimensional (static) bathymetry data to a netCDF-File (-> constructor) with respect to the boundary sizes.
 *
 * @param i_b bathymetry values.
 * @param i_leftBoundarySize size of the left boundary.
 * @param i_rightBoundarySize size of the right boundary.
 */
void io::NetCdfWriter::writeBathymetryOneDimensional( const std::vector<double> &i_b,
                                                      const int &i_leftBoundarySize,
                                                      const int &i_rightBoundarySize) {
  Float2D l_bTwoDimensional( i_b.size(), 1);

  for(unsigned int i = 0; i < i_b.size(); i++) {
    l_bTwoDimensional[i][0] = i_b[i];
  }

  int l_boundarySizeTwoDimensional[4];
  l_boundarySizeTwoDimensional[0] = i_leftBoundarySize;
  l_boundarySizeTwoDimensional[1] = i_rightBoundarySize;
  l_boundarySizeTwoDimensional[2] = l_boundarySizeTwoDimensional[3] = 0;

  writeBathymetry(l_bTwoDimensional, l_boundarySizeTwoDimensional);
}
