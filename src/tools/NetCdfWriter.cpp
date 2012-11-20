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

#include "NetCdfWriter.hh"
#include <string>
#include <vector>
#include <iostream>
#include <cassert>

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
	nc_close(dataFile);
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
	int status;

	//create a netCDF-file, an existing file will be replaced
	status = nc_create(fileName.c_str(), NC_NETCDF4, &dataFile);

  //check if the netCDF-file creation constructor succeeded.
	if (status != NC_NOERR) {
		assert(false);
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
	int l_timeDim, l_xDim, l_yDim;
	nc_def_dim(dataFile, "time", NC_UNLIMITED, &l_timeDim);
	nc_def_dim(dataFile, "x", nX, &l_xDim);
	nc_def_dim(dataFile, "y", nY, &l_yDim);

	//variables (TODO: add rest of CF-1.5)
	int l_xVar, l_yVar;

	nc_def_var(dataFile, "time", NC_FLOAT, 1, &l_timeDim, &timeVar);
	ncPutAttText(timeVar, "long_name", "Time");
	ncPutAttText(timeVar, "units", "seconds since simulation start"); // the word "since" is important for the paraview reader

	nc_def_var(dataFile, "x", NC_FLOAT, 1, &l_xDim, &l_xVar);
	nc_def_var(dataFile, "y", NC_FLOAT, 1, &l_yDim, &l_yVar);

	//variables, fastest changing index is on the right (C syntax), will be mirrored by the library
	int dims[] = {l_timeDim, l_yDim, l_xDim};
	nc_def_var(dataFile, "h",  NC_FLOAT, 3, dims, &hVar);
	nc_def_var(dataFile, "hu", NC_FLOAT, 3, dims, &huVar);
	nc_def_var(dataFile, "hv", NC_FLOAT, 3, dims, &hvVar);
	nc_def_var(dataFile, "b",  NC_FLOAT, 2, &dims[1], &bVar);

	//set attributes to match CF-1.5 convention
	ncPutAttText(NC_GLOBAL, "Conventions", "CF-1.5");
	ncPutAttText(NC_GLOBAL, "title", "Computed tsunami solution");
	ncPutAttText(NC_GLOBAL, "history", "SWE");
	ncPutAttText(NC_GLOBAL, "institution", "Technische Universitaet Muenchen, Department of Informatics, Chair of Scientific Computing");
	ncPutAttText(NC_GLOBAL, "source", "Bathymetry and displacement data.");
	ncPutAttText(NC_GLOBAL, "references", "http://www5.in.tum.de/SWE");
	ncPutAttText(NC_GLOBAL, "comment", "SWE is free software and licensed under the GNU General Public License. Remark: In general this does not hold for the used input data.");

	//setup grid size
	float gridPosition = i_originX + (float).5 * i_dX;
	for(size_t i = 0; i < nX; i++) {
		nc_put_var1_float(dataFile, l_xVar, &i, &gridPosition);

		gridPosition += i_dX;
	}

	gridPosition = i_originY + (float).5 * i_dY;
	for(size_t j = 0; j < nY; j++) {
		nc_put_var1_float(dataFile, l_yVar, &j, &gridPosition);

    	gridPosition += i_dY;
	}
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
 * @param i_ncVariable time dependent netCDF-variable to which the output is written to.
 */
void io::NetCdfWriter::writeVarTimeDependent( const Float2D &i_matrix,
                                              const int i_boundarySize[4],
                                              int i_ncVariable ) {
	//write col wise, necessary to get rid of the boundary
	//storage in Float2D is col wise
	//read carefully, the dimensions are confusing
	size_t start[] = {timeStep, 0, 0};
	size_t count[] = {1, nY, 1};
	for(int col = 0; col < nX; col++) {
		start[2] = col; //select col (dim "x")
		nc_put_vara_float(dataFile, i_ncVariable, start, count,
				&i_matrix[col+i_boundarySize[0]][i_boundarySize[2]]); //write col
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
 * @param i_ncVariable time independent netCDF-variable to which the output is written to.
 */
void io::NetCdfWriter::writeVarTimeIndependent( const Float2D &i_matrix,
                                                const int i_boundarySize[4],
                                                int i_ncVariable ) {
	//write col wise, necessary to get rid of the boundary
	//storage in Float2D is col wise
	//read carefully, the dimensions are confusing
	size_t start[] = {0, 0};
	size_t count[] = {nY, 1};
	for(int col = 0; col < nX; col++) {
		start[1] = col; //select col (dim "x")
		nc_put_vara_float(dataFile, i_ncVariable, start, count,
				&i_matrix[col+i_boundarySize[0]][i_boundarySize[2]]); //write col
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
	//write i_time
	nc_put_var1_float(dataFile, timeVar, &timeStep, &i_time);

	//write water height
	writeVarTimeDependent(i_h, i_boundarySize, hVar);

	//write momentum in x-direction
	writeVarTimeDependent(i_hu, i_boundarySize, huVar);

	//write momentum in y-direction
	writeVarTimeDependent(i_hv, i_boundarySize, hvVar);

	// Increment timeStep for next call
	timeStep++;
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
