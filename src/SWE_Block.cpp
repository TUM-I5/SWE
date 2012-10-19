/**
 * @file
 * This file is part of SWE.
 *
 * @author Michael Bader, Kaveh Rahnema, Tobias Schnabel
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
 * TODO
 */

#include "SWE_Block.hh"
#include "tools/help.hh"
#include <cmath>
#include <iostream>
#include <cassert>
#include <limits>


// define static variables:
  // block size: numer of cells in x and y direction:
  int SWE_Block::nx = 10;
  int SWE_Block::ny = 10;
  // grid sizes dx and dy:
  float SWE_Block::dx = 0.1f;
  float SWE_Block::dy = 0.1f;

  // gravitational acceleration
  const float SWE_Block::g = 9.81f;

/**
 * Constructor: allocate variables for simulation
 *
 * unknowns h (water height), hu,hv (discharge in x- and y-direction), 
 * and b (bathymetry) are defined on grid indices [0,..,nx+1]*[0,..,ny+1]
 * -> computational domain is [1,..,nx]*[1,..,ny]
 * -> plus ghost cell layer
 *
 * The constructor is protected: no instances of SWE_Block can be 
 * generated.
 *
 */
SWE_Block::SWE_Block()
: h(nx+2,ny+2), hu(nx+2,ny+2), hv(nx+2,ny+2), b(nx+2,ny+2)
{
  // set WALL as default boundary condition
  for(int i=0; i<4;i++) {
     boundary[i] = PASSIVE;
     neighbour[i] = NULL;
  };
  
}

/**
 * Destructor: de-allocate all variables
 */
SWE_Block::~SWE_Block() {
}

//==================================================================
// methods for external read/write to main variables h, hu, hv, and b
// Note: temporary and non-local variables depending on the main 
// variables are synchronised before/after their update or read
//==================================================================

/**
 * Initializes the unknowns and bathymetry in all grid cells according to the given SWE_Scenario.
 *
 * In the case of multiple SWE_Blocks at this point, it is not clear how the boundary conditions
 * should be set. This is because an isolated SWE_Block doesn't have any in information about the grid.
 * Therefore the calling routine, which has the information about multiple blocks, has to take care about setting
 * the right boundary conditions.
 * 
 * @param i_scenario scenario, which is used during the setup.
 * @param i_multipleBlocks are the multiple SWE_blocks?
 */
void SWE_Block::initScenario( float _offsetX, float _offsetY,
							  SWE_Scenario &i_scenario,
                              const bool i_multipleBlocks ) {
	offsetX = _offsetX;
	offsetY = _offsetY;

  // initialize water height and discharge
  for(int i=1; i<=nx; i++)
    for(int j=1; j<=ny; j++) {
      float x = offsetX + (i-0.5f)*dx;
      float y = offsetY + (j-0.5f)*dy;
      h[i][j] =  i_scenario.getWaterHeight(x,y);
      hu[i][j] = i_scenario.getVeloc_u(x,y) * h[i][j];
      hv[i][j] = i_scenario.getVeloc_v(x,y) * h[i][j]; 
    };

  // initialize bathymetry
  for(int i=0; i<=nx+1; i++) {
    for(int j=0; j<=ny+1; j++) {
      b[i][j] = i_scenario.getBathymetry( offsetX + (i-0.5f)*dx,
                                          offsetY + (j-0.5f)*dy );
    }
  }

  // in the case of multiple blocks the calling routine takes care about proper boundary conditions.
  if( i_multipleBlocks == false ) {
    // obtain boundary conditions for all four edges from scenario
    setBoundaryType(BND_LEFT, i_scenario.getBoundaryType(BND_LEFT));
    setBoundaryType(BND_RIGHT, i_scenario.getBoundaryType(BND_RIGHT));
    setBoundaryType(BND_BOTTOM, i_scenario.getBoundaryType(BND_BOTTOM));
    setBoundaryType(BND_TOP, i_scenario.getBoundaryType(BND_TOP));
  }

  // perform update after external write to variables 
  synchAfterWrite();

}

/**
 * set water height h in all interior grid cells (i.e. except ghost layer) 
 * to a uniform value
 */
void SWE_Block::setWaterHeight(float _h) {

  for(int i=1; i<=nx; i++)
    for(int j=1; j<=ny; j++) {
      h[i][j] = _h;
    };

  synchWaterHeightAfterWrite();
}

/**
 * set water height h in all interior grid cells (i.e. except ghost layer) 
 * to values specified by parameter function _h
 */
void SWE_Block::setWaterHeight(float (*_h)(float, float)) {

  for(int i=1; i<=nx; i++)
    for(int j=1; j<=ny; j++) {
      h[i][j] =  _h(offsetX + (i-0.5f)*dx, offsetY + (j-0.5f)*dy);
    };

  synchWaterHeightAfterWrite();
}

/**
 * set discharge unknowns in all interior grid cells (i.e. except ghost layer) 
 * to a uniform value
 * Note: unknowns hu and hv represent momentum, while parameters u and v are velocities!
 *    (_u and _v are multiplied with the resp. cell-local value of h) 
 */
void SWE_Block::setDischarge(float _u, float _v) {

  for(int i=1; i<=nx; i++)
    for(int j=1; j<=ny; j++) {
      hu[i][j] = h[i][j] * _u;
      hv[i][j] = h[i][j] * _v;
    };

  synchDischargeAfterWrite();
}

/**
 * set discharge in all interior grid cells (i.e. except ghost layer) 
 * to values specified by parameter functions
 * Note: unknowns hu and hv represent momentum, while parameters u and v are velocities! 
 */
void SWE_Block::setDischarge(float (*_u)(float, float), float (*_v)(float, float)) {

  for(int i=1; i<=nx; i++)
    for(int j=1; j<=ny; j++) {
      float x = offsetX + (i-0.5f)*dx;
      float y = offsetY + (j-0.5f)*dy;
      hu[i][j] = _u(x,y) * h[i][j];
      hv[i][j] = _v(x,y) * h[i][j]; 
    };

  synchDischargeAfterWrite();
}

/**
 * set Bathymetry b in all grid cells (incl. ghost/boundary layers)
 * to a uniform value
 * bathymetry source terms are re-computed
 */
void SWE_Block::setBathymetry(float _b) {

  for(int i=0; i<=nx+1; i++)
    for(int j=0; j<=ny+1; j++)
      b[i][j] = _b;

  synchBathymetryAfterWrite();
}

/**
 * set Bathymetry b in all grid cells (incl. ghost/boundary layers)
 * using the specified bathymetry function;
 * bathymetry source terms are re-computed
 */
void SWE_Block::setBathymetry(float (*_b)(float, float)) {

  for(int i=0; i<=nx+1; i++)
    for(int j=0; j<=ny+1; j++)
      b[i][j] = _b(offsetX + (i-0.5f)*dx, offsetY + (j-0.5f)*dy);

  synchBathymetryAfterWrite();
}

// /** 
// 	Restores values for h, v, and u from file data
// 	@param _h		array holding h-values in sequence
// 	@param _v		array holding v-values in sequence
// 	@param _u		array holding u-values in sequence
// */
// void SWE_Block::setInitValues(float* _h, float* _u, float* _v) {
// 	/* Corresponding output code
// 	for (int j=1; j<ny+1;j++)
// 	for (int i=1;i<nx+1;i++)
// 		Vtk_file <<(h[i][j]+b[i][j])<<endl;
// 	*/
// 	int i, j;	
// 	for(int k=0; k<nx*ny; k++) {
// 		i = (k % ny) + 1;
// 		j = (k / ny) + 1;
// 		h[i][j] = _h[k];
// 		hu[i][j] = _h[k] * _u[k];
// 		hv[i][j] = _h[k] * _v[k]; 
// 	};
// 
// 	synchWaterHeightAfterWrite();
// 	synchDischargeAfterWrite();
// }

// /** 
// 	Restores values for h, v, and u from file data
// 	@param _b		array holding b-values in sequence
// */
// void SWE_Block::setBathymetry(float* _b) {
// 	// Set all inner cells to the value available
// 	int i, j;	
// 	for(int k=0; k<nx*ny; k++) {
// 		i = (k % ny) + 1;
// 		j = (k / ny) + 1;
// 		b[i][j] = _b[k];
// 	};
// 
// 	// Set ghost cells values such that normals = 0
// 	// Boundaries
// 	for(int i=1; i<=nx; i++) {
// 		b[0][i] = b[1][i];
// 		b[nx+1][i] = b[nx][i];
// 		b[i][0] = b[i][1];
// 		b[i][nx+1] = b[i][nx];
// 	}
// 	// Corners
// 	b[0][0] = b[1][1];
// 	b[0][ny+1] = b[1][ny];
// 	b[nx+1][0] = b[nx][1];
// 	b[nx+1][ny+1] = b[nx][ny];
// 
// 	synchBathymetryAfterWrite();
// }

/**
 * return reference to water height unknown h
 */
const Float2D& SWE_Block::getWaterHeight() { 
  synchWaterHeightBeforeRead();
  return h; 
};

/**
 * return reference to discharge unknown hu
 */
const Float2D& SWE_Block::getDischarge_hu() { 
  synchDischargeBeforeRead();
  return hu; 
};

/**
 * return reference to discharge unknown hv
 */
const Float2D& SWE_Block::getDischarge_hv() { 
  synchDischargeBeforeRead();
  return hv;
};

/**
 * return reference to bathymetry unknown b
 */
const Float2D& SWE_Block::getBathymetry() { 
  synchBathymetryBeforeRead();
  return b; 
};

//==================================================================
// methods for simulation
//==================================================================

/**
 * set wall boundary tpye for the four block boundaries
 */
void SWE_Block::setWallBoundaries() {
  
  boundary[BND_LEFT]   = WALL;
  boundary[BND_RIGHT]  = WALL;
  boundary[BND_BOTTOM] = WALL;
  boundary[BND_TOP]    = WALL;
}

/**
 * set outflow boundary tpye for the four block boundaries
 */
void SWE_Block::setOutflowBoundaries() {
  
  boundary[BND_LEFT]   = OUTFLOW;
  boundary[BND_RIGHT]  = OUTFLOW;
  boundary[BND_BOTTOM] = OUTFLOW;
  boundary[BND_TOP]    = OUTFLOW;
}

/**
 * Set the boundary type for specific block boundary.
 *
 * @param i_edge location of the edge relative to the SWE_block.
 * @param i_boundaryType type of the boundary condition.
 * @param i_inflow pointer to an SWE_Block1D, which specifies the inflow (should be NULL for WALL or OUTFLOW boundary)
 */
void SWE_Block::setBoundaryType( const BoundaryEdge i_edge,
                                 const BoundaryType i_boundaryType,
                                 const SWE_Block1D* i_inflow) {
  boundary[i_edge] = i_boundaryType;
  neighbour[i_edge] = i_inflow;

  // set bathymetry values in the ghost layer, if necessary
  for(int j=0; j<=ny+1; j++) {
    if( boundary[BND_LEFT] == OUTFLOW || boundary[BND_LEFT] == WALL ) {
      b[0][j] = b[1][j];
    }
    if( boundary[BND_RIGHT] == OUTFLOW || boundary[BND_RIGHT] == WALL ) {
      b[nx+1][j] = b[nx][j];
    }
  }
  for(int i=0; i<=nx+1; i++) {
    if( boundary[BND_BOTTOM] == OUTFLOW || boundary[BND_BOTTOM] == WALL ) {
      b[i][0] = b[i][1];
    }
    if( boundary[BND_TOP] == OUTFLOW || boundary[BND_TOP] == WALL ) {
      b[i][ny+1] = b[i][ny];
    }
  }

  // synchronize after an external update of the bathymetry
  synchBathymetryAfterWrite();
}

// /**
//  * define a CONNECT boundary:
//  * the block boundary with index egde (see method setBoundaryType()) 
//  * is connect to the boundary neighEdge of block neighBlock
//  */
// void SWE_Block::connectBoundaries(BoundaryEdge edge, SWE_Block &neighBlock, BoundaryEdge neighEdge) {
// 
//   boundary[edge] = CONNECT;
//   neighBlock.boundary[neighEdge] = CONNECT;
// 
//   neighbour[edge] = &neighBlock;
//   neighBlock.neighbour[neighEdge] = this;
// }
// 
/**
 * register the row or column layer next to a boundary as a "copy layer",
 * from which values will be copied into the ghost layer or a neighbour;
 * @return	a SWE_Block1D object that contains row variables h, hu, and hv
 */
SWE_Block1D* SWE_Block::registerCopyLayer(BoundaryEdge edge){

  switch (edge) {
    case BND_LEFT:
      return new SWE_Block1D( h.getColProxy(1), hu.getColProxy(1), hv.getColProxy(1) );
    case BND_RIGHT:
      return new SWE_Block1D( h.getColProxy(nx), hu.getColProxy(nx), hv.getColProxy(nx) );
    case BND_BOTTOM:
      return new SWE_Block1D( h.getRowProxy(1), hu.getRowProxy(1), hv.getRowProxy(1));
    case BND_TOP:
      return new SWE_Block1D( h.getRowProxy(ny), hu.getRowProxy(ny), hv.getRowProxy(ny));
  };
  return NULL;
}

/**
 * "grab" the ghost layer at the specific boundary in order to set boundary values 
 * in this ghost layer externally. 
 * The boundary conditions at the respective ghost layer is set to PASSIVE, 
 * such that the grabbing program component is responsible to provide correct 
 * values in the ghost layer, for example by receiving data from a remote 
 * copy layer via MPI communication. 
 * @param	specified edge
 * @return	a SWE_Block1D object that contains row variables h, hu, and hv
 */
SWE_Block1D* SWE_Block::grabGhostLayer(BoundaryEdge edge){

  boundary[edge] = PASSIVE;
  switch (edge) {
    case BND_LEFT:
      return new SWE_Block1D( h.getColProxy(0), hu.getColProxy(0), hv.getColProxy(0) );
    case BND_RIGHT:
      return new SWE_Block1D( h.getColProxy(nx+1), hu.getColProxy(nx+1), hv.getColProxy(nx+1) );
    case BND_BOTTOM:
      return new SWE_Block1D( h.getRowProxy(0), hu.getRowProxy(0), hv.getRowProxy(0));
    case BND_TOP:
      return new SWE_Block1D( h.getRowProxy(ny+1), hu.getRowProxy(ny+1), hv.getRowProxy(ny+1));
  };
  return NULL;
}


/**
 * set the values of all ghost cells depending on the specifed 
 * boundary conditions;
 * if the ghost layer replicates the variables of a remote SWE_Block, 
 * the values are copied
 */
void SWE_Block::setGhostLayer() {

#ifdef DBG
  cout << "Set simple boundary conditions " << endl << flush;
#endif
  // call to virtual function to set ghost layer values 
  setBoundaryConditions();

  // for a CONNECT boundary, data will be copied from a neighbouring
  // SWE_Block (via a SWE_Block1D proxy object)
  // -> these copy operations cannot be executed in GPU/accelerator memory, e.g.
  //    setBoundaryConditions then has to take care that values are copied.
  
#ifdef DBG
  cout << "Set CONNECT boundary conditions in main memory " << endl << flush;
#endif
  // left boundary
  if (boundary[BND_LEFT] == CONNECT) {
     for(int j=0; j<=ny+1; j++) {
       h[0][j] = neighbour[BND_LEFT]->h[j];
       hu[0][j] = neighbour[BND_LEFT]->hu[j];
       hv[0][j] = neighbour[BND_LEFT]->hv[j];
      };
  };
  
  // right boundary
  if(boundary[BND_RIGHT] == CONNECT) {
     for(int j=0; j<=ny+1; j++) {
       h[nx+1][j] = neighbour[BND_RIGHT]->h[j];
       hu[nx+1][j] = neighbour[BND_RIGHT]->hu[j];
       hv[nx+1][j] = neighbour[BND_RIGHT]->hv[j];
      };
  };

  // bottom boundary
  if(boundary[BND_BOTTOM] == CONNECT) {
     for(int i=0; i<=nx+1; i++) {
       h[i][0] = neighbour[BND_BOTTOM]->h[i];
       hu[i][0] = neighbour[BND_BOTTOM]->hu[i];
       hv[i][0] = neighbour[BND_BOTTOM]->hv[i];
      };
  };

  // top boundary
  if(boundary[BND_TOP] == CONNECT) {
     for(int i=0; i<=nx+1; i++) {
       h[i][ny+1] = neighbour[BND_TOP]->h[i];
       hu[i][ny+1] = neighbour[BND_TOP]->hu[i];
       hv[i][ny+1] = neighbour[BND_TOP]->hv[i];
     }
  };

#ifdef DBG
  cout << "Synchronize ghost layers (for heterogeneous memory) " << endl << flush;
#endif
  // synchronize the ghost layers (for PASSIVE and CONNECT conditions)
  // with accelerator memory
  synchGhostLayerAfterWrite();
}

/**
 * Compute the largest allowed time step for the current grid block
 * (reference implementation) depending on the current values of 
 * variables h, hu, and hv, and store this time step size in member 
 * variable maxTimestep.
 *
 * @param i_dryTol dry tolerance (dry cells do not affect the time step).
 * @param i_cflNumber CFL number of the used method.
 */
void SWE_Block::computeMaxTimestep( const float i_dryTol,
                                    const float i_cflNumber ) {
  
  // initialize the maximum wave speed
  float l_maximumWaveSpeed = (float) 0;

  // compute the maximum wave speed within the grid
  for(int i=1; i <= nx; i++) {
    for(int j=1; j <= ny; j++) {
      if( h[i][j] > i_dryTol ) {
        float l_momentum = std::max( std::abs( hu[i][j] ),
                                     std::abs( hv[i][j] ) );

        float l_particleVelocity = l_momentum / h[i][j];
        
        // approximate the wave speed
        float l_waveSpeed = l_particleVelocity + std::sqrt( g * h[i][j] );
        
        l_maximumWaveSpeed = std::max( l_maximumWaveSpeed, l_waveSpeed );
      }
    }
  }
  
  float l_minimumCellLength = std::min( dx, dy );

  // set the maximum time step variable
  maxTimestep = l_minimumCellLength / l_maximumWaveSpeed;

  // apply the CFL condition
  maxTimestep *= i_cflNumber;
}


//==================================================================
// protected member functions for simulation
// (to provide a reference implementation)
//==================================================================

/**
 * set the values of all ghost cells depending on the specifed 
 * boundary conditions
 * - set boundary conditions for typs WALL and OUTFLOW
 * - derived classes need to transfer ghost layers
 */
void SWE_Block::setBoundaryConditions() {

  // CONNECT boundary conditions are set in the calling function setGhostLayer
  // PASSIVE boundary conditions need to be set by the component using SWE_Block

  // Linker Rand
  switch(boundary[BND_LEFT]) {
    case WALL:
    {
      for(int j=1; j<=ny; j++) {
        h[0][j] = h[1][j];
        hu[0][j] = -hu[1][j];
        hv[0][j] = hv[1][j];
      };
      break;
    }
    case OUTFLOW:
    {
      for(int j=1; j<=ny; j++) {
        h[0][j] = h[1][j];
        hu[0][j] = hu[1][j];
        hv[0][j] = hv[1][j];
      };
      break;
    }
    case CONNECT:
    case PASSIVE:
      break;
    default:
      assert(false);
      break;
  };

  // Rechter Rand
  switch(boundary[BND_RIGHT]) {
    case WALL:
    {
      for(int j=1; j<=ny; j++) {
        h[nx+1][j] = h[nx][j];
        hu[nx+1][j] = -hu[nx][j];
        hv[nx+1][j] = hv[nx][j];
      };
      break;
    }
    case OUTFLOW:
    {
      for(int j=1; j<=ny; j++) {
        h[nx+1][j] = h[nx][j];
        hu[nx+1][j] = hu[nx][j];
        hv[nx+1][j] = hv[nx][j];
      };
      break;
    }
    case CONNECT:
    case PASSIVE:
      break;
    default:
      assert(false);
      break;
  };

  // Unterer Rand
  switch(boundary[BND_BOTTOM]) {
    case WALL:
    {
      for(int i=1; i<=nx; i++) {
        h[i][0] = h[i][1];
        hu[i][0] = hu[i][1];
        hv[i][0] = -hv[i][1];
      };
      break;
    }
    case OUTFLOW:
    {
      for(int i=1; i<=nx; i++) {
        h[i][0] = h[i][1];
        hu[i][0] = hu[i][1];
        hv[i][0] = hv[i][1];
      };
      break;
    }
    case CONNECT:
    case PASSIVE:
      break;
    default:
      assert(false);
      break;
  };

  // Oberer Rand
  switch(boundary[BND_TOP]) {
    case WALL:
    {
      for(int i=1; i<=nx; i++) {
        h[i][ny+1] = h[i][ny];
        hu[i][ny+1] = hu[i][ny];
        hv[i][ny+1] = -hv[i][ny];
      };
      break;
    }
    case OUTFLOW:
    {
      for(int i=1; i<=nx; i++) {
        h[i][ny+1] = h[i][ny];
        hu[i][ny+1] = hu[i][ny];
        hv[i][ny+1] = hv[i][ny];
      };
      break;
    }
    case CONNECT:
    case PASSIVE:
      break;
    default:
      assert(false);
      break;
  };

  // only required for visualisation: set values in corner ghost cells
  h[0][0] = h[1][1];
  h[0][ny+1] = h[1][ny];
  h[nx+1][0] = h[nx][1];
  h[nx+1][ny+1] = h[nx][ny];

}


//==================================================================
// protected member functions for memory model: 
// in case of temporary variables (especial in non-local memory, for 
// example on accelerators), the main variables h, hu, hv, and b 
// are not necessarily updated after each time step.
// The following methods are called to synchronise before or after 
// external read or write to the variables.
//==================================================================

/**
 * Update all temporary and non-local (for heterogeneous computing) variables
 * after an external update of the main variables h, hu, hv, and b.
 */
void SWE_Block::synchAfterWrite() {
   synchWaterHeightAfterWrite();
   synchDischargeAfterWrite();
   synchBathymetryAfterWrite();
}

/**
 * Update temporary and non-local (for heterogeneous computing) variables
 * after an external update of the water height h
 */
void SWE_Block::synchWaterHeightAfterWrite() {}

/**
 * Update temporary and non-local (for heterogeneous computing) variables
 * after an external update of the discharge variables hu and hv
 */
void SWE_Block::synchDischargeAfterWrite() {}

/**
 * Update temporary and non-local (for heterogeneous computing) variables
 * after an external update of the bathymetry b
 */
void SWE_Block::synchBathymetryAfterWrite() {}

/**
 * Update the ghost layers (only for CONNECT and PASSIVE boundary conditions)
 * after an external update of the main variables h, hu, hv, and b in the 
 * ghost layer.
 */
void SWE_Block::synchGhostLayerAfterWrite() {}

/**
 * Update all temporary and non-local (for heterogeneous computing) variables
 * before an external access to the main variables h, hu, hv, and b.
 */
void SWE_Block::synchBeforeRead() {
   synchWaterHeightBeforeRead();
   synchDischargeBeforeRead();
   synchBathymetryBeforeRead();
}

/**
 * Update temporary and non-local (for heterogeneous computing) variables
 * before an external access to the water height h
 */
void SWE_Block::synchWaterHeightBeforeRead() {}

/**
 * Update temporary and non-local (for heterogeneous computing) variables
 * before an external access to the discharge variables hu and hv
 */
void SWE_Block::synchDischargeBeforeRead() {}

/**
 * Update temporary and non-local (for heterogeneous computing) variables
 * before an external access to the bathymetry b
 */
void SWE_Block::synchBathymetryBeforeRead() {}

/**
 * Update (for heterogeneous computing) variables in copy layers
 * before an external access to the unknowns
 */
void SWE_Block::synchCopyLayerBeforeRead() {}



//==================================================================
// methods for VTK output (i.e., for visualisation)
//==================================================================

/**
 * Write a VTK file (using XML format) for visualisation using ParaView
 * -> writes unknowns h, u, v, and b of a single SWE_Block
 *    as a RECTILINEAR grid for ParaView
 */
void SWE_Block::writeVTKFileXML(string FileName, int offsetX, int offsetY) {

   synchBeforeRead();

   // VTK HEADER
   Vtk_file.open(FileName.c_str());
   Vtk_file <<"<?xml version=\"1.0\"?>"<<endl;
   Vtk_file << "<VTKFile type=\"StructuredGrid\">"<<endl;
   Vtk_file << "<StructuredGrid WholeExtent=\""
            <<offsetX<<" "<<offsetX+nx<<" "<<offsetY<<" "<<offsetY+ny<<" 0"<<" 0\">"<<endl;
   Vtk_file << "<Piece Extent=\""<<offsetX<<" "<<offsetX+nx<<" "
                                 <<offsetY<<" "<<offsetY+ny<<" 0"<<" 0\">"<<endl;
   Vtk_file << "<Points>"<<endl;
   Vtk_file <<"<DataArray NumberOfComponents=\"3\" type=\"Float64\" format=\"ascii\">"<<endl;

   //GITTER PUNKTE
   for (int j=0; j<ny+1;j++)
      for (int i=0;i<nx+1;i++)
	  Vtk_file << (offsetX+i)*dx <<" "<< (offsetY+j)*dy <<" 0"<<endl;
   Vtk_file <<"</DataArray>"<<endl;
   Vtk_file << "</Points>"<<endl;
   Vtk_file <<endl;
   Vtk_file << "<CellData>"<<endl;

   // Water surface height (h+b)
   Vtk_file << "<DataArray Name=\"H\" type=\"Float64\" format=\"ascii\">"<<endl;
   for (int j=1; j<ny+1;j++)
      for (int i=1;i<nx+1;i++)
	 Vtk_file << h[i][j]+b[i][j] << endl;
   Vtk_file << "</DataArray>"<<endl;
   // Velocities
   Vtk_file << "<DataArray Name=\"U\" type=\"Float64\" format=\"ascii\">"<<endl;
   for (int j=1; j<ny+1;j++)
      for (int i=1;i<nx+1;i++)
	 Vtk_file << ((h[i][j]>0) ? hu[i][j]/h[i][j] : 0.0 ) <<endl;
   Vtk_file << "</DataArray>"<<endl;
   Vtk_file << "<DataArray Name=\"V\" type=\"Float64\" format=\"ascii\">"<<endl;
   for (int j=1; j<ny+1;j++)
      for (int i=1;i<nx+1;i++)
	 Vtk_file << ((h[i][j]>0) ? hv[i][j]/h[i][j] : 0.0 ) <<endl;
   Vtk_file << "</DataArray>"<<endl;
   // Bathymetry
   Vtk_file << "<DataArray Name=\"B\" type=\"Float64\" format=\"ascii\">"<<endl;
   for (int j=1; j<ny+1;j++)	
      for (int i=1;i<nx+1;i++)
	 Vtk_file << b[i][j]<<endl;
   Vtk_file << "</DataArray>"<<endl;
   Vtk_file << "</CellData>"<<endl;
   Vtk_file << "</Piece>"<<endl;
   Vtk_file << "</StructuredGrid>"<<endl;
   Vtk_file << "</VTKFile>"<<endl;
   Vtk_file.close();

}


/**
 * Write a VTK file for visualisation using ParaView
 * -> writes h, u, and v unknowns of a single SWE_Block
 *    as a RECTILINEAR grid for ParaView
 */
void SWE_Block::writeVTKFile(string FileName) {

	synchBeforeRead();
	
        // VTK HEADER
	Vtk_file.open(FileName.c_str());
	Vtk_file <<"# vtk DataFile Version 2.0"<<endl;
	Vtk_file << "HPC Tutorials: Michael Bader, Kaveh Rahnema, Oliver Meister"<<endl;
	Vtk_file << "ASCII"<<endl;
	Vtk_file << "DATASET RECTILINEAR_GRID"<<endl;
	Vtk_file << "DIMENSIONS "<< nx+1<<" "<<ny+1<<" "<<"1"<<endl;
	Vtk_file <<"X_COORDINATES "<< nx+1 <<" double"<<endl;
	//GITTER PUNKTE
	for (int i=0;i<nx+1;i++)
				Vtk_file << i*dx<<endl;
	Vtk_file <<"Y_COORDINATES "<< ny+1 <<" double"<<endl;
	//GITTER PUNKTE
	for (int i=0;i<ny+1;i++)
				Vtk_file << i*dy<<endl;
	Vtk_file <<"Z_COORDINATES 1 double"<<endl;
	Vtk_file <<"0"<<endl;
	Vtk_file << "CELL_DATA "<<ny*nx<<endl;
	Vtk_file << "SCALARS H double 1"<<endl;
	Vtk_file << "LOOKUP_TABLE default"<<endl;
	//DOFS
	for (int j=1; j<ny+1;j++)
		for (int i=1;i<nx+1;i++)
			Vtk_file <<(h[i][j]+b[i][j])<<endl;
	Vtk_file << "SCALARS U double 1"<<endl;
	Vtk_file << "LOOKUP_TABLE default"<<endl;
	for (int j=1; j<ny+1;j++)	
		for (int i=1;i<nx+1;i++)
			Vtk_file << ((h[i][j]>0) ? hu[i][j]/h[i][j] : 0.0 ) <<endl;
	Vtk_file << "SCALARS V double 1"<<endl;
	Vtk_file << "LOOKUP_TABLE default"<<endl;
	for (int j=1; j<ny+1;j++)	
		for (int i=1;i<nx+1;i++)
			Vtk_file << ((h[i][j]>0) ? hv[i][j]/h[i][j] : 0.0 ) <<endl;
	Vtk_file << "SCALARS B double 1"<<endl;
	Vtk_file << "LOOKUP_TABLE default"<<endl;
	for (int j=1; j<ny+1;j++)	
		for (int i=1;i<nx+1;i++)
			Vtk_file <<b[i][j]<<endl;
	Vtk_file.close();

}

/**
 * Write a VTK file for visualisation using ParaView
 * -> writes h, u, and v unknowns of a single SWE_Block
 *    as a STRUCTURED grid for ParaView
 *    (allows 3D effect for water surface)
 */
void SWE_Block::writeVTKFile3D(string FileName) {

	synchBeforeRead();
	
	// VTK HEADER
	Vtk_file.open(FileName.c_str());
	Vtk_file <<"# vtk DataFile Version 2.0"<<endl;
	Vtk_file << "HPC Tutorials: Michael Bader, Kaveh Rahnema, Oliver Meister"<<endl;
	Vtk_file << "ASCII"<<endl;
	Vtk_file << "DATASET STRUCTURED_GRID"<<endl;
	Vtk_file << "DIMENSIONS "<< nx+1<<" "<<ny+1<<" "<<"1"<<endl;
	Vtk_file << "POINTS "<<(nx+1)*(ny+1)<<" double"<<endl;
	//GITTER PUNKTE
	for (int j=0; j<ny+1;j++)
			for (int i=0;i<nx+1;i++)
				Vtk_file << i*dx<<" "<<j*dy<<" "
				         << 0.25*(h[i][j]+h[i+1][j]+h[i][j+1]+h[i+1][j+1]
					         +b[i][j]+b[i+1][j]+b[i][j+1]+b[i+1][j+1]) 
					 <<endl;
	Vtk_file <<endl;
	Vtk_file << "CELL_DATA "<<ny*nx<<endl;
	Vtk_file << "SCALARS H double 1"<<endl;
	Vtk_file << "LOOKUP_TABLE default"<<endl;
	//DOFS
	for (int j=1; j<ny+1;j++)
		for (int i=1;i<nx+1;i++)
			//Vtk_file <<(h[i][j])<<endl;
			Vtk_file <<(h[i][j]+b[i][j])<<endl;
	Vtk_file << "SCALARS U double 1"<<endl;
	Vtk_file << "LOOKUP_TABLE default"<<endl;
	for (int j=1; j<ny+1;j++)	
		for (int i=1;i<nx+1;i++)
			Vtk_file << ((h[i][j]>0) ? hu[i][j]/h[i][j] : 0.0 ) <<endl;
	Vtk_file << "SCALARS V double 1"<<endl;
	Vtk_file << "LOOKUP_TABLE default"<<endl;
	for (int j=1; j<ny+1;j++)	
		for (int i=1;i<nx+1;i++)
			Vtk_file << ((h[i][j]>0) ? hv[i][j]/h[i][j] : 0.0 ) <<endl;
	Vtk_file << "SCALARS B double 1"<<endl;
	Vtk_file << "LOOKUP_TABLE default"<<endl;
	for (int j=1; j<ny+1;j++)	
		for (int i=1;i<nx+1;i++)
			Vtk_file <<b[i][j]<<endl;
	Vtk_file.close();

}


/*
 * end of single patch output function
 */

//==================================================================
// external class-related methods
//==================================================================

/**
 * overloads the operator << such that you can use statements as 
 * cout << swe to print the unknowns of the SWE_Block swe to 
 * the output stream cout.
 * @param os	output stream
 * @param swe	SWE_Block to print
 * @return	reference to the parameter os
 */
ostream& operator<<(ostream& os, const SWE_Block& swe) {
  
  os << "Gitterzellen: " << swe.nx << "x" << swe.ny << endl;

  cout << "Wellenhoehe:" << endl;
  for(int i=0; i<=swe.nx+1; i++) {
    for(int j=0; j<=swe.ny+1; j++) {
      os << swe.h[i][j] << "  ";
    };
    os << endl;
  };

  cout << "Geschwindigkeit in x-Richtung:" << endl;
  for(int i=0; i<=swe.nx+1; i++) {
    for(int j=0; j<=swe.ny+1; j++) {
      os << swe.hu[i][j] << "  ";
    };
    os << endl;
  };

  cout << "Geschwindigkeit in y-Richtung:" << endl;
  for(int i=0; i<=swe.nx-1; i++) {
    for(int j=0; j<=swe.ny-1; j++) {
      os << swe.hv[i][j] << "  ";
    };
    os << endl;
  };

  os << flush;

  return os;
}

