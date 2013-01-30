/**
 * @file
 * This file is part of SWE.
 *
 * @author Michael Bader, Kaveh Rahnema, Tobias Schnabel
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

#include "SWE_RusanovBlock.hh"
#include <math.h>

/**
 * Constructor: allocate variables for simulation
 *
 * unknowns h,hu,hv,b are defined on grid indices [0,..,nx+1]*[0,..,ny+1]
 * -> computational domain is [1,..,nx]*[1,..,ny]
 * -> plus ghost cell layer
 *
 * flux terms are defined for edges with indices [0,..,nx]*[1,..,ny]
 * or [1,..,nx]*[0,..,ny] (for horizontal/vertical edges)
 * Flux term with index (i,j) is located on the edge between 
 * cells with index (i,j) and (i+1,j) or (i,j+1)
 *
 * bathymetry source terms are defined for cells with indices [1,..,nx]*[1,..,ny]
 *
 * @ param _offsetX	x coordinate of block origin
 * @ param _offsetY	y coordinate of block origin
 */
SWE_RusanovBlock::SWE_RusanovBlock(float _offsetX, float _offsetY) 
: SWE_Block(_offsetX,_offsetY),
  Fh(nx+1,ny+1), Fhu(nx+1,ny+1), Fhv(nx+1,ny+1),
  Gh(nx+1,ny+1), Ghu(nx+1,ny+1), Ghv(nx+1,ny+1),
  Bx(nx+1,ny+1), By(nx+1,ny+1)
{
}

/**
 * Destructor: de-allocate all variables
 */
SWE_RusanovBlock::~SWE_RusanovBlock() {
}


/**
 * implements interface function updateUnknowns:
 * based on the (Rusanov) fluxes computed on each edge
 * (and stored in the variables Fh, Gh, etc.);
 * compute the balance terms for each cell, and update the 
 * unknowns according to an Euler time step.
 * @param dt	size of the time step.
 */
void SWE_RusanovBlock::updateUnknowns(float dt) {

  for(int i=1; i<=nx; i++)
    for(int j=1; j<=ny; j++) {
      h[i][j] -= dt *( (Fh[i][j]-Fh[i-1][j])/dx + (Gh[i][j]-Gh[i][j-1])/dy );
      
      if (h[i][j] > 0) {
	 hu[i][j] -= dt *( (Fhu[i][j]-Fhu[i-1][j])/dx + (Ghu[i][j]-Ghu[i][j-1])/dy  
		            + Bx[i][j]/dx );
	 hv[i][j] -= dt *( (Fhv[i][j]-Fhv[i-1][j])/dx + (Ghv[i][j]-Ghv[i][j-1])/dy  
		            + By[i][j]/dy );
      } else {
         // set all unknowns to 0, if h turns out non-positive
         h[i][j] = 0.0;
         hu[i][j] = 0.0;
         hv[i][j] = 0.0;
      };
    };
}

/** 
 * Depending on the current values of h, hu, hv (incl. ghost layers)  
 * update these unknowns in each grid cell
 * (ghost layers and bathymetry are not updated).
 * The Rusanov implementation of simulateTimestep 
 * subsequently calls the functions computeNumericalFluxes (to compute 
 * all fluxes on grid edges), and updateUnknowns (to update the variables 
 * according to flux values, typically according to an Euler time step). 
 * @param dt	size of the time step
 */
void SWE_RusanovBlock::simulateTimestep(float dt) {

  // compute the numerical fluxes on all edges
  computeNumericalFluxes();

  // update the unknowns according to an Euler timestep 
  updateUnknowns(dt);

}

/**
 * implements interface function simulate:
 * perform forward-Euler time steps, starting with simulation time tStart,:
 * until simulation time tEnd is reached; 
 * boundary conditions and bathymetry source terms are computed for each 
 * timestep as required -
 * intended as main simulation loop between two checkpoints
 */
float SWE_RusanovBlock::simulate(float tStart, float tEnd) {
  float t = tStart;
  do {
     // set values in ghost cells:
     setGhostLayer();
     
     // execute Euler time step:
     simulateTimestep(maxTimestep);

     t += maxTimestep; cout << "Simulation at time " << t << endl << flush;

     // calculate and set largest allowed time step:
     computeMaxTimestep();
     
  } while(t < tEnd);

  return t;
}


/**
 * compute the flux terms on all edges;
 * before the computation, computeBathymetrySources is called
 */
void SWE_RusanovBlock::computeNumericalFluxes() {

  // compute bathymetry source terms (depend on h)
  computeBathymetrySources();

  // fluxes in x direction:
  for(int i=0; i<=nx; i++)
    for(int j=1; j<=ny; j++) {
      
//       if ( h[i+1][j]==0 && h[i][j]+b[i][j] < b[i+1][j] ) {
//          Fh[i][j] = 0.0;
// 	 Fhu[i][j] = 0.5*g*h[i][j]*(h[i][j]+b[i][j]);
// 	 Fhv[i][j] = 0.0;
// 	 continue;
//       };
//       if ( h[i][j]==0 && b[i][j] > h[i+1][j]+b[i+1][j] ) {
//          Fh[i][j] = 0.0;
// 	 Fhu[i][j] = 0.5*g*h[i+1][j]*(h[i+1][j]+b[i+1][j]);
// 	 Fhv[i][j] = 0.0;
// 	 continue;
//       };
     
      double uij = (h[i][j] > 0) ? hu[i][j]/h[i][j] : 0.0;
      double uip1j = (h[i+1][j] > 0) ? hu[i+1][j]/h[i+1][j] : 0.0;
     
      double llf = computeLocalSV(i,j,'x');
      double upwind = ( fabs(uij) > fabs(uip1j) ) ? fabs(uij) : fabs(uip1j);
      
      // h-Komponente
      Fh[i][j] = computeFlux( hu[i][j], hu[i+1][j], h[i][j], h[i+1][j], upwind );

//       double hhpb = (h[i+1][j] > 0) ? h[i+1][j]*(h[i+1][j]+b[i+1][j]) : h[i][j]*(h[i][j]+b[i][j]);
// 
//       // hu-Komponente
//       Fhu[i][j] = computeFlux(hu[i][j]*uij + 0.5*g*h[i][j]*(h[i][j]+b[i][j]),
//                               hu[i+1][j]*uip1j + 0.5*g*hhpb,
// 			      hu[i][j], 
// 			      hu[i+1][j], llf );

      Fhu[i][j] = computeFlux(hu[i][j]*uij + 0.5*g*h[i][j]*h[i][j],
                              hu[i+1][j]*uip1j + 0.5*g*h[i+1][j]*h[i+1][j],
			      hu[i][j], 
			      hu[i+1][j], llf );
      Fhv[i][j] = computeFlux( uij*hv[i][j],uip1j*hv[i+1][j], hv[i][j], hv[i+1][j], llf );
    };
  
  // fluxes in y direction:
  for(int j=0; j<=ny; j++)
    for(int i=1; i<=nx; i++) {

//       if ( h[i][j+1]==0 && h[i][j]+b[i][j] < b[i][j+1] ) {
//          Gh[i][j] = 0.0;
// 	 Ghu[i][j] = 0.0;
// 	 Ghv[i][j] = 0.5*g*h[i][j]*(h[i][j]+b[i][j]);
// 	 continue;
//       };
//       if ( h[i][j]==0 && b[i][j] > h[i][j+1]+b[i][j+1] ) {
//          Gh[i][j] = 0.0;
// 	 Ghu[i][j] = 0.0;
// 	 Ghv[i][j] = 0.5*g*h[i][j]*(h[i][j+1]+b[i][j+1]);
// 	 continue;
//       };
     
      double vij = (h[i][j] > 0) ? hv[i][j]/h[i][j] : 0.0;
      double vijp1 = (h[i][j+1] > 0) ? hv[i][j+1]/h[i][j+1] : 0.0;
     
      double llf = computeLocalSV(i,j,'y');
      double upwind = ( fabs(vij) > fabs(vijp1) ) ? fabs(vij) : fabs(vijp1);

      // h-Konponente
      Gh[i][j] = computeFlux( hv[i][j], hv[i][j+1], h[i][j], h[i][j+1], upwind );

//       double hhpb = (h[i][j+1] > 0) ? h[i][j+1]*(h[i][j+1]+b[i][j+1]) : h[i][j]*(h[i][j]+b[i][j]);
//       Ghu[i][j] = computeFlux(hu[i][j]*vij,hu[i][j+1]*vijp1, hu[i][j], hu[i][j+1], llf );
//       Ghv[i][j] = computeFlux(hv[i][j]*vij + 0.5*g*h[i][j]*(h[i][j]+b[i][j]),
//                               hv[i][j+1]*vijp1 + 0.5*g*hhpb,
// 			      hv[i][j], 
// 			      hv[i][j+1], llf );

      Ghu[i][j] = computeFlux(hu[i][j]*vij,hu[i][j+1]*vijp1, hu[i][j], hu[i][j+1], llf );
      Ghv[i][j] = computeFlux(hv[i][j]*vij + 0.5*g*h[i][j]*h[i][j],
                              hv[i][j+1]*vijp1 + 0.5*g*h[i][j+1]*h[i][j+1],
			      hv[i][j], 
			      hv[i][j+1], llf );
    };
}

/**
 * compute the flux term on a given edge 
 * (acc. to local Lax-Friedrich method aka Rusanov flux):
 * fLow and fHigh contain the values of the flux function in the two 
 * adjacent grid cells 
 * xiLow and xiHigh are the values of the unknowns in the two 
 * adjacent grid cells 
 * "Low" represents the cell with lower i/j index ("High" for larger index). 
 * llf should contain the local signal velocity (as compute by computeLocalSV)
 * for llf=dx/dt (or dy/dt), we obtain the standard Lax Friedrich method
 */
float SWE_RusanovBlock::computeFlux(float fLow, float fHigh, float xiLow, float xiHigh,
                                    float llf) {
  // Rusanov / local Lax-Friedrich
  return 0.5*(fLow+fHigh) - 0.5*llf*(xiHigh-xiLow);
}

/**
 * computes the local signal velocity in x- or y-direction
 * for two adjacent cells with indices (i,j) and (i+1,j)
 * (if dir='x') or (i,j+1) (if dir='y'
 */
float SWE_RusanovBlock::computeLocalSV(int i, int j, char dir) {
  float sv1, sv2;
  
  if (dir=='x') {
     sv1 = (h[i][j] > 0.0) ? (fabs(hu[i][j]/h[i][j]) + sqrt(g*h[i][j])) : 0.0;
     sv2 = (h[i+1][j] > 0.0) ? (fabs(hu[i+1][j]/h[i+1][j]) + sqrt(g*h[i+1][j])) : 0.0;
  } else {
     sv1 = (h[i][j] > 0.0) ? (fabs(hv[i][j]/h[i][j]) + sqrt(g*h[i][j])) : 0.0;
     sv2 = (h[i][j+1] > 0.0) ? (fabs(hv[i][j+1]/h[i][j+1]) + sqrt(g*h[i][j+1])) : 0.0;
  };
  
  return (sv1 > sv2) ? sv1 : sv2;
}

/**
 * compute the bathymetry source terms in all cells
 */
void SWE_RusanovBlock::computeBathymetrySources() {

  for(int i=1; i<=nx; i++)
    for(int j=1; j<=ny; j++) {

      if ( h[i][j] <= 0 ) {
         Bx[i][j] = 0.0; By[i][j] = 0.0;
	 continue;
      };
// 
//       if ( h[i-1][j]==0 && h[i][j]+b[i][j] < b[i-1][j] )
//          Bx[i][j] = g * 0.25f*h[i+1][j]*b[i+1][j];
//          // Bx[i][j] = g * 0.5f*(h[i+1][j]+h[i-1][j]) * (b[i+1][j] - b[i-1][j]);
//       else 
//          Bx[i][j] = g * 0.5f*(h[i+1][j]+h[i-1][j]) * 0.5f*(b[i+1][j] - b[i-1][j]);
//       // Bx[i][j] = g * h[i][j] * 0.5f*(b[i+1][j] - b[i-1][j]);
// 
//       if ( h[i+1][j]==0 && h[i][j]+b[i][j] < b[i+1][j] )
//          Bx[i][j] -= g * 0.25f*h[i-1][j]*b[i-1][j];
//          // Bx[i][j] = g * 0.5f*(h[i][j]+h[i-1][j]) * (b[i][j] - b[i-1][j]);
//       else 
//           Bx[i][j] = g * 0.5f*(h[i+1][j]+h[i-1][j]) * 0.5f*(b[i+1][j] - b[i-1][j]);
     

      if ( h[i+1][j]==0 && h[i][j]+b[i][j] < b[i+1][j] )
         Bx[i][j] = - g * 0.25f*h[i-1][j]*b[i-1][j];
         // Bx[i][j] = g * 0.5f*(h[i][j]+h[i-1][j]) * (b[i][j] - b[i-1][j]);
      else if ( h[i-1][j]==0 && h[i][j]+b[i][j] < b[i-1][j] )
         Bx[i][j] = g * 0.25f*h[i+1][j]*b[i+1][j];
         // Bx[i][j] = g * 0.5f*(h[i+1][j]+h[i-1][j]) * (b[i+1][j] - b[i-1][j]);
      else 
         Bx[i][j] = g * 0.5f*(h[i+1][j]+h[i-1][j]) * 0.5f*(b[i+1][j] - b[i-1][j]);
      // Bx[i][j] = g * h[i][j] * 0.5f*(b[i+1][j] - b[i-1][j]);
      
      if ( h[i][j+1]==0 && h[i][j]+b[i][j] < b[i][j+1] )
         By[i][j] = g * 0.25f*h[i][j-1]*b[i][j-1];
         // By[i][j] = g * 0.5f*(h[i][j]+h[i][j-1]) * (b[i][j] - b[i][j-1]);
      else if ( h[i][j-1]==0 && h[i][j]+b[i][j] < b[i][j-1] )
         By[i][j] = g * 0.25f*h[i][j+1]*b[i][j+1];
         // By[i][j] = g * 0.5f*(h[i][j+1]+h[i][j]) * (b[i][j+1] - b[i][j]);
      else 
         By[i][j] = g * 0.5f*(h[i][j+1]+h[i][j-1]) * 0.5f*(b[i][j+1] - b[i][j-1]);
      // By[i][j] = g * h[i][j] * 0.5f*(b[i][j+1] - b[i][j-1]);
  };
}



//==================================================================

/**
 * overload operator<< such that data can be written via cout <<
 * -> needs to be declared as friend to be allowed to access private data
 */
ostream& operator<<(ostream& os, const SWE_RusanovBlock& swe) {
  
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

  cout << "Fluss - Wellenhoehe:" << endl;
  for(int i=0; i<=swe.nx; i++) {
    for(int j=1; j<=swe.ny; j++) {
      os << swe.Fh[i][j] << "  ";
    };
    os << endl;
  };

  cout << "Fluss - Durchfluss in x-Richtung:" << endl;
  for(int i=0; i<=swe.nx; i++) {
    for(int j=1; j<=swe.ny; j++) {
      os << swe.Fhu[i][j] << "  ";
    };
    os << endl;
  };

  cout << "Fluss - Durchfluss in y-Richtung:" << endl;
  for(int i=0; i<=swe.nx; i++) {
    for(int j=1; j<=swe.ny; j++) {
      os << swe.Fhv[i][j] << "  ";
    };
    os << endl;
  };

  cout << "Fluss - Wellenhoehe:" << endl;
  for(int i=1; i<=swe.nx; i++) {
    for(int j=0; j<=swe.ny; j++) {
      os << swe.Gh[i][j] << "  ";
    };
    os << endl;
  };

  cout << "Fluss - Durchfluss in x-Richtung:" << endl;
  for(int i=1; i<=swe.nx; i++) {
    for(int j=0; j<=swe.ny; j++) {
      os << swe.Ghu[i][j] << "  ";
    };
    os << endl;
  };

  cout << "Fluss - Durchfluss in y-Richtung:" << endl;
  for(int i=1; i<=swe.nx; i++) {
    for(int j=0; j<=swe.ny; j++) {
      os << swe.Ghv[i][j] << "  ";
    };
    os << endl;
  };
  
  os << flush;

  return os;
}

