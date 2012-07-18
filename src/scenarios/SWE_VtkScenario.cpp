#include "SWE_VtkScenario.h"
#include <fstream>
#include <algorithm>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

/**
 *  Protected Constructor for class SWE_VtkScenario:
 *  called by readVtkFile to create an object of correct size
 *  NOTE: the provided number of cells in x and y direction includes ghost cells
 *  @param nx	cells in x direction
 *  @param ny	cells in y direction
 */
SWE_VtkScenario::SWE_VtkScenario(int nx, int ny) 
: h(nx,ny),u(nx,ny),v(nx,ny),b(nx,ny)
{
cout << "Create VtkScenario of size " << nx << 'x' << ny << endl; 
   // set illegal value for dx and dy
   // -> to be overwritten by readVtkFile
   dx = -1.0; dy = -1.0;
}

/**
 *  get water height at specified position (piecewise constant interpolation);
 *  @param x	position x-coordinate
 *  @param y	position y-coordinate
 */
float SWE_VtkScenario::getWaterHeight(float x, float y)
{
   if (x<0 || x >1) std::cout << "ERROR in x coordinate: " << x << endl << flush;
   if (y<0 || y >1) std::cout << "ERROR in y coordinate: " << y << endl << flush;
   int i = int((x-bndPos[BND_LEFT])/dx);
   int j = int((y-bndPos[BND_BOTTOM])/dy);
   if (i<0 || i >h.getCols()) 
      std::cout << "ERROR in i coordinate: " << i 
                << " with x coordinate: " << x 
                << " offset " << (x-bndPos[BND_LEFT]) << '/' << dx << endl << flush;
   if (j<0 || j >h.getRows()) std::cout << "ERROR in j coordinate: " << j 
                << " with y coordinate: " << y << endl << flush;
   return h[i][j];
}

/**
 *  get velocity (x-direction) at specified position (piecewise constant interpolation);
 *  @param x	position x-coordinate
 *  @param y	position y-coordinate
 */
float SWE_VtkScenario::getVeloc_u(float x, float y)
{
   int i = int((x-bndPos[BND_LEFT])/dx);
   int j = int((y-bndPos[BND_BOTTOM])/dy);
   return u[i][j];
}

/**
 *  get velocity (x-direction) at specified position (piecewise constant interpolation);
 *  @param x	position x-coordinate
 *  @param y	position y-coordinate
 */
float SWE_VtkScenario::getVeloc_v(float x, float y)
{
   int i = int((x-bndPos[BND_LEFT])/dx);
   int j = int((y-bndPos[BND_BOTTOM])/dy);
   return v[i][j];
}

/**
 *  get bathymetry at specified position (piecewise constant interpolation);
 *  @param x	position x-coordinate
 *  @param y	position y-coordinate
 */
float SWE_VtkScenario::getBathymetry(float x, float y)
{
   int i,j;
   
   if (x<0) {
      std::cout << "request left of boundary" << x << endl << flush;
      i = 0;
   } else if (x >1) {
      std::cout << "request right of boundary" << x << endl << flush;
      i = b.getCols()-1;
   } else {
      i = int((x-bndPos[BND_LEFT])/dx);
   };
   
   if (y<0) {
      std::cout << "request below boundary" << y << endl << flush;
      j = 0;
   } else if (y>1) {
      std::cout << "request above boundary" << y << endl << flush;
      j = b.getRows()-1;
   } else {
      j = int((y-bndPos[BND_BOTTOM])/dy);
   };

   return b[i][j];
}

/**
 *  Read VTK file which was written by SWE_Block.
 *  Use finite state automaton to parse input data.
 *  @param filename	input filename
 */
SWE_VtkScenario* SWE_VtkScenario::readVtkFile(char const * filename)
{
   // Helper variables
   string line;
   int state = 0;
   int lineCount = 0;
   int totalLines = 0;
   int offsetCount = 0;
   int dimx, dimy;
   SWE_VtkScenario* scene;
   float* h = NULL;
   float* u = NULL;
   float* v = NULL;
   float* b = NULL;
   int i,j;
   int lastElem = 0;
   size_t found;
   vector<string> values;

   // Open file
   ifstream file(filename, ios::in);
   if (!file.good()) {
      std::cout << "Couldn't open file: " << filename << "\n";
      return NULL;
   }
   // Count lines
   totalLines = std::count(std::istreambuf_iterator<char>(file), 
            std::istreambuf_iterator<char>(), '\n');
   file.seekg (0, ios::beg);

   // Read line by line
   while (std::getline(file,line)) {
      // Skip empty lines
      lineCount++;
      if (! line.length()) continue;
      
      // Convert to lower case
      std::transform(line.begin(), line.end(),line.begin(), ::tolower);
      
      // Output progress
      if (lineCount % 10000 == 0 || lineCount == totalLines) {
         std::cout << "\r" << "Reading input file: " << 
            (lineCount*100)/totalLines << "% completed.       " << flush;
      }

      switch(state)
      {
      // Header part
      case 0: 
         found = line.find("# vtk datafile");
         if (found!=string::npos) {
            state = 1;
         }
         break;
      case 1:
         state = 2;
         break;
      case 2:
         if (line.find("ascii") == 0) {
            state = 3;
         }
         break;
      case 3:
         if (line.find("dataset structured_grid") == 0) {
            state = 4;
         }
         break;
      case 4:
         found = line.find("dimensions");
         if (found!=string::npos) {
            stringExplode(line," ", &values);
            if (values.size() != 4) {
               state = -1;
            } else {
               dimx = atoi(values[1].c_str()) - 1;
               dimy = atoi(values[2].c_str()) - 1;
cout << "Read grid dimensions " << dimx << 'x' << dimy << endl;
               scene = new SWE_VtkScenario(dimx,dimy);
	       h = scene->h.elemVector();
	       u = scene->u.elemVector();
	       v = scene->v.elemVector();
	       b = scene->b.elemVector();
               lastElem = dimx*dimy-1;
               state = 5;
            }
         } 
         break;
      // Data part
      case 5:
         // read grid data stored in VTK file
         found = line.find(" 0 ");
         if (found != string::npos && scene->dx <= 0.0f) {
            stringExplode(line," ", &values);
            scene->dx = (float) atof(values[0].c_str());
	    // TODO: 
	    // --> this will change, if VTK files can have other origins than (0,0)
	    scene->bndPos[BND_BOTTOM] = 0.0;
	    scene->bndPos[BND_TOP]    = dimx*scene->dx;
         }
         if (line.substr(0,2) == "0 " && scene->dy <= 0.0f) {
            stringExplode(line," ", &values);
            scene->dy = (float) atof(values[1].c_str());
	    // TODO: 
	    // --> this will change, if VTK files can have other origins than (0,0)
	    scene->bndPos[BND_LEFT]  = 0.0;
	    scene->bndPos[BND_RIGHT] = dimy*scene->dy;
          }

         found = line.find("scalars h");
         if (found != string::npos) {
            state = 10;
         } 
         break;
      case 10:
         found = line.find("lookup_table");
         if (found!=string::npos) {
            state = 11;
            offsetCount = 0;
         } 
         break;
      case 11:
         // Read h-values
         i = offsetCount%dimx; j= offsetCount/dimx;
         h[i*dimy+j] = (float) atof(line.c_str());
         if (offsetCount < lastElem) 
            offsetCount++;
         else
            state = 12;
         break;
      case 12:
         found = line.find("scalars u");
         if (found!=string::npos) {
            state = 20;
         } 
         break;
      case 20:
         found = line.find("lookup_table");
         if (found!=string::npos) {
            state = 21;
            offsetCount = 0;
         } 
         break;
      case 21:
         // Read u-values
         i = offsetCount%dimx; j= offsetCount/dimx;
         u[i*dimy+j] = (float) atof(line.c_str());
         if (offsetCount < lastElem) 
            offsetCount++;
         else
            state = 22;
         break;
      case 22:
         found = line.find("scalars v");
         if (found!=string::npos) {
            state = 30;
         } 
         break;
      case 30:
         found = line.find("lookup_table");
         if (found!=string::npos) {
            state = 31;
            offsetCount = 0;
         } 
         break;
      case 31:
         // Read v-values
         i = offsetCount%dimx; j= offsetCount/dimx;
         v[i*dimy+j] = (float) atof(line.c_str());
         if (offsetCount < lastElem) 
            offsetCount++;
         else
            state = 32;
         break;
      case 32:
         found = line.find("scalars b");
         if (found!=string::npos) {
            state = 40;
         } 
         break;
      case 40:
         found = line.find("lookup_table");
         if (found!=string::npos) {
            state = 41;
            offsetCount = 0;
         } 
         break;
      case 41:
         // Read b-values
         i = offsetCount%dimx; j= offsetCount/dimx;
         i = i*dimy+j;
         b[i] = (float) atof(line.c_str());
         h[i] = h[i] - b[i];
         if (offsetCount < lastElem) 
            offsetCount++;
         else
            state = 42;   
         break;
      default:
         break;
      }

   }
   std::cout << std::endl;

   // Close file
   file.close();
   // delete allocated scene object
   if (state >= 5 && state < 42)
      delete scene;
   // return NULL pointer in case of error
   if (state < 42) {
      std::cout << "Parsing Error " << state << " -> Omitting input file." << std::endl;
      return NULL;
   }

cout << "Read input data into scenario: " << endl 
     << (*scene) << endl;  
   scene->computeHeightAtRest();
   
   // return pointer to scene, if everything went OK
   return scene;
}

/**
 *  Computes an estimate of the water height at rest 
 *  (considering water height h plus bathymetry b)
 *  -> average value of h+b (considering only cells with h>0).
 *  Value is stored in member variable heightAtRest
 */
void SWE_VtkScenario::computeHeightAtRest()
{
   float havg = 0.0;
   int cnt = 0;
   for(int i=0;i<h.getRows();i++)
      for(int j=0;j<h.getCols();j++)
         if (h[i][j] > 0.0f) {
            havg += h[i][j] + b[i][j];
	    cnt++;
         };

   heightAtRest = havg/cnt;
}

/**
    Split a string by  seperator into substrings and return them as
	an array
	@param str			string to split
	@param separator	delimiter string
	@return	results		array containing all substrings 

*/	
void SWE_VtkScenario::stringExplode(string str, string separator, vector<string>* results){
    size_t found;
	results->clear();
    found = str.find_first_of(separator);
    while(found != string::npos){
        if(found > 0){
            results->push_back(str.substr(0,found));
        }
        str = str.substr(found+1);
        found = str.find_first_of(separator);
    }
    if(str.length() > 0){
        results->push_back(str);
    }
}

//==================================================================
// external class-related methods
//==================================================================

/**
 * overload operator<< such that data can be written via cout <<
 * -> needs to be declared as friend to be allowed to access private data
 */
ostream& operator<<(ostream& os, SWE_VtkScenario& scene) {
  
  int nx = scene.h.getCols();
  int ny = scene.h.getRows();
  
  os << "Gitterzellen: " << nx << "x" << ny << endl;

  cout << "Wellenhoehe:" << endl;
  for(int i=0; i<nx; i++) {
    for(int j=0; j<ny; j++) {
      os << scene.h[i][j] << "  ";
    };
    os << endl;
  };

  cout << "Geschwindigkeit in x-Richtung:" << endl;
  for(int i=0; i<nx; i++) {
    for(int j=0; j<ny; j++) {
      os << scene.u[i][j] << "  ";
    };
    os << endl;
  };

  cout << "Geschwindigkeit in y-Richtung:" << endl;
  for(int i=0; i<nx; i++) {
    for(int j=0; j<ny; j++) {
      os << scene.v[i][j] << "  ";
    };
    os << endl;
  };

  cout << "Bathymetry:" << endl;
  for(int i=0; i<nx; i++) {
    for(int j=0; j<ny; j++) {
      os << scene.b[i][j] << "  ";
    };
    os << endl;
  };

  os << flush;

  return os;
}

