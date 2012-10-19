#include "SWE_VtkScenarioVisInfo.h"
#include <fstream>
#include <algorithm>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#include <math.h>

using namespace std;

/**
 *  Protected Constructor for class SWE_VtkScenarioVisInfo:
 *  called by readVtkFile to create an object of correct size;
 *  only calls constructor of SWE_VtkScenario
 *  NOTE: the provided number of cells in x and y direction includes ghost cells
 *  @param nx	cells in x direction
 *  @param ny	cells in y direction
 */
SWE_VtkScenarioVisInfo::SWE_VtkScenarioVisInfo(int nx, int ny) 
: SWE_VtkScenario(nx,ny)
{
}

/**
 *  Read VTK file which was written by SWE_Block.
 *  Use finite state automaton to parse input data.
 *  @param filename	input filename
 */
SWE_VtkScenarioVisInfo* SWE_VtkScenarioVisInfo::readVtkFile(char const * filename)
{
   // Helper variables
   string line;
   int state = 0;
   int lineCount = 0;
   int totalLines = 0;
   int offsetCount = 0;
   int dimx, dimy;
   SWE_VtkScenarioVisInfo* scene;
   float* h = NULL;
   float* u = NULL;
   float* v = NULL;
   float* b = NULL;
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
      offsetCount++;
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
               dimx = atoi(values[1].c_str());
               dimy = atoi(values[2].c_str());
               scene = new SWE_VtkScenarioVisInfo(dimx,dimy);
               h = scene->h.elemVector();
               u = scene->u.elemVector();
               v = scene->v.elemVector();
               b = scene->b.elemVector();
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
	    scene->bndPos[BND_LEFT]  = 0.0;
	    scene->bndPos[BND_RIGHT] = dimx*scene->dx;
         }
         if (line.substr(0,2) == "0 " && scene->dy <= 0.0f) {
            stringExplode(line," ", &values);
            scene->dy = (float) atof(values[1].c_str());
	    // TODO: 
	    // --> this will change, if VTK files can have other origins than (0,0)
	    scene->bndPos[BND_BOTTOM]  = 0.0;
	    scene->bndPos[BND_TOP]     = dimy*scene->dy;
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
         h[offsetCount - 1] = (float) atof(line.c_str());
         if (offsetCount >= (dimx - 1)*(dimy - 1)) {
            state = 12;
         }
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
         u[offsetCount - 1] = (float) atof(line.c_str());
         if (offsetCount >= (dimx - 1)*(dimy - 1)) {
            state = 22;
         }
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
         v[offsetCount - 1] = (float) atof(line.c_str());
         if (offsetCount >= (dimx - 1)*(dimy - 1)) {
            state = 32;
         }
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
         b[offsetCount - 1] = (float) atof(line.c_str());
         h[offsetCount - 1] = h[offsetCount - 1] - b[offsetCount - 1];
         if (offsetCount >= (dimx - 1)*(dimy - 1)) {
            state = 42;   
         }
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
      std::cout << "Parsing Error. Omitting input file." << std::endl;
      return NULL;
   }
   
   // std::cout << "Compute scaling Info for read file\n" << flush;
   scene->computeHeightAtRest();
   scene->computeScalingApproximation();
   
   // return pointer to scene, if everything went OK
   return scene;
}

/**
    Computes a first approximation of the scaling values needed
    for visualization.
	Gets called after variables have been read from VTK file;
	determines the average, mininimum and maximum values of the 
	bathymetry and water surface data.
	Uses latter values to estimate the scaling factors.
*/	
void SWE_VtkScenarioVisInfo::computeScalingApproximation() {
   // Averages of B and (B + H)
   float avgB, avgH; 
   // Minimum values
   float minB, minH;
   // Maximum values
   float maxB, maxH;
        
   int nx = h.getRows();
   int ny = h.getCols();
   int maxDim = (nx > ny) ? nx : ny;

   // first, compute averages: 
   avgH = 0.0f;
   avgB = 0.0f;
   int cnt = 0;
   for(int i=0; i<nx; i++)
      for(int j=0; j<ny; j++) 
         if (h[i][j]> 0) {
            // Update averages
            avgH += (h[i][j] + b[i][j]);
            avgB += b[i][j];
	    cnt++;
         };

   bAverage = avgB/cnt;
   hAverage = avgH/cnt;

   // then, conpute min and max values, using average values as initial guess
   minB = bAverage;
   maxB = bAverage;
   minH = hAverage;
   maxH = hAverage;
   for(int i=0; i<nx; i++)
      for(int j=0; j<ny; j++) 
         if (h[i][j]> 0) {
            float bij = b[i][j];
	    float hbij = h[i][j] + bij;
	    // Update minima
	    minH = (hbij < minH) ? hbij : minH;
            minB = (bij < minB) ? bij : minB;
            // Update maxima
	    maxH = (hbij > maxH) ? hbij : maxH;
            maxB = (bij > maxB) ? bij : maxB;
         };

   std::cout << "Computed min and max values for visualisation info: \n"
             << "h: max: " << maxH << ", min: " << minH << endl
             << "b: max: " << maxB << ", min: " << minB << endl
             << flush; 

   if (fabs(maxB - minB) > 0.0001f) {
      bScale = (maxDim/20.0f)/(maxB - minB);
   } else {
      bScale = (maxDim/20.0f);
   }
   bOffset = bScale*(bAverage - minB) + maxDim/15.0f;
   // alternative values
   bScale = minB;
   bOffset = bAverage;

   if ((maxH - minH) < 0.0001f) {
      hScale = 1.0f/(maxH- minH);
   } else {
      hScale = 1.0f;
   }
   hScale = maxDim*(hScale);
   hOffset = bOffset + (maxDim/5.0f);
   
   bScale = minH;
   hOffset = hAverage;

   std::cout << "Computed values for visualisation info: "
             << "average h: " << hAverage << "( bOffset:" << hOffset << ", scale: " << hScale << ");" 
             << "average b: " << bAverage << "( bOffset:" << bOffset << ", scale: " << bScale << ");"
             << endl << flush; 
}
