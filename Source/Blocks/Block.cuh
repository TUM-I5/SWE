/**
 * @file
 * This file is part of SWE.
 *
 * @author Michael Bader, Kaveh Rahnema, Tobias Schnabel
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
 * TODO
 */

#ifndef __SWE_BLOCKCUDA_HH
#define __SWE_BLOCKCUDA_HH

#include <cuda_runtime.h>
#include <fstream>
#include <iostream>

#include "HelperFunctions.cuh"

#include "Blocks/Block.hpp"
#include "Tools/RealType.hpp"

// constexpr int kTileSize = 32;
// constexpr int kTileSize = 16;
constexpr int kTileSize = 8;

namespace cuda
{
#if defined(ENABLE_CUDA_UMA_ALLOCATOR)
#define H_  h_.data_
#define HU_ hu_.data_
#define HV_ hv_.data_
#define B_  b_.data_
#else
#define H_ hd_
#define HU_ hud_
#define HV_ hvd_
#define B_ bd_
#endif // defined(ENABLE_CUDA_UMA_ALLOCATOR)
#define ADDRESS(x) &x

    /**
     * Block extends the base class Block towards
     * a base class for a CUDA implementation of the shallow water equations.
     * It adds the respective variables in GPU memory, and provides
     * methods for data transfer between main and GPU memory.
     */
    class Block: public Blocks::Block
    {

      public:
        // Constructor und Destructor
        Block(int nx, int ny, RealType dx, RealType dy);
        virtual ~Block();

        // object methods
        virtual RealType simulate(RealType tStart, RealType tEnd) override = 0;

        virtual void updateUnknowns(RealType dt) override = 0;

        // ---> COULD BE IMPLEMENTED TO PROVIDE A DEFAULT IMPLEMENTATION
        //     // determine maximum possible time step
        //     virtual RealType getMaxTimestep() override;

        // deliver a pointer to proxy class that represents
        // the layer that is copied to an external ghost layer
        virtual Blocks::Block1D* registerCopyLayer(BoundaryEdge edge) override;
        // "grab" the ghost layer in order to set these values externally
        virtual Blocks::Block1D* grabGhostLayer(BoundaryEdge edge) override;

        // access to CUDA variables
        /**
         *  @return	pointer to the array #hd_ (water height) in device memory
         */
        const RealType* getCUDA_waterHeight() { return hd_; };
        /**
         *  @return	pointer to the array #hb (bathymetry) in device memory
         */
        const RealType* getCUDA_bathymetry() { return bd_; };

      protected:
        // synchronisation Methods
        virtual void synchAfterWrite() override;
        virtual void synchWaterHeightAfterWrite() override;
        virtual void synchDischargeAfterWrite() override;
        virtual void synchBathymetryAfterWrite() override;
        virtual void synchGhostLayerAfterWrite() override;

        virtual void synchBeforeRead() override;
        virtual void synchWaterHeightBeforeRead() override;
        virtual void synchDischargeHuBeforeRead() override;
        virtual void synchDischargeHvBeforeRead() override;
        virtual void synchBathymetryBeforeRead() override;
        virtual void synchCopyLayerBeforeRead() override;

        // set boundary conditions in ghost layers (set boundary conditions)
        virtual void setBoundaryConditions() override;

        // define arrays for main unknowns in CUDA global memory:
        // hd_, hud_, hvd_, and bd_ are CUDA arrays corresp. to h, hu, hv, and b
        RealType* hd_;
        RealType* hud_;
        RealType* hvd_;
        RealType* bd_;

      private:
        // separate memory to hold bottom and top ghost and copy layer
        // in main memory allowing non-strided access
        RealType*        bottomLayer_;
        RealType*        topLayer_;
        Blocks::Block1D* bottomGhostLayer;
        Blocks::Block1D* bottomCopyLayer;
        Blocks::Block1D* topGhostLayer;
        Blocks::Block1D* topCopyLayer;
        // and resp. memory on the CUDA device:
        RealType* bottomLayerDevice;
        RealType* topLayerDevice;

        // helper arrays: store maximum height and velocities to determine time step
        RealType* maxhd;
        RealType* maxvd;

        const int full_grid_size_in_bytes_{-1};

      public:
        // print information about the CUDA device
        static void printDeviceInformation();

        /**
         *  Initializes the cuda device
         *  Has to be called once at the beginning.
         */
        static void init(int i_cudaDevice = 0);
        /** Cleans up the cuda device */
        static void finalize();
    };

    /**
        Return index of hd_[i][j] in linearised array
        @param i,j		x- and y-coordinate of grid cell
        @param ny		grid size in y-direction (without ghost layers)
    */
    inline __device__ int getCellCoord(int x, int y, int ny) { return x * (ny + 2) + y; }

    /**
        Return index of edge-data Fhd[i][j] or Ghd[i][j] in linearised array
        @param i,j		x- and y-coordinate of grid cell
        @param ny		grid size in y-direction (without ghost layers)
    */
    inline __device__ int getEdgeCoord(int x, int y, int ny) { return x * (ny + 1) + y; }

    /**
        Return index of a specific element in the arrays of bathymetry source terms
        @param i,j		x- and y-coordinate of grid cell
        @param ny		grid size in y-direction (without ghost layers)
    */
    inline __device__ int getBathyCoord(int x, int y, int ny) { return x * ny + y; }

    Blocks::Block* getCudaBlockInstance(RealType nx, RealType ny, RealType dx, RealType dy);
} // namespace cuda

#endif
