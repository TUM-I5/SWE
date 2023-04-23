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

#include "Block.cuh"
#include "BlockKernels.cuh"

#if defined(WITH_SOLVER_FWAVE) || defined(WITH_SOLVER_AUGRIE) || defined(WITH_SOLVER_HLLE)
#include "WavePropagationBlock.cuh"
#elif defined(WITH_SOLVER_RUSANOV)
#include "SWE_RusanovBlockCUDA.hh"
#endif

#include <cassert>
#include <cmath>
#include <cstdlib>

#include "Tools/Logger.hpp"
#include "Tools/RealType.hpp"

using namespace std;

// #include "SWE_BlockCUDA_kernels.cu"

namespace cuda
{

    Blocks::Block* getCudaBlockInstance(RealType nx_, RealType ny_, RealType dx, RealType dy)
    {
#if defined(WITH_SOLVER_FWAVE) || defined(WITH_SOLVER_AUGRIE)
        Blocks::Block* block = new cuda::WavePropagationBlock(nx_, ny_, dx, dy);
#else
#error "Currently, only the fWave solver is supported!"
#endif
        return block;
    }

    /**
     * Constructor: allocate variables for simulation
     *
     * unknowns h_,hu_,hv_,b are defined on grid indices [0,..,nx_+1]*[0,..,ny_+1]
     * -> computational domain is [1,..,nx_]*[1,..,ny_]
     * -> plus ghost cell layer
     *
     * flux terms are defined for edges with indices [0,..,nx_]*[1,..,ny_]
     * or [1,..,nx_]*[0,..,ny_] (for horizontal/vertical edges)
     * Flux term with index (i,j) is located on the edge between
     * cells with index (i,j) and (i+1,j) or (i,j+1)
     *
     * bathymetry source terms are defined for cells with indices [1,..,nx_]*[1,..,ny_]
     *
     */
    Block::Block(int nx, int ny, RealType dx, RealType dy):
        Blocks::Block(nx, ny, dx, dy),
        full_grid_size_in_bytes_((nx_ + 2) * (ny_ + 2) * sizeof(RealType))
    {
        if (nx_ % kTileSize != 0)
        {
            LOG_WAR << "WARNING: nx_ not a multiple of kTileSize  -> will lead to crashes!" << endl << flush;
        };
        if (ny_ % kTileSize != 0)
        {
            LOG_WAR << "WARNING: ny_ not a multiple of kTileSize  -> will lead to crashes!" << endl << flush;
        };

       // allocate CUDA memory for unknows h_,hu_,hv_ and bathymetry b
        CUDA_MALLOC((void**)ADDRESS(H_), full_grid_size_in_bytes_);
        CUDA_MALLOC((void**)ADDRESS(HU_), full_grid_size_in_bytes_);
        CUDA_MALLOC((void**)ADDRESS(HV_), full_grid_size_in_bytes_);
        CUDA_MALLOC((void**)ADDRESS(B_), full_grid_size_in_bytes_);

        int size = nx_ + 2;
        // allocate consecutive memory for 2 columns with three unknowns each
        // (h_, hu_, hv_, excluding b) for copy/ghost layer at bottom/top boundary_
        bottomLayer_     = new RealType[6 * size];
        bottomGhostLayer = new Blocks::Block1D(bottomLayer_, bottomLayer_ + size, bottomLayer_ + (2 * size), size);
        bottomCopyLayer  = new Blocks::Block1D(
            bottomLayer_ + (3 * size), bottomLayer_ + (4 * size), bottomLayer_ + (5 * size), size
        );

        // same for top boundary_:
        topLayer_     = new RealType[6 * size];
        topGhostLayer = new Blocks::Block1D(topLayer_, topLayer_ + size, topLayer_ + (2 * size), size);
        topCopyLayer  = new Blocks::Block1D(
            topLayer_ + (3 * size), topLayer_ + (4 * size), topLayer_ + (5 * size), size
        );

        // allocate resp. memory on the CUDA device
        CUDA_MALLOC((void**)&bottomLayerDevice, 6 * size * sizeof(RealType));
        CUDA_MALLOC((void**)&topLayerDevice, 6 * size * sizeof(RealType));

        init(0); // FIXME[epic=SWE,seq=54] print for all gpus
    }

    /**
     * Destructor: de-allocate all variables
     */
    Block::~Block()
    {
        CUDA_FREE(H_);
        CUDA_FREE(HU_);
        CUDA_FREE(HV_);
        CUDA_FREE(B_);
        //  CUDA_FREE(maxhd); CUDA_FREE(maxvd);

        // free memory for copy/ghost layers
        delete bottomLayer_;
        delete topLayer_;

        CUDA_FREE(bottomLayerDevice);
        CUDA_FREE(topLayerDevice);
    }

    //==================================================================
    // methods for simulation
    //==================================================================

    /**
     * set the values of all ghost cells depending on the specifed
     * boundary_ conditions
     */
    void Block::setBoundaryConditions()
    {
#ifdef DBG
        cout << "Call kernel to compute h_ in ghost layer corner (for visualisation only) " << flush << endl;
#endif
        // Fill ghost layer corner cells
        kernelHdBufferEdges<<<1, 1>>>(H_, nx_, ny_);
        kernelHdBufferEdges<<<1, 1>>>(HU_, nx_, ny_);
        kernelHdBufferEdges<<<1, 1>>>(HV_, nx_, ny_);

#ifdef DBG
        cout << "Call kernel to compute left/right boundaries " << flush << endl;
#endif

#ifdef ENABLE_MPI
        // Sync from device to send to neighbors
        synchAfterWrite();

#if !defined(ENABLE_CUDA_UMA_ALLOCATOR)

        synchGhostLayerAfterWrite();
#endif // ENABLE_CUDA_UMA_ALLOCATOR

#endif // ENABLE_MPI
        if (boundary_[BoundaryEdge::Left] == BoundaryType::Passive
            || boundary_[BoundaryEdge::Left] == BoundaryType::Connect)
        {
            // nothing to be done:
            // ghost values are copied by Block::synchGhostLayerAfterWrite(...)
        }
        else
        {
            dim3 dimBlock(1, kTileSize);
            dim3 dimGrid(1, ny_ / kTileSize);
            kernelLeftBoundary<<<dimGrid, dimBlock>>>(H_, HU_, HV_, nx_, ny_, boundary_[BoundaryEdge::Left]);
        };

        if (boundary_[BoundaryEdge::Right] == BoundaryType::Passive
            || boundary_[BoundaryEdge::Right] == BoundaryType::Connect)
        {
            // nothing to be done:
            // ghost values are copied by Block::synchGhostLayerAfterWrite(...)
        }
        else
        {
            dim3 dimBlock(1, kTileSize);
            dim3 dimGrid(1, ny_ / kTileSize);
            kernelRightBoundary<<<dimGrid, dimBlock>>>(H_, HU_, HV_, nx_, ny_, boundary_[BoundaryEdge::Right]);
        };

#ifdef DBG
        cout << "Call kernel to compute bottom/top boundaries " << flush << endl;
#endif
        switch (boundary_[BoundaryEdge::Bottom])
        {
        case BoundaryType::Connect:
        {
            // set ghost layer data in auxiliary data structure for ghost layer:
            for (int i = 0; i <= nx_ + 1; i++)
            {
                bottomGhostLayer->h[i]  = neighbour_[BoundaryEdge::Bottom]->h[i];
                bottomGhostLayer->hu[i] = neighbour_[BoundaryEdge::Bottom]->hu[i];
                bottomGhostLayer->hv[i] = neighbour_[BoundaryEdge::Bottom]->hv[i];
            };
        }
        case BoundaryType::Passive: /* also executed for case BoundaryType::Connect */
        {
            // copy ghost layer data from buffer bottomLayerDevice
            // into bottom ghost layer of unknowns
            dim3 dimBlock(kTileSize, 1);
            dim3 dimGrid(nx_ / kTileSize, 1);
            kernelBottomGhostBoundary<<<dimGrid, dimBlock>>>(H_, HU_, HV_, bottomLayerDevice, nx_, ny_);
            break;
        }
        default:
        {
            // set simple boundary_ conditions (BoundaryType::Outflow, BoundaryType::Wall) by resp. kernel:
            dim3 dimBlock(kTileSize, 1);
            dim3 dimGrid(nx_ / kTileSize, 1);
            kernelBottomBoundary<<<dimGrid, dimBlock>>>(H_, HU_, HV_, nx_, ny_, boundary_[BoundaryEdge::Bottom]);
        }
        };

        switch (boundary_[BoundaryEdge::Top])
        {
        case BoundaryType::Connect:
        {
            // set ghost layer data in auxiliary data structure for ghost layer:
            for (int i = 0; i <= nx_ + 1; i++)
            {
                topGhostLayer->h[i]  = neighbour_[BoundaryEdge::Top]->h[i];
                topGhostLayer->hu[i] = neighbour_[BoundaryEdge::Top]->hu[i];
                topGhostLayer->hv[i] = neighbour_[BoundaryEdge::Top]->hv[i];
            };
        }
        case BoundaryType::Passive: /* also executed for case BoundaryType::Connect */
        {
            // copy ghost layer data from buffer bottomLayerDevice
            // into bottom ghost layer of unknowns
            dim3 dimBlock(kTileSize, 1);
            dim3 dimGrid(nx_ / kTileSize, 1);
            kernelTopGhostBoundary<<<dimGrid, dimBlock>>>(H_, HU_, HV_, topLayerDevice, nx_, ny_);
            break;
        }
        default:
        {
            // set simple boundary_ conditions (BoundaryType::Outflow, BoundaryType::Wall) by resp. kernel:
            dim3 dimBlock(kTileSize, 1);
            dim3 dimGrid(nx_ / kTileSize, 1);
            kernelTopBoundary<<<dimGrid, dimBlock>>>(H_, HU_, HV_, nx_, ny_, boundary_[BoundaryEdge::Top]);
        }
        };
    }

    /**
     * synchronise the ghost layer content of h_, hu_, and hv_ in main memory
     * with device memory and auxiliary data structures, i.e. transfer
     * memory from main/auxiliary memory into device memory
     */
    void Block::synchGhostLayerAfterWrite()
    {

        if (boundary_[BoundaryEdge::Left] == BoundaryType::Passive
            || boundary_[BoundaryEdge::Left] == BoundaryType::Connect)
        {
#ifdef DBG
            cout << "Set left passive boundary_" << endl;
#endif
            // transfer h_, hu_, and hv_ from left ghost layer to resp. device memory
            CUDA_MEM_CPY(
                H_,
                h_[0],
                (ny_ + 2) * sizeof(RealType),
                cudaMemcpyHostToDevice,
                "Transferring left ghost layer to device"
            );
            CUDA_MEM_CPY(
                HU_,
                hu_[0],
                (ny_ + 2) * sizeof(RealType),
                cudaMemcpyHostToDevice,
                "Transferring left ghost layer to device"
            );
            CUDA_MEM_CPY(
                HV_,
                hv_[0],
                (ny_ + 2) * sizeof(RealType),
                cudaMemcpyHostToDevice,
                "Transferring left ghost layer to device"
            );
        };
        if (boundary_[BoundaryEdge::Right] == BoundaryType::Passive
            || boundary_[BoundaryEdge::Right] == BoundaryType::Connect)
        {
#ifdef DBG
            cout << "Set right passive boundary_" << endl;
#endif
            // transfer h_, hu_, and hv_ from right ghost layer to resp. device memory
            CUDA_MEM_CPY(
                H_ + ((nx_ + 1) * (ny_ + 2)),
                h_[nx_ + 1],
                (ny_ + 2) * sizeof(RealType),
                cudaMemcpyHostToDevice,
                "Transferring right ghost layer to device"
            );
            CUDA_MEM_CPY(
                HU_ + ((nx_ + 1) * (ny_ + 2)),
                hu_[nx_ + 1],
                (ny_ + 2) * sizeof(RealType),
                cudaMemcpyHostToDevice,
                "Transferring right ghost layer to device"
            );
            CUDA_MEM_CPY(
                HV_ + ((nx_ + 1) * (ny_ + 2)),
                hv_[nx_ + 1],
                (ny_ + 2) * sizeof(RealType),
                cudaMemcpyHostToDevice,
                "Transferring right ghost layer to device"
            );
        };

        // bottom and top boundary_ ghost layers are replicated (for efficiency reasons)
        // in the  memory regions starting from bottomLayer_ and topLayer_
        // -> these need to be transfered to device memory
        if (boundary_[BoundaryEdge::Bottom] == BoundaryType::Passive
            || boundary_[BoundaryEdge::Bottom] == BoundaryType::Connect)
        {
#ifdef DBG
            cout << "Set bottom passive boundary_" << endl;
            cout << " * transfer " << 3 * (nx_ + 2) * sizeof(RealType) << " bytes." << endl;
#endif
            // transfer bottom ghost layer
            // (3 arrays - h_, hu_, hv_ - of nx_+2 floats, consecutive in memory)
            CUDA_MEM_CPY(
                bottomLayerDevice,
                bottomLayer_,
                3 * (nx_ + 2) * sizeof(RealType),
                cudaMemcpyHostToDevice,
                "Transferring bottom ghost layer to device"
            );
        };

        if (boundary_[BoundaryEdge::Top] == BoundaryType::Passive
            || boundary_[BoundaryEdge::Top] == BoundaryType::Connect)
        {
#ifdef DBG
            cout << "Set top passive boundary_" << endl;
            cout << " * transfer " << 3 * (nx_ + 2) * sizeof(RealType) << " bytes." << endl;
#endif
            // transfer top ghost layer
            // (3 arrays - h_, hu_, hv_ - of nx_+2 floats, consecutive in memory)
            CUDA_MEM_CPY(
                topLayerDevice,
                topLayer_,
                3 * (nx_ + 2) * sizeof(RealType),
                cudaMemcpyHostToDevice,
                "Transferring top ghost layer to device"
            );
        };
    }

    /**
     * Update (for heterogeneous computing) variables h_, hu_, hv_ in copy layers
     * before an external access to the unknowns
     * (only required for BoundaryType::Passive and BoundaryType::Connect boundaries)
     * - copy (up-to-date) content from device memory into main memory
     */
    void Block::synchCopyLayerBeforeRead()
    {
#ifdef DBG
        cout << "Transfer copy layers from device to main memory" << flush << endl;
#endif

#if !defined(ENABLE_CUDA_UMA_ALLOCATOR)
        // copy values in copy layers to main memory
        // left copy layer:
        int offset = ny_ + 2;
        CUDA_MEM_CPY(
            h_[1],
            H_ + offset,
            (ny_ + 2) * sizeof(RealType),
            cudaMemcpyDeviceToHost,
            "Transferring left copy layer to device"
        );
        CUDA_MEM_CPY(
            hu_[1],
            HU_ + offset,
            (ny_ + 2) * sizeof(RealType),
            cudaMemcpyDeviceToHost,
            "Transferring left copy layer to device"
        );
        CUDA_MEM_CPY(
            hv_[1],
            HV_ + offset,
            (ny_ + 2) * sizeof(RealType),
            cudaMemcpyDeviceToHost,
            "Transferring left copy layer to device"
        );
        // right copy layer
        offset = nx_ * (ny_ + 2);
        CUDA_MEM_CPY(
            h_[nx_],
            H_ + offset,
            (ny_ + 2) * sizeof(RealType),
            cudaMemcpyDeviceToHost,
            "Transferring right copy layer to device"
        );
        CUDA_MEM_CPY(
            hu_[nx_],
            HU_ + offset,
            (ny_ + 2) * sizeof(RealType),
            cudaMemcpyDeviceToHost,
            "Transferring right copy layer to device"
        );
        CUDA_MEM_CPY(
            hv_[nx_],
            HV_ + offset,
            (ny_ + 2) * sizeof(RealType),
            cudaMemcpyDeviceToHost,
            "Transferring right copy layer to device"
        );

        int size = 3 * (nx_ + 2);
        // bottom copy layer
        if (boundary_[BoundaryEdge::Bottom] == BoundaryType::Passive
            || boundary_[BoundaryEdge::Bottom] == BoundaryType::Connect)
        {
            dim3 dimBlock(kTileSize, 1);
            dim3 dimGrid(nx_ / kTileSize, 1);
            kernelBottomCopyLayer<<<dimGrid, dimBlock>>>(H_, HU_, HV_, bottomLayerDevice + size, nx_, ny_);

            CUDA_MEM_CPY(
                bottomLayer_ + size,
                bottomLayerDevice + size,
                size * sizeof(RealType),
                cudaMemcpyDeviceToHost,
                "Transferring bottom copy layer to device"
            );
        };

        // top copy layer
        if (boundary_[BoundaryEdge::Top] == BoundaryType::Passive
            || boundary_[BoundaryEdge::Top] == BoundaryType::Connect)
        {
            dim3 dimBlock(kTileSize, 1);
            dim3 dimGrid(nx_ / kTileSize, 1);
            kernelTopCopyLayer<<<dimGrid, dimBlock>>>(H_, HU_, HV_, topLayerDevice + size, nx_, ny_);

            CUDA_MEM_CPY(
                topLayer_ + size,
                topLayerDevice + size,
                size * sizeof(RealType),
                cudaMemcpyDeviceToHost,
                "Transferring top copy layer to device"
            );
        };
#endif // ENABLE_CUDA_UMA_ALLOCATOR
    }

    /**
     * register the row or column layer next to a boundary_ as a "copy layer",
     * from which values will be copied into the ghost layer or a neighbour_;
     * @return	a Blocks::Block1D object that contains row variables h_, hu_, and hv_
     */
    Blocks::Block1D* Block::registerCopyLayer(BoundaryEdge edge)
    {

        // for TOP and BOTTOM layer, the implementation is identical to that in Blocks::Block
        // for LEFT and RIGHT layer, separate layers are used that avoid strided copies
        // when transferring memory between host and device memory
        switch (edge)
        {
        case BoundaryEdge::Left:
            return new Blocks::Block1D(h_.getColProxy(1), hu_.getColProxy(1), hv_.getColProxy(1));
        case BoundaryEdge::Right:
            return new Blocks::Block1D(h_.getColProxy(nx_), hu_.getColProxy(nx_), hv_.getColProxy(nx_));
        case BoundaryEdge::Bottom:
            // transfer bottom ghost and copy layer to extra Blocks::Block1D
            for (int i = 0; i <= nx_ + 1; i++)
            {
                bottomGhostLayer->h[i]  = h_[i][0];
                bottomGhostLayer->hu[i] = hu_[i][0];
                bottomGhostLayer->hv[i] = hv_[i][0];
                bottomCopyLayer->h[i]   = h_[i][1];
                bottomCopyLayer->hu[i]  = hu_[i][1];
                bottomCopyLayer->hv[i]  = hv_[i][1];
            };
            return bottomCopyLayer;
        case BoundaryEdge::Top:
            // transfer bottom ghost and copy layer to extra Blocks::Block1D
            for (int i = 0; i <= nx_ + 1; i++)
            {
                topGhostLayer->h[i]  = h_[i][ny_ + 1];
                topGhostLayer->hu[i] = hu_[i][ny_ + 1];
                topGhostLayer->hv[i] = hv_[i][ny_ + 1];
                topCopyLayer->h[i]   = h_[i][ny_];
                topCopyLayer->hu[i]  = hu_[i][ny_];
                topCopyLayer->hv[i]  = hv_[i][ny_];
            };
            return topCopyLayer;
        };
        return NULL;
    }

    /**
     * "grab" the ghost layer at the specific boundary_ in order to set boundary_ values
     * in this ghost layer externally.
     * The boundary_ conditions at the respective ghost layer is set to BoundaryType::Passive,
     * such that the grabbing program component is responsible to provide correct
     * values in the ghost layer, for example by receiving data from a remote
     * copy layer via MPI communication.
     * @param	specified edge
     * @return	a Blocks::Block1D object that contains row variables h_, hu_, and hv_
     */
    Blocks::Block1D* Block::grabGhostLayer(BoundaryEdge edge)
    {

        // for TOP and BOTTOM layer, the implementation is identical to that in Blocks::Block
        // for LEFT and RIGHT layer, separate layers are used that avoid strided copies
        // when transferring memory between host and device memory
        boundary_[edge] = BoundaryType::Passive;
        switch (edge)
        {
        case BoundaryEdge::Left:
            return new Blocks::Block1D(h_.getColProxy(0), hu_.getColProxy(0), hv_.getColProxy(0));
        case BoundaryEdge::Right:
            return new Blocks::Block1D(h_.getColProxy(nx_ + 1), hu_.getColProxy(nx_ + 1), hv_.getColProxy(nx_ + 1));
        case BoundaryEdge::Bottom:
            return bottomGhostLayer;
        case BoundaryEdge::Top:
            return topGhostLayer;
        };
        return NULL;
    }

    /**
     * Print some available information about the CUDA devices.
     */
    void Block::printDeviceInformation()
    {
        Tools::Logger::logger.getDefaultOutputStream() << "Printing device information";

        //! id of the CUDA device.
        int l_deviceId;
        cudaGetDevice(&l_deviceId);

        //! total number of CUDA devices on this host.
        int l_deviceCount;
        cudaGetDeviceCount(&l_deviceCount);

        //! drive and runtime version
        int l_driverVersion, l_runtimeVersion;
        cudaDriverGetVersion(&l_driverVersion);
        cudaRuntimeGetVersion(&l_runtimeVersion);

        //! device properties
        cudaDeviceProp l_deviceProperty;
        cudaGetDeviceProperties(&l_deviceProperty, l_deviceId);

        // print information about the current device

        Tools::Logger::logger.getDefaultOutputStream(
        ) << "Current CUDA device (relative to host): "
          << l_deviceId << " ( " << l_deviceCount << " in total)" << std::endl;

        Tools::Logger::logger.getDefaultOutputStream(
        ) << "CUDA device properties: "
          << l_deviceProperty.name << " (name), " << l_driverVersion << "/" << l_runtimeVersion
          << " (driver/runtime version), " << l_deviceProperty.major << "." << l_deviceProperty.minor
          << " (compute capability)" << std::endl;
    }

    //==================================================================
    // protected member functions for memory model:
    // in case of temporary variables (especial in non-local memory, for
    // example on accelerators), the main variables h_, hu_, hv_, and b
    // are not necessarily updated after each time step.
    // The following methods are called to synchronise before or after
    // external read or write to the variables.
    //==================================================================

    /**
     * Update all temporary and non-local (for heterogeneous computing) variables
     * after an external update of the main variables h_, hu_, hv_, and b.
     */
    void Block::synchAfterWrite()
    {
        // update h_, hu_, hv_, b in device memory
        synchWaterHeightAfterWrite();
        synchDischargeAfterWrite();
        synchBathymetryAfterWrite();

        // update the auxiliary data structures for copy and ghost layers
        // at bottom (and top, see below) boundaries
        // -> only required for BoundaryType::Passive and BoundaryType::Connect boundaries
        if (boundary_[BoundaryEdge::Bottom] == BoundaryType::Passive
            || boundary_[BoundaryEdge::Bottom] == BoundaryType::Connect)
        {
            // transfer bottom ghost and copy layer to extra Blocks::Block1D
            for (int i = 0; i <= nx_ + 1; i++)
            {
                bottomGhostLayer->h[i]  = h_[i][0];
                bottomGhostLayer->hu[i] = hu_[i][0];
                bottomGhostLayer->hv[i] = hv_[i][0];
                bottomCopyLayer->h[i]   = h_[i][1];
                bottomCopyLayer->hu[i]  = hu_[i][1];
                bottomCopyLayer->hv[i]  = hv_[i][1];
            };
        };

        if (boundary_[BoundaryEdge::Top] == BoundaryType::Passive
            || boundary_[BoundaryEdge::Top] == BoundaryType::Connect)
        {
            // transfer bottom ghost and copy layer to extra Blocks::Block1D
            for (int i = 0; i <= nx_ + 1; i++)
            {
                topGhostLayer->h[i]  = h_[i][ny_ + 1];
                topGhostLayer->hu[i] = hu_[i][ny_ + 1];
                topGhostLayer->hv[i] = hv_[i][ny_ + 1];
                topCopyLayer->h[i]   = h_[i][ny_];
                topCopyLayer->hu[i]  = hu_[i][ny_];
                topCopyLayer->hv[i]  = hv_[i][ny_];
            };
        };
    }

    /**
     * Update temporary and non-local (for heterogeneous computing) variables
     * after an external update of the water height h_
     */
    void Block::synchWaterHeightAfterWrite()
    {
#ifdef DBG
        cout << "Load water height h_ into device memory" << flush << endl;
#endif
#if !defined(ENABLE_CUDA_UMA_ALLOCATOR)
        CUDA_MEM_CPY(H_, h_.getData(), full_grid_size_in_bytes_, cudaMemcpyHostToDevice); // TODO review th, "memory of
#endif                                                                                    // h_ not transferred"is
    }

    /**
     * Update temporary and non-local (for heterogeneous computing) variables
     * after an external update of the discharge variables hu_ and hv_
     */
    void Block::synchDischargeAfterWrite()
    {
#ifdef DBG
        cout << "Load discharge hu_ and hv_ into device memory" << flush << endl;
#endif
#if !defined(ENABLE_CUDA_UMA_ALLOCATOR)
        CUDA_MEM_CPY(HU_, hu_.getData(), full_grid_size_in_bytes_, cudaMemcpyHostToDevice); // TODO review th, "memory
                                                                                            // of hu_ not transferred"is
        CUDA_MEM_CPY(HV_, hv_.getData(), full_grid_size_in_bytes_, cudaMemcpyHostToDevice); // TODO review th, "memory
                                                                                            // of hv_ not transferred"is
#endif
    }

    /**
     * Update temporary and non-local (for heterogeneous computing) variables
     * after an external update of the bathymetry b
     */
    void Block::synchBathymetryAfterWrite()
    {
#ifdef DBG
        cout << "Load bathymetry unknowns into device memory" << flush << endl;
#endif
#if !defined(ENABLE_CUDA_UMA_ALLOCATOR)
        CUDA_MEM_CPY(B_, b_.getData(), full_grid_size_in_bytes_, cudaMemcpyHostToDevice); // TODO review th, "memory of
                                                                                          // b not transferred"is
#endif
        //  computeBathymetrySources();
    }

    /**
     * Update the main variables h_, hu_, hv_, and b before an external read access:
     * copies current content of the respective device variables H_, HU_, HV_, B_
     */
    void Block::synchBeforeRead()
    {
        synchWaterHeightBeforeRead();
        synchDischargeHuBeforeRead();
        synchDischargeHvBeforeRead();
        synchBathymetryBeforeRead();

        /* --- the following part is excluded:
           --- it should not be necessary to update the auxiliary data structures
           --- for top/bottom copy/ghost layers in main memory: these need to be
           --- kept consistent by class Block (in particular, they cannot be
           --- changed externally).
           --- */
        //   if (boundary_[BoundaryEdge::Bottom] == BoundaryType::Passive || boundary_[BoundaryEdge::Bottom] ==
        //   BoundaryType::Connect) {
        //      // transfer bottom ghost and copy layer to extra Blocks::Block1D
        //      for(int i=0; i<=nx_+1; i++) {
        //        h_[i][0] = bottomGhostLayer->h[i];
        //        hu_[i][0]= bottomGhostLayer->hu[i];
        //        hv_[i][0]= bottomGhostLayer->hv[i];
        //        h_[i][1] = bottomCopyLayer->h[i];
        //        hu_[i][1]= bottomCopyLayer->hu[i];
        //        hv_[i][1]= bottomCopyLayer->hv[i];
        //      };
        //   };
        //
        //   if (boundary_[BoundaryEdge::Top] == BoundaryType::Passive || boundary_[BoundaryEdge::Top] ==
        //   BoundaryType::Connect) {
        //      // transfer bottom ghost and copy layer to extra Blocks::Block1D
        //      for(int i=0; i<=nx_+1; i++) {
        //        h_[i][ny_+1]  = topGhostLayer->h[i];
        //        hu_[i][ny_+1] = topGhostLayer->hu[i];
        //        hv_[i][ny_+1] = topGhostLayer->hv[i];
        //        h_[i][ny_]  = topCopyLayer->h[i];
        //        hu_[i][ny_] = topCopyLayer->hu[i];
        //        hv_[i][ny_] = topCopyLayer->hv[i];
        //      };
        //   };
    }

    /**
     * Update temporary and non-local (for heterogeneous computing) variables
     * before an external access to the water height h_
     */
    void Block::synchWaterHeightBeforeRead()
    {
#ifdef DBG
        cout << "Copy water height h_ from device" << flush << endl;
#endif
#if !defined(ENABLE_CUDA_UMA_ALLOCATOR)
        CUDA_MEM_CPY(h_.getData(), H_, full_grid_size_in_bytes_, cudaMemcpyDeviceToHost); // TODO review th, "memory of
#endif                                                                                    // h_ not transferred"is

#ifdef DBG
        cout << "Set water height in ghost-layer corner cells" << flush << endl;
#endif
        // only required for visualisation: set values in corner ghost cells
        h_[0][0]             = h_[1][1];
        h_[0][ny_ + 1]       = h_[1][ny_];
        h_[nx_ + 1][0]       = h_[nx_][1];
        h_[nx_ + 1][ny_ + 1] = h_[nx_][ny_];
    }

    /**
     * Update temporary and non-local (for heterogeneous computing) variables
     * before an external access to the discharge variables hu_
     */
    void Block::synchDischargeHuBeforeRead()
    {
#ifdef DBG
        cout << "Copy discharge hu_ and hv_ from device" << flush << endl;
#endif
#if !defined(ENABLE_CUDA_UMA_ALLOCATOR)
        CUDA_MEM_CPY(hu_.getData(), HU_, full_grid_size_in_bytes_, cudaMemcpyDeviceToHost);
#endif
    }

    /**
     * Update temporary and non-local (for heterogeneous computing) variables
     * before an external access to the discharge variables hv_
     */
    void Block::synchDischargeHvBeforeRead()
    {
#ifdef DBG
        cout << "Copy discharge hu_ and hv_ from device" << flush << endl;
#endif
#if !defined(ENABLE_CUDA_UMA_ALLOCATOR)
        CUDA_MEM_CPY(hv_.getData(), HV_, full_grid_size_in_bytes_, cudaMemcpyDeviceToHost);
#endif
    }

    /**
     * Update temporary and non-local (for heterogeneous computing) variables
     * before an external access to the bathymetry b
     */
    void Block::synchBathymetryBeforeRead()
    {
#ifdef DBG
        cout << "Copy water bathymetry b from device" << flush << endl;
#endif
#if !defined(ENABLE_CUDA_UMA_ALLOCATOR)
        CUDA_MEM_CPY(b_.getData(), B_, full_grid_size_in_bytes_, cudaMemcpyDeviceToHost); // TODO review th, "memory of
#endif                                                                                    // b not transferred"is
    }

    void Block::init(int i_cudaDevice)
    {
        std::cout << std::endl;

        Tools::Logger::logger.setProcessRank(i_cudaDevice);

// REVIEW[epic=SWE,seq=55] enable cude device setting?
// cudaSetDevice(i_cudaDevice);

// check for a valid CUDA device id
#ifndef NDEBUG
        int l_deviceCount;
        cudaGetDeviceCount(&l_deviceCount);
        assert((i_cudaDevice >= 0) && (i_cudaDevice < l_deviceCount));
#endif

        printDeviceInformation();

        // Make sure the cuda device is reset at exit
        atexit(Block::finalize);
    }

    void Block::finalize()
    {
        // reset the cuda device
        Tools::Logger::logger.getDefaultOutputStream() << "Resetting the CUDA devices";
        cudaDeviceReset();
    }

} // namespace cuda