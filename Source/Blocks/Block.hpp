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

#pragma once

#include "Scenarios/Scenario.hpp"
#include "Tools/Float2D.hpp"
#include "Tools/RealType.hpp"
#include "Types/BoundaryEdge.hpp"
#include "Types/BoundaryType.hpp"

namespace Blocks
{
    /**
     * Blocks::Block1D is a simple struct that can represent a single line or row of
     * Blocks::Block unknowns (using the Tools::Float1D proxy class).
     * It is intended to unify the implementation of inflow and periodic boundary
     * conditions, as well as the ghost/copy-layer connection between several Blocks::Block
     * grids.
     */
    struct Block1D
    {
        Block1D(
            const Tools::Float1D<RealType>& h_, const Tools::Float1D<RealType>& hu_, const Tools::Float1D<RealType>& hv_
        ):
            h(h_),
            hu(hu_),
            hv(hv_){};
        Block1D(RealType* h_, RealType* hu_, RealType* hv_, int size_, int stride_ = 1):
            h(h_, size_, stride_),
            hu(hu_, size_, stride_),
            hv(hv_, size_, stride_){};

        Tools::Float1D<RealType> h;
        Tools::Float1D<RealType> hu;
        Tools::Float1D<RealType> hv;
    };

    /**
     * Blocks::Block is the main data structure to compute our shallow water model
     * on a single Cartesian grid block:
     * Blocks::Block is an abstract class (and interface) that should be extended
     * by respective implementation classes.
     *
     * <h3>Cartesian Grid for Discretization:</h3>
     *
     * Blocks::Block uses a regular Cartesian grid of size #nx by #ny, where each
     * grid cell carries three unknowns:
     * - the water level #h
     * - the momentum components #hu and #hv (in x- and y- direction, resp.)
     * - the bathymetry #b
     *
<<<<<<< master
     * Each of the components is stored as a 2D array, implemented as a Tools::Float2D object,
     * and are defined on grid indices [0,..,#nx+1]*[0,..,#ny+1].
     * The computational domain is indexed with [1,..,#nx]*[1,..,#ny].
=======
     */
    Block(int nx, int ny, RealType dx, RealType dy);
    Block(
      int nx, int ny, RealType dx, RealType dy,
      Tools::Float2D<RealType>& h,
      Tools::Float2D<RealType>& hu,
      Tools::Float2D<RealType>& hv
    );
    
    /**
     * Sets the bathymetry on BoundaryType::Outflow or BoundaryType::Wall.
     * Should be called very time a boundary is changed to a BoundaryType::Outflow or
     * BoundaryType::Wall <b>or</b> the bathymetry changes.
     */
    void setBoundaryBathymetry();

    // Synchronization Methods
    /**
     * Updates all temporary and non-local (for heterogeneous computing) variables
     * after an external update of the main variables h, hu, hv, and b.
     */
    virtual void synchAfterWrite();
    virtual void synchWaterHeightAfterWrite();
    virtual void synchDischargeAfterWrite();
    virtual void synchBathymetryAfterWrite();

    /**
     * Updates the ghost layers (only for BoundaryType::Connect and BoundaryType::Passive conditions)
     * after an external update of the main variables h, hu, hv, and b in the
     * ghost layer.
     */
    virtual void synchGhostLayerAfterWrite();

    /**
     * Updates all temporary and non-local (for heterogeneous computing) variables
     * before an external access to the main variables h, hu, hv, and b.
     */
    virtual void synchBeforeRead();
    virtual void synchWaterHeightBeforeRead();
    virtual void synchDischargeBeforeRead();
    virtual void synchBathymetryBeforeRead();
    virtual void synchCopyLayerBeforeRead();

    /// Sets boundary conditions in ghost layers (set boundary conditions)
    /**
     * Sets the values of all ghost cells depending on the specifed
     * boundary conditions
     * - set boundary conditions for types BoundaryType::Wall and BoundaryType::Outflow
     * - derived classes need to transfer ghost layers
     */
    virtual void setBoundaryConditions();

  public:
    /**
     * Destructor: de-allocate all variables
     */
    virtual ~Block() = default;

    static Block* getBlockInstance(int nx, int ny, RealType dx, RealType dy);
    static Block* getBlockInstance(
      int nx, int ny, RealType dx, RealType dy,
      Tools::Float2D<RealType>& h,
      Tools::Float2D<RealType>& hu,
      Tools::Float2D<RealType>& hv
    );

    /// Initialises unknowns to a specific scenario
    /**
     * Initialises the unknowns and bathymetry in all grid cells according to the given Scenarios::Scenario.
>>>>>>> master
     *
     * The mesh sizes of the grid in x- and y-direction are stored in static variables
     * #dx and #dy. The position of the Cartesian grid in space is stored via the
     * coordinates of the left-bottom corner of the grid, in the variables
     * #offsetX and #offsetY.
     *
     * <h3>Ghost layers:</h3>
     *
     * To implement the behaviour of the fluid at boundaries and for using
     * multiple block in serial and parallel settings, Blocks::Block adds an
     * additional layer of so-called ghost cells to the Cartesian grid,
     * as illustrated in the following figure.
     * Cells in the ghost layer have indices 0 or #nx+1 / #ny+1.
     *
     * \image html ghost_cells.gif
     *
     * <h3>Memory Model:</h3>
     *
     * The variables #h, #hu, #hv for water height and momentum will typically be
     * updated by classes derived from Blocks::Block. However, it is not assumed that
     * such and updated will be performed in every time step.
     * Instead, subclasses are welcome to update #h, #hu, and #hv in a lazy fashion,
     * and keep data in faster memory (incl. local memory of acceleration hardware,
     * such as GPGPUs), instead.
     *
     * It is assumed that the bathymetry data #b is not changed during the algorithm
     * (up to the exceptions mentioned in the following).
     *
     * To force a synchronization of the respective data structures, the following
     * methods are provided as part of Blocks::Block:
     * - synchAfterWrite() to synchronize #h, #hu, #hv, and #b after an external update
     *   (reading a file, e.g.);
     * - synchWaterHeightAfterWrite(), synchDischargeAfterWrite(), synchBathymetryAfterWrite():
     *   to synchronize only #h or momentum (#hu and #hv) or bathymetry #b;
     * - synchGhostLayerAfterWrite() to synchronize only the ghost layers
     * - synchBeforeRead() to synchronize #h, #hu, #hv, and #b before an output of the
     *   variables (writing a visualization file, e.g.)
     * - synchWaterHeightBeforeRead(), synchDischargeHuBeforeRead(), synchDischargeHvBeforeRead(),
     * synchBathymetryBeforeRead(): as synchBeforeRead(), but only for the specified variables
     * - synchCopyLayerBeforeRead(): synchronizes the copy layer only (i.e., a layer that
     *   is to be replicated in a neighbouring Blocks::Block.
     *
     * <h3>Derived Classes</h3>
     *
     * As Blocks::Block just provides an abstract base class together with the most
     * important data structures, the implementation of concrete models is the
     * job of respective derived classes (see the class diagram at the top of this
     * page). Similar, parallel implementations that are based on a specific
     * parallel programming model (such as OpenMP) or parallel architecture
     * (such as GPU/CUDA) should form subclasses of their own.
     * Please refer to the documentation of these classes for more details on the
     * model and on the parallelisation approach.
     */
    class Block
    {
      protected:
        // Grid size: number of cells (incl. ghost layer in x and y direction):
        int nx_; ///< Size of Cartesian arrays in x-direction
        int ny_; ///< Size of Cartesian arrays in y-direction
        // Mesh size dx and dy:
        RealType dx_; ///< Mesh size of the Cartesian grid in x-direction
        RealType dy_; ///< Mesh size of the Cartesian grid in y-direction

        // Define arrays for unknowns:
        // h (water level) and u, v (velocity in x and y direction)
        // hd, ud, and vd are respective CUDA arrays on GPU
        Tools::Float2D<RealType> h_;  ///< Array that holds the water height for each element
        Tools::Float2D<RealType> hu_; ///< Array that holds the x-component of the momentum for each element (water
                                      ///< height h multiplied by velocity in x-direction)
        Tools::Float2D<RealType> hv_; ///< Array that holds the y-component of the momentum for each element (water
                                      ///< height h multiplied by velocity in y-direction)
        Tools::Float2D<RealType> b_; ///< Array that holds the bathymetry data (sea floor elevation) for each element

        /// Type of boundary conditions at Left, Right, Top, and Bottom boundary
        BoundaryType boundary_[4];
        /// For Connect boundaries: pointer to connected neighbour block
        const Block1D* neighbour_[4];

        /// Maximum time step allowed to ensure stability of the method
        /**
         * maxTimeStep_ can be updated as part of the methods computeNumericalFluxes
         * and updateUnknowns (depending on the numerical method)
         */
        RealType maxTimeStep_;

        // Offset of current block
        RealType offsetX_; ///< x-coordinate of the origin (left-bottom corner) of the Cartesian grid
        RealType offsetY_; ///< y-coordinate of the origin (left-bottom corner) of the Cartesian grid

        /**
         * Constructor: allocate variables for simulation
         *
         * unknowns h (water height), hu,hv (discharge in x- and y-direction),
         * and b (bathymetry) are defined on grid indices [0,..,nx+1]*[0,..,ny+1]
         * -> computational domain is [1,..,nx]*[1,..,ny]
         * -> plus ghost cell layer
         *
         * The constructor is protected: no instances of Blocks::Block can be
         * generated.
         *
         */
        Block(int nx, int ny, RealType dx, RealType dy);

        /**
         * Sets the bathymetry on BoundaryType::Outflow or BoundaryType::Wall.
         * Should be called very time a boundary is changed to a BoundaryType::Outflow or
         * BoundaryType::Wall <b>or</b> the bathymetry changes.
         */
        void setBoundaryBathymetry();

        // Synchronization Methods
        /**
         * Updates all temporary and non-local (for heterogeneous computing) variables
         * after an external update of the main variables h, hu, hv, and b.
         */
        virtual void synchAfterWrite();
        virtual void synchWaterHeightAfterWrite();
        virtual void synchDischargeAfterWrite();
        virtual void synchBathymetryAfterWrite();

        /**
         * Updates the ghost layers (only for BoundaryType::Connect and BoundaryType::Passive conditions)
         * after an external update of the main variables h, hu, hv, and b in the
         * ghost layer.
         */
        virtual void synchGhostLayerAfterWrite();

        /**
         * Updates all temporary and non-local (for heterogeneous computing) variables
         * before an external access to the main variables h, hu, hv, and b.
         */
        virtual void synchBeforeRead();
        virtual void synchWaterHeightBeforeRead();
        virtual void synchDischargeHuBeforeRead();
        virtual void synchDischargeHvBeforeRead();
        virtual void synchBathymetryBeforeRead();
        virtual void synchCopyLayerBeforeRead();

        /// Sets boundary conditions in ghost layers (set boundary conditions)
        /**
         * Sets the values of all ghost cells depending on the specifed
         * boundary conditions
         * - set boundary conditions for types BoundaryType::Wall and BoundaryType::Outflow
         * - derived classes need to transfer ghost layers
         */
        virtual void setBoundaryConditions();

      public:
        /**
         * Destructor: de-allocate all variables
         */
        virtual ~Block() = default;

        static Block* getBlockInstance(int nx, int ny, RealType dx, RealType dy);

        /// Initialises unknowns to a specific scenario
        /**
         * Initialises the unknowns and bathymetry in all grid cells according to the given Scenarios::Scenario.
         *
         * In the case of multiple Blocks::Block at this point, it is not clear how the boundary conditions
         * should be set. This is because an isolated Blocks::Block doesn't have any information about the grid.
         * Therefore the calling routine, which has the information about multiple blocks, has to take care about
         * setting the right boundary conditions.
         *
         * @param scenario Scenarios::Scenario, which is used during the setup.
         * @param useMultipleBlocks Are there multiple blocks?
         */
        void initialiseScenario(
            RealType offsetX, RealType offsetY, Scenarios::Scenario& scenario, const bool useMultipleBlocks = false
        );

        /// Sets the water height according to a given function
        /**
         * Sets water height h in all interior grid cells (i.e. except ghost layer)
         * to values specified by parameter function h.
         */
        void setWaterHeight(RealType (*h)(RealType, RealType));

        /// Sets the momentum/discharge according to the provided functions
        /**
         * Sets discharge in all interior grid cells (i.e. except ghost layer)
         * to values specified by parameter functions.
         * Note: unknowns hu and hv represent momentum, while parameters u and v are velocities!
         */
        void setDischarge(RealType (*u)(RealType, RealType), RealType (*v)(RealType, RealType));

        /// Sets the bathymetry to a uniform value
        /**
         * Sets Bathymetry b in all grid cells (incl. ghost/boundary layers)
         * to a uniform value;
         * bathymetry source terms are re-computed.
         */
        void setBathymetry(RealType b);

        /// Sets the bathymetry according to a given function
        /**
         * Sets Bathymetry b in all grid cells (incl. ghost/boundary layers)
         * using the specified bathymetry function;
         * bathymetry source terms are re-computed.
         */
        void setBathymetry(RealType (*b)(RealType, RealType));

        // Read access to arrays of unknowns
        const Tools::Float2D<RealType>& getWaterHeight();
        const Tools::Float2D<RealType>& getDischargeHu();
        const Tools::Float2D<RealType>& getDischargeHv();
        const Tools::Float2D<RealType>& getBathymetry();

        /**
         * Sets the boundary type for specific block boundary.
         *
         * @param edge Location of the edge relative to the Blocks::Block.
         * @param boundaryType Type of the boundary condition.
         * @param inflow Pointer to an Blocks::Block1D, which specifies the inflow (should be nullptr for
         * BoundaryType::Wall or BoundaryType::Outflow).
         */
        void setBoundaryType(BoundaryEdge edge, BoundaryType boundaryType, const Block1D* inflow = nullptr);

        /// Returns a pointer to proxy class to access the copy layer
        /**
         * Registers the row or column layer next to a boundary as a "copy layer",
         * from which values will be copied into the ghost layer or a neighbour;
         * @return a Blocks::Block1D object that contains row variables h, hu, and hv.
         */
        virtual Block1D* registerCopyLayer(BoundaryEdge edge);

        /**
         * "Grab" the ghost layer at the specific boundary in order to set boundary values
         * in this ghost layer externally.
         * The boundary conditions at the respective ghost layer is set to BoundaryType::Passive,
         * such that the grabbing program component is responsible to provide correct
         * values in the ghost layer, for example by receiving data from a remote
         * copy layer via MPI communication.
         * @param specified edge.
         * @return a Blocks::Block1D object that contains row variables h, hu, and hv.
         */
        virtual Block1D* grabGhostLayer(BoundaryEdge edge);

        /**
         * Sets the values of all ghost cells depending on the specifed
         * boundary conditions;
         * if the ghost layer replicates the variables of a remote Blocks::Block,
         * the values are copied.
         */
        void setGhostLayer();

        /**
         * Computes the largest allowed time step for the current grid block
         * (reference implementation) depending on the current values of
         * variables h, hu, and hv, and store this time step size in member
         * variable maxTimestep.
         *
         * @param dryTol Dry tolerance (dry cells do not affect the time step).
         * @param cfl CFL number of the used method.
         */
        void computeMaxTimeStep(const RealType dryTol = 0.1f, const RealType cfl = 0.4f);

        /// Returns maximum size of the time step to ensure stability of the method
        RealType getMaxTimeStep() const;

        /// Executes a single time step (with fixed time step size) of the simulation
        virtual void simulateTimeStep(RealType dt);

        /// Performs the simulation starting with simulation time tStart, until simulation time tEnd is reached
        /**
         * Implements the main simulation loop between two checkpoints;
         * Note: this implementation can only be used, if you only use a single Blocks::Block
         *       and only apply simple boundary conditions!
         *       In particular, Blocks::Block::simulate can not trigger calls to exchange values
         *       of copy and ghost layers between blocks!
         * @param tStart Time where the simulation is started
         * @param tEnd Time of the next checkpoint
         * @return Actual end time reached
         */
        virtual RealType simulate(RealType tStart, RealType tEnd);

        /// Computes the numerical fluxes for each edge of the Cartesian grid
        /**
         * The computation of fluxes strongly depends on the chosen numerical
         * method. Hence, this purely virtual function has to be implemented
         * in the respective derived classes.
         */
        virtual void computeNumericalFluxes() = 0;

        /// Computes the new values of the unknowns h, hu, and hv in all grid cells
        /**
         * Based on the numerical fluxes (computed by computeNumericalFluxes)
         * and the specified time step size dt, an Euler time step is executed.
         * As the computational fluxes will depend on the numerical method,
         * this purely virtual function has to be implemented separately for
         * each specific numerical model (and parallelisation approach).
         * @param dt Size of the time step
         */
        virtual void updateUnknowns(RealType dt) = 0;

        // Access methods to grid sizes
        /// Returns #nx, i.e. the grid size in x-direction
        int getNx() const;
        /// Returns #ny, i.e. the grid size in y-direction
        int getNy() const;
    };

} // namespace Blocks