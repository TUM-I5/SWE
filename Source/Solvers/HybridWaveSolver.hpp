/**
 * Hybrid.hpp
 *
 ****
 **** Hybrid Riemann Solver for the Shallow Water Equation
 ****
 *
 *  Created on: Jan 04, 2012
 *  Last Update: Feb 20, 2012
 *
 ****
 *
 *  Author: Alexander Breuer
 *    Homepage: http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer
 *    E-Mail: breuera AT in.tum.de
 *
 ****
 *
 * (Main) Literature:
 *
 *   @webpage{levequeclawpack,
 *            Author = {LeVeque, R. J.},
 *            Lastchecked = {January, 05, 2011},
 *            Title = {Clawpack Sofware},
 *            Url = {https://github.com/dlgeorge/clawpack-4.x/blob/master/geoclaw/2d/lib/rpn2ez_fast_geo.f}}
 *
 *   TODO: Where to find literature about hybrid-solvers?
 *
 ****
 *
 * Acknowledgments:
 *   Special thanks go to R.J. LeVeque and D.L. George for publishing their code
 *   and the corresponding documentation (-> Literature).
 */

#ifndef HYBRID_HPP_
#define HYBRID_HPP_

#include "AugRie.hpp"
#include "FWave.hpp"

#ifndef NDEBUG
#include <iostream>
#endif

namespace solver {
  template <typename T>
  class Hybrid;
} // namespace solver

/**
 * Hybrid Wave Propagation Solver.
 * Combines the "simple" f-Wave approach with the more complex Augmented Riemann Solver.
 *
 * T should be double or float.
 */
template <typename T>
class solver::Hybrid {
private:
  //! gravity constant.
  const T g;
  //! numerical definition of "dry".
  const T dryTol;
  //! f-Wave solver for "simple" problems.
  solver::FWave<T> fWaveSolver;
  //! Augmented Riemann-solver for more complex problems.
  solver::AugRie<T> augRieSolver;

#ifndef NDEBUG
  //! how often was the f-Wave solver used.
  long counterFWave;

  //! how often was the AugRie solver used.
  long counterAugRie;
#endif

public:
  /**
   * Constructor of the hybrid solver with optional parameters.
   *
   * @param i_dryTolerance numerical definition of "dry".
   * @param i_gravity gravity constant.
   * @param i_newtonTolerance numerical definition of "convergence" (used in the AugRie solver only).
   * @param i_maxNumberOfNewtonIterations maximum steps for the Newton-Raphson method (used in the AugRie solver only).
   * @param i_zeroTolerance numerical definition of zero.
   */
  Hybrid(
    T   i_dryTolerance                = (T)0.01,
    T   i_gravity                     = (T)9.81,
    T   i_newtonTolerance             = (T)0.000001,
    int i_maxNumberOfNewtonIterations = 10,
    T   zeroTolerance                 = (T)0.000000001
  ):
    g(i_gravity),
    dryTol(i_dryTolerance),
    fWaveSolver(i_dryTolerance, i_gravity, zeroTolerance),
    augRieSolver(i_dryTolerance, i_gravity, i_newtonTolerance, i_maxNumberOfNewtonIterations, zeroTolerance) {
#ifndef NDEBUG
    // set the counters to zero.
    counterFWave = counterAugRie = 0;

#ifndef SUPPRESS_SOLVER_DEBUG_OUTPUT
    // print some information about the used solver.
    std::cout
      << "  *** solver::Hybrid created" << std::endl
      << "    zeroTolerance=" << zeroTolerance << std::endl
      << "    gravity=" << i_gravity << std::endl
      << "    dryTolerance=" << i_dryTolerance << std::endl
      << "    newtonTolerance=" << i_newtonTolerance << std::endl
      << "    maxNumberNumberOfNewtonIterations=" << i_maxNumberOfNewtonIterations << std::endl
      << "  ***\n\n";
#endif
#endif
  };

  /**
   * Compute net updates for the cell on the left/right side of the edge.
   * The "correct" solver is used automatically.
   *
   * @param i_hLeft height on the left side of the edge.
   * @param i_hRight height on the right side of the edge.
   * @param i_huLeft momentum on the left side of the edge.
   * @param i_huRight momentum on the right side of the edge.
   * @param i_bLeft bathymetry on the left side of the edge.
   * @param i_bRight bathymetry on the right side of the edge.
   *
   * @param o_hUpdateLeft will be set to: Net-update for the height of the cell on the left side of the edge.
   * @param o_hUpdateRight will be set to: Net-update for the height of the cell on the right side of the edge.
   * @param o_huUpdateLeft will be set to: Net-update for the momentum of the cell on the left side of the edge.
   * @param o_huUpdateRight will be set to: Net-update for the momentum of the cell on the right side of the edge.
   * @param o_maxWaveSpeed will be set to: Maximum (linearized) wave speed -> Should be used in the CFL-condition.
   */
  void computeNetUpdates(
    const T& i_hLeft,
    const T& i_hRight,
    const T& i_huLeft,
    const T& i_huRight,
    const T& i_bLeft,
    const T& i_bRight,

    T& o_hUpdateLeft,
    T& o_hUpdateRight,
    T& o_huUpdateLeft,
    T& o_huUpdateRight,
    T& o_maxWaveSpeed
  ) {
    if (useAugmentedRiemannSolver(i_hLeft, i_hRight, i_huLeft, i_huRight, i_bLeft, i_bRight, o_hUpdateLeft, o_hUpdateRight, o_huUpdateLeft, o_huUpdateRight, o_maxWaveSpeed) == true) {

      augRieSolver.computeNetUpdates(
        i_hLeft,
        i_hRight,
        i_huLeft,
        i_huRight,
        i_bLeft,
        i_bRight,

        o_hUpdateLeft,
        o_hUpdateRight,
        o_huUpdateLeft,
        o_huUpdateRight,
        o_maxWaveSpeed
      );

#ifndef NDEBUG
      counterAugRie++;
#endif

    }
#ifndef NDEBUG
    else {
      counterFWave++;
    }
#endif
  }

  /**
   * Should the Augmented Riemann Solver be used?
   *  Side effect: If its sufficient to use the f-Wave solver the net-updates
   *               will be computed already (hybrid f-Wave version).
   *
   * @param i_hLeft height on the left side of the edge.
   * @param i_hRight height on the right side of the edge.
   * @param i_huLeft momentum on the left side of the edge.
   * @param i_huRight momentum on the right side of the edge.
   * @param i_bLeft bathymetry on the left side of the edge.
   * @param i_bRight bathymetry on the right side of the edge.
   *
   * @param o_hUpdateLeft will be set to: Net-update for the height of the cell on the left side of the edge.
   * @param o_hUpdateRight will be set to: Net-update for the height of the cell on the right side of the edge.
   * @param o_huUpdateLeft will be set to: Net-update for the momentum of the cell on the left side of the edge.
   * @param o_huUpdateRight will be set to: Net-update for the momentum of the cell on the right side of the edge.
   * @param o_maxWaveSpeed will be set to: Maximum (linearized) wave speed -> Should be used in the CFL-condition.
   *
   * @return A boolean whether the Augmented Riemann Solver should be used or not.
   */
  bool useAugmentedRiemannSolver(
    const T& i_hLeft,
    const T& i_hRight,
    const T& i_huLeft,
    const T& i_huRight,
    const T& i_bLeft,
    const T& i_bRight,

    T& o_hUpdateLeft,
    T& o_hUpdateRight,
    T& o_huUpdateLeft,
    T& o_huUpdateRight,
    T& o_maxWaveSpeed
  ) {
    // use the augmented solver by default
    bool useAugRieSolver = true;

    if (i_hLeft > dryTol && i_hRight > dryTol) { // both cells wet?

      // compute the velocities
      T uLeft  = i_huLeft / i_hLeft;
      T uRight = i_huRight / i_hRight;

      // definition of a "simple" Riemann-problem (GeoClaw)
      // TODO: Literature?
      if( std::fabs(i_hLeft - i_hRight) < (T)0.001*std::min(i_hLeft, i_hRight) && //"small" jump in height
            std::max(uLeft*uLeft, uRight*uRight)/std::max(g*i_hLeft, g*i_hRight) < (T)0.01 ) { //"small" jump in velocity
        useAugRieSolver = false;

        // compute simplified wave speeds (Geoclaw)
        // TODO: Literature?
        T waveSpeeds[2];
        waveSpeeds[0] = (uLeft + uRight) * (T).5 - std::sqrt(g * (i_hLeft + i_hRight) * (T).5);
        waveSpeeds[1] = (uLeft + uRight) * (T).5 + std::sqrt(g * (i_hLeft + i_hRight) * (T).5);

        // compute the f-waves
        fWaveSolver.computeNetUpdatesHybrid(
          i_hLeft,
          i_hRight,
          i_huLeft,
          i_huRight,
          i_bLeft,
          i_bRight,
          uLeft,
          uRight,
          waveSpeeds,

          o_hUpdateLeft,
          o_hUpdateRight,
          o_huUpdateLeft,
          o_huUpdateRight,
          o_maxWaveSpeed
        );
      }
    }
    // DryDry case: Don't use any solver and set everything to zero
    else if (i_hLeft < dryTol && i_hRight < dryTol) {
      o_hUpdateLeft = o_hUpdateRight = o_huUpdateLeft = o_huUpdateRight = o_maxWaveSpeed = 0;
      useAugRieSolver                                                                    = false;
    }

    return useAugRieSolver;
  }

#ifndef NDEBUG
  /**
   * Prints the current statistics of the hybrid solver.
   */
  void printStats() {
    std::cout
      << "  *** solver::Hybrid printing statistics" << std::endl
      << "    times the f-Wave solver was used: " << counterFWave << std::endl
      << "    times the AugRie solver was used: " << counterAugRie << std::endl
      << "  ***\n\n";
  }

  /**
   * Get the current statistics of the hybrid solver.
   *
   * @param o_counterFWave will be set to: times the f-Wave solver was used.
   * @param o_counterAugRie will be set to: times the Augmented Riemann solver was used.
   */
  void getStats(long& o_counterFWave, long& o_counterAugRie) {
    o_counterFWave  = counterFWave;
    o_counterAugRie = counterAugRie;
  }
#endif
};

#endif /* HYBRID_HPP_ */
