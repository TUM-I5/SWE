/**
 * AugRieSolver.hpp
 *
 ****
 **** Approximate Augmented Riemann Solver for the Shallow Water Equations
 ****
 *
 *  Created on: Sep 12, 2011
 *  Last Update: Feb 18, 2012
 *
 ****
 *
 *  Author: Alexander Breuer
 *    Homepage: http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer
 *    E-Mail: breuera AT in.tum.de
 *
 *  Some optimizations: Martin Schreiber
 *    Homepage: http://www5.in.tum.de/wiki/index.php/Martin_Schreiber
 *    E-Mail: schreibm AT in.tum.de
 *
 ****
 *
 * (Main) Literature:
 *
 *   @phdthesis{george2006finite,
 *              Author = {George, D.L.},
 *              Title = {Finite volume methods and adaptive refinement for tsunami propagation and inundation},
 *              Year = {2006}}
 *
 *   @article{george2008augmented,
 *            Author = {George, D.L.},
 *            Journal = {Journal of Computational Physics},
 *            Number = {6},
 *            Pages = {3089--3113},
 *            Publisher = {Elsevier},
 *            Title = {Augmented Riemann solvers for the shallow water equations over variable topography with steady
 *states and inundation}, Volume = {227}, Year = {2008}}
 *
 *   @book{leveque2002finite,
 *         Author = {LeVeque, R. J.},
 *         Date-Added = {2011-09-13 14:09:31 +0000},
 *         Date-Modified = {2011-10-31 09:46:40 +0000},
 *         Publisher = {Cambridge University Press},
 *         Title = {Finite Volume Methods for Hyperbolic Problems},
 *         Volume = {31},
 *         Year = {2002}}
 *
 *   @webpage{levequeclawpack,
 *            Author = {LeVeque, R. J.},
 *            Lastchecked = {January, 05, 2011},
 *            Title = {Clawpack Sofware},
 *            Url = {https://github.com/clawpack/clawpack-4.x/blob/master/geoclaw/2d/lib}}
 *
 ****
 *
 * Acknowledgments:
 *   Special thanks go to R.J. LeVeque and D.L. George for publishing their code
 *   and the corresponding documentation (-> Literature).
 */

/*
 * TODO: store nLow/nHigh variables in array[2]
 */

#pragma once

#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "WavePropagationSolver.hpp"

/** Switch features of the solver on/off
 *
 *  The solver is not strictly positivity preserving with correctrarefactions.
 *
 *  Solver seems to fail in the "Single wave on a simple beach"-benchmark
 *  if correctrarefactions or complexsteadystatewave is used.
 *
 *  TODO: Further investigation is recommended.
 */
//#define CORRECT_RARE_FACTIONS
//#define COMPLEX_STEADY_STATE_WAVE

namespace Solvers {

  /**
   * Approximate Augmented Riemann Solver for the Shallow Water Equations.
   *
   * T should be double or float.
   */
  template <class T>
  class AugRieSolver: public WavePropagationSolver<T> {
  private:
    // Use nondependent names (template base class)
    using WavePropagationSolver<T>::dryTol_;
    using WavePropagationSolver<T>::gravity_;
    using WavePropagationSolver<T>::zeroTol_;

    using WavePropagationSolver<T>::hLeft_;
    using WavePropagationSolver<T>::hRight_;
    using WavePropagationSolver<T>::huLeft_;
    using WavePropagationSolver<T>::huRight_;
    using WavePropagationSolver<T>::bLeft_;
    using WavePropagationSolver<T>::bRight_;
    using WavePropagationSolver<T>::uLeft_;
    using WavePropagationSolver<T>::uRight_;

    using WavePropagationSolver<T>::wetDryState_;

    using WavePropagationSolver<T>::DryDry;
    using WavePropagationSolver<T>::WetWet;
    using WavePropagationSolver<T>::WetDryInundation;
    using WavePropagationSolver<T>::WetDryWall;
    using WavePropagationSolver<T>::WetDryWallInundation;
    using WavePropagationSolver<T>::DryWetInundation;
    using WavePropagationSolver<T>::DryWetWall;
    using WavePropagationSolver<T>::DryWetWallInundation;

    //! Newton-tolerance (exit Newton-Raphson-method, if we are close enough to the root)
    const T newtonTol_;
    //! Maximum number of Newton-Raphson-Iterations
    const int maxNumberOfNewtonIterations_;

    //! Height of our homogeneous Riemann-problem at middle state (computed by determineMiddleState)
    T hMiddle_;

    //! Shock or inner rarefaction speeds of our homogeneous Riemann-problem (computed by computeMiddleState)
    T middleStateSpeeds_[2];

    /**
     * The Riemann-struture of the homogeneous Riemann-problem.
     */
    enum RiemannStructure {
      DrySingleRarefaction,  /**< 1st wave family: contact discontinuity; 2nd wave family: rarefaction. */
      SingleRarefactionDry,  /**< 1st wave family: rarefaction; 2nd wave family: contact discontinuity. */
      ShockShock,            /**< 1st wave family: shock; 2nd wave family: shock. */
      ShockRarefaction,      /**< 1st wave family: shock; 2nd wave family: rarefaction. */
      RarefactionShock,      /**< 1st wave family: rarefaction; 2nd wave family: shock. */
      RarefactionRarefaction /**< 1st wave family: rarefaction; 2nd wave family: rarefaction. */
    };

    //! Riemann-structure of our homogeneous Riemann-problem (determined by determineRiemannStructure)
    RiemannStructure riemannStructure_;

    // Declare variables which are used over and over again
#if 0
    T sqrtGh_[2];
    T sqrth_[2];
#define sqrtGhLeft_ (sqrtGh_[0])
#define sqrtGhRight_ (sqrtGh_[1])

#define sqrthLeft_ (sqrth_[0])
#define sqrthRight_ (sqrth_[1])
#else
    T sqrtGhLeft_;
    T sqrtGhRight_;

    T sqrthLeft_;
    T sqrthRight_;
#endif
    T sqrtG_;

  public:
    /**
     * Constructor of the Augmented Riemann solver with optional parameters.
     *
     * @param dryTolerance numerical definition of "dry".
     * @param gravity gravity constant.
     * @param newtonTolerance numerical definition of "convergence" (used in the AugRie solver only).
     * @param maxNumberOfNewtonIterations maximum steps for the Newton-Raphson method (used in the AugRie solver only).
     * @param zeroTolerance numerical definition of zero.
     */
    AugRieSolver(
      T   dryTolerance                = T(0.01),
      T   gravity                     = T(9.81),
      T   newtonTolerance             = T(0.000001),
      int maxNumberOfNewtonIterations = 10,
      T   zeroTolerance               = T(0.00001)
    ):
      WavePropagationSolver<T>(dryTolerance, gravity, zeroTolerance),
      newtonTol_(newtonTolerance),
      maxNumberOfNewtonIterations_(maxNumberOfNewtonIterations) {}

    ~AugRieSolver() override = default;

    /**
     * Compute net updates for the left/right cell of the edge.
     *
     * maxWaveSpeed will be set to the maximum (linearized) wave speed => CFL
     */
    void computeNetUpdates(
      const T& hLeft,
      const T& hRight,
      const T& huLeft,
      const T& huRight,
      const T& bLeft,
      const T& bRight,

      T& o_hUpdateLeft,
      T& o_hUpdateRight,
      T& o_huUpdateLeft,
      T& o_huUpdateRight,
      T& o_maxWaveSpeed
#if ENABLE_AUGMENTED_RIEMANN_EIGEN_COEFFICIENTS
      ,
      T o_eigenCoefficients[3]
#endif
    ) {
      // Store parameters to member variables
      WavePropagationSolver<T>::storeParameters(hLeft, hRight, huLeft, huRight, bLeft, bRight);

      // Set speeds to zero (will be determined later)
      uLeft_ = uRight_ = T(0.0);

      // Reset net updates and the maximum wave speed
      o_hUpdateLeft = o_hUpdateRight = o_huUpdateLeft = o_huUpdateRight = T(0.0);
      o_maxWaveSpeed                                                    = T(0.0);

#if ENABLE_AUGMENTED_RIEMANN_EIGEN_COEFFICIENTS
      // Reset eigen coefficients
      o_eigenCoefficients[0] = o_eigenCoefficients[1] = o_eigenCoefficients[2] = 0;
#endif

      // Determine the wet/dry state and compute local variables correspondingly
      determineWetDryState();

      if (wetDryState_ == WavePropagationSolver<T>::WetDryState::DryDry) {
        return; // Nothing to do in a dry region
      }

      // Precompute some terms which are fixed during
      // the computation after some specific point
      sqrtG_      = std::sqrt(gravity_);
      sqrthLeft_  = std::sqrt(hLeft);
      sqrthRight_ = std::sqrt(hRight);

      sqrtGhLeft_  = sqrtG_ * sqrthLeft_;
      sqrtGhRight_ = sqrtG_ * sqrthRight_;

      // Where to store the three waves
      T fWaves[3][2];
      // and their speeds
      T waveSpeeds[3];

      // Compute the augmented decomposition
      // (thats the place where the computational work is done..)
      computeWaveDecomposition(
        fWaves,
        waveSpeeds
#if ENABLE_AUGMENTED_RIEMANN_EIGEN_COEFFICIENTS
        ,
        o_eigenCoefficients
#endif
      );

      // Compute the updates from the three propagating waves
      // A^-\delta Q = \sum{s[i]<0} \beta[i] * r[i] = A^-\delta Q = \sum{s[i]<0} Z^i
      // A^+\delta Q = \sum{s[i]>0} \beta[i] * r[i] = A^-\delta Q = \sum{s[i]<0} Z^i
      for (int waveNumber = 0; waveNumber < 3; waveNumber++) {
        if (waveSpeeds[waveNumber] < -zeroTol_) { // Left going
          o_hUpdateLeft += fWaves[waveNumber][0];
          o_huUpdateLeft += fWaves[waveNumber][1];
        }

        else if (waveSpeeds[waveNumber] > zeroTol_) { // Right going
          o_hUpdateRight += fWaves[waveNumber][0];
          o_huUpdateRight += fWaves[waveNumber][1];
        } else { // TODO: this case should not happen mathematically, but it does. Where is the bug? Machine accuracy
                 // only?
          o_hUpdateLeft += T(0.5) * fWaves[waveNumber][0];
          o_huUpdateLeft += T(0.5) * fWaves[waveNumber][1];

          o_hUpdateRight += T(0.5) * fWaves[waveNumber][0];
          o_huUpdateRight += T(0.5) * fWaves[waveNumber][1];
        }

        // No wave speeds => zero strength fWaves
        //        assert(std::fabs(fWaves[waveNumber][0]) < zeroTol_);
        //        assert(std::fabs(fWaves[waveNumber][1]) < zeroTol_);
      }

#ifndef NDEBUG
      if (wetDryState_ == WavePropagationSolver<T>::WetDryState::DryWetWall) {
        assert(std::fabs(o_hUpdateLeft) < zeroTol_ && std::fabs(o_huUpdateLeft) < zeroTol_);
      } else if (wetDryState == WavePropagationSolver<T>::WetDryState::WetDryWall) {
        assert(std::fabs(o_hUpdateRight) < zeroTol_ && std::fabs(o_huUpdateRight) < zeroTol_);
      }
#endif

      // Compute maximum wave speed (-> CFL-condition)
      waveSpeeds[0] = std::fabs(waveSpeeds[0]);
      waveSpeeds[1] = std::fabs(waveSpeeds[1]);
      waveSpeeds[2] = std::fabs(waveSpeeds[2]);

      o_maxWaveSpeed = std::max(waveSpeeds[0], waveSpeeds[1]);
      o_maxWaveSpeed = std::max(o_maxWaveSpeed, waveSpeeds[2]);
    }

  private:
    /**
     * Determine the wet/dry state and set member variables accordingly.
     */
    void determineWetDryState() override {
      // Compute speeds or set them to zero (dry cells)
      if (hLeft_ > dryTol_) {
        uLeft_ = huLeft_ / hLeft_;
      } else {
        bLeft_ += hLeft_;
        hLeft_ = huLeft_ = uLeft_ = 0;
      }

      if (hRight_ > dryTol_) {
        uRight_ = huRight_ / hRight_;
      } else {
        bRight_ += hRight_;
        hRight_ = huRight_ = uRight_ = 0;
      }

      // Test for simple wet/wet case since this is most probably the
      // most frequently executed branch
      if (hLeft_ >= dryTol_ && hRight_ >= dryTol_) {
        wetDryState_ = WavePropagationSolver<T>::WetDryState::WetWet;
      }

      // Check for the dry/dry-case
      else if (hLeft_ < dryTol_ && hRight_ < dryTol_) {
        wetDryState_ = WavePropagationSolver<T>::WetDryState::DryDry;
      }

      // We have a shoreline: one cell dry, one cell wet

      // Check for simple inundation problems
      //  (=> dry cell lies lower than the wet cell)
      else if (hLeft_ < dryTol_ && hRight_ + bRight_ > bLeft_) {
        wetDryState_ = WavePropagationSolver<T>::WetDryState::DryWetInundation;
      }

      else if (hRight_ < dryTol_ && hLeft_ + bLeft_ > bRight_) {
        wetDryState_ = WavePropagationSolver<T>::WetDryState::WetDryInundation;
      }

      // Dry cell lies higher than the wet cell
      else if (hLeft_ < dryTol_) {
        // Lets check if the momentum is able to overcome the difference in height
        //   => solve homogeneous Riemann-problem to determine the middle state height
        //      which would arise if there is a wall (wall-boundary-condition)
        //        \cite[ch. 6.8.2]{george2006finite})
        //        \cite[ch. 5.2]{george2008augmented}
        computeMiddleState(hRight_, hRight_, -uRight_, uRight_, -huRight_, huRight_, maxNumberOfNewtonIterations_);

        if (hMiddle_ + bRight_ > bLeft_) {
          // Momentum is large enough, continue with the original values
          //          bLeft = hMiddle + bRight;
          wetDryState_ = WavePropagationSolver<T>::WetDryState::DryWetWallInundation;
        } else {
          // Momentum is not large enough, use wall-boundary-values
          hLeft_  = hRight_;
          uLeft_  = -uRight_;
          huLeft_ = -huRight_;
          bLeft_ = bRight_ = T(0.0);
          wetDryState_     = WavePropagationSolver<T>::WetDryState::DryWetWall;
        }
      }

      else if (hRight_ < dryTol_) {
        // Lets check if the momentum is able to overcome the difference in height
        //   => solve homogeneous Riemann-problem to determine the middle state height
        //      which would arise if there is a wall (wall-boundary-condition)
        //        \cite[ch. 6.8.2]{george2006finite})
        //        \cite[ch. 5.2]{george2008augmented}
        computeMiddleState(hLeft_, hLeft_, uLeft_, -uLeft_, huLeft_, -huLeft_, maxNumberOfNewtonIterations_);

        if (hMiddle_ + bLeft_ > bRight_) {
          // Momentum is large enough, continue with the original values
          //          bRight = hMiddle + bLeft;
          wetDryState_ = WavePropagationSolver<T>::WetDryState::WetDryWallInundation;
        } else {
          hRight_  = hLeft_;
          uRight_  = -uLeft_;
          huRight_ = -huLeft_;
          bRight_ = bLeft_ = T(0.0);
          wetDryState_     = WavePropagationSolver<T>::WetDryState::WetDryWall;
        }
      }
      // Done with all cases
      else {
        assert(false);
      }

      // Limit the effect of the source term if there is a "wall"
      //\cite[end of ch. 5.2?]{george2008augmented}
      //\cite[rpn2ez_fast_geo.f][levequeclawpack]
      if (wetDryState_ == WavePropagationSolver<T>::WetDryState::DryWetWallInundation) {
        bLeft_ = hRight_ + bRight_;
      }

      else if (wetDryState_ == WavePropagationSolver<T>::WetDryState::WetDryWallInundation) {
        bRight_ = hLeft_ + bLeft_;
      }
    }

    /**
     * Compute the augmented wave decomposition.
     *
     * @param o_fWaves will be set to: Decomposition into f-Waves.
     * @param o_waveSpeeds will be set to: speeds of the linearized waves (eigenvalues).
     */
    void computeWaveDecomposition(
      T o_fWaves[3][2],
      T o_waveSpeeds[3]
#ifdef ENABLE_AUGMENTED_RIEMANN_EIGEN_COEFFICIENTS
      ,
      T o_eigenCoefficients[3]
#endif
    ) {
      // Compute eigenvalues of the jacobian matrices in states Q_{i-1} and Q_{i} (char. speeds)
      T characteristicSpeeds[2]{};
      characteristicSpeeds[0] = uLeft_ - sqrtGhLeft_;
      characteristicSpeeds[1] = uRight_ + sqrtGhRight_;

      // Compute "Roe speeds"
      T hRoe = T(0.5) * (hRight_ + hLeft_);
      T uRoe = uLeft_ * sqrthLeft_ + uRight_ * sqrthRight_;
      uRoe /= sqrthLeft_ + sqrthRight_;

      T roeSpeeds[2]{};
      // Optimization for dumb compilers
      T sqrtGhRoe  = std::sqrt(gravity_ * hRoe);
      roeSpeeds[0] = uRoe - sqrtGhRoe;
      roeSpeeds[1] = uRoe + sqrtGhRoe;

      // Compute the middle state of the homogeneous Riemann-Problem
      if (wetDryState_ != WavePropagationSolver<T>::WetDryState::WetDryWall && wetDryState_ != WavePropagationSolver<T>::WetDryState::DryWetWall) {
        // Case WDW and DWW was computed in
        // determineWetDryState already
        computeMiddleState(hLeft_, hRight_, uLeft_, uRight_, huLeft_, huRight_);
      }

      // Compute extended Eindfeldt speeds (Einfeldt speeds + middle state speeds)
      //   \cite[ch. 5.2]{george2008augmented}, \cite[ch. 6.8]{george2006finite}
      T extEinfeldtSpeeds[2] = {T(0), T(0)};
      if (wetDryState_ == WavePropagationSolver<T>::WetDryState::WetWet || wetDryState_ == WavePropagationSolver<T>::WetDryState::WetDryWall || wetDryState_ == WavePropagationSolver<T>::WetDryState::DryWetWall) {
        extEinfeldtSpeeds[0] = std::min(characteristicSpeeds[0], roeSpeeds[0]);
        extEinfeldtSpeeds[0] = std::min(extEinfeldtSpeeds[0], middleStateSpeeds_[1]);

        extEinfeldtSpeeds[1] = std::max(characteristicSpeeds[1], roeSpeeds[1]);
        extEinfeldtSpeeds[1] = std::max(extEinfeldtSpeeds[1], middleStateSpeeds_[0]);
      } else if (hLeft_ < dryTol_) { // Ignore undefined speeds
        extEinfeldtSpeeds[0] = std::min(roeSpeeds[0], middleStateSpeeds_[1]);
        extEinfeldtSpeeds[1] = std::max(characteristicSpeeds[1], roeSpeeds[1]);

        assert(middleStateSpeeds[0] < extEinfeldtSpeeds[1]);
      } else if (hRight_ < dryTol_) { // Ignore undefined speeds
        extEinfeldtSpeeds[0] = std::min(characteristicSpeeds[0], roeSpeeds[0]);
        extEinfeldtSpeeds[1] = std::max(roeSpeeds[1], middleStateSpeeds_[0]);

        assert(middleStateSpeeds_[1] > extEinfeldtSpeeds[0]);
      } else {
        assert(false);
      }

      // HLL middle state
      //   \cite[theorem 3.1]{george2006finite}, \cite[ch. 4.1]{george2008augmented}
      T hLLMiddleHeight = huLeft_ - huRight_ + extEinfeldtSpeeds[1] * hRight_ - extEinfeldtSpeeds[0] * hLeft_;
      hLLMiddleHeight /= extEinfeldtSpeeds[1] - extEinfeldtSpeeds[0];
      hLLMiddleHeight = std::max(hLLMiddleHeight, T(0.0));

      // Define eigenvalues
      T eigenValues[3]{};
      eigenValues[0] = extEinfeldtSpeeds[0];
      // eigenValues[1] --> corrector wave
      eigenValues[2] = extEinfeldtSpeeds[1];

      // Define eigenvectors
      T eigenVectors[3][3]{};

      // Set first and third eigenvector
      eigenVectors[0][0] = T(1);
      eigenVectors[0][2] = T(1);

      eigenVectors[1][0] = eigenValues[0];
      eigenVectors[1][2] = eigenValues[2];

      eigenVectors[2][0] = eigenValues[0] * eigenValues[0];
      eigenVectors[2][2] = eigenValues[2] * eigenValues[2];

      // Compute rarefaction corrector wave
      //   \cite[ch. 6.7.2]{george2006finite}, \cite[ch. 5.1]{george2008augmented}
      bool strongRarefaction = false;
#ifdef CORRECT_RARE_FACTIONS
      if ((riemannStructure_ == ShockRarefaction || riemannStructure_ == RarefactionShock || riemannStructure_ == RarefactionRarefaction) && hMiddle_ > dryTol_) { // Limit to ensure non-negative depth

        // TODO: GeoClaw, riemann_aug_JCP; No explicit boundaries for "strong rarefaction" in literature?
        T rareBarrier[2] = {0.5, 0.9};

        // Lower rare barrier in the case of a transsonic rarefaction
        if ((riemannStructure_ == RarefactionShock || riemannStructure_ == RarefactionRarefaction) && eigenValues[0] * middleStateSpeeds_[0] < T(0.0)) { // Transsonic rarefaction, first wave family
          rareBarrier[0] = T(0.2);
        } else if ((riemannStructure_ == ShockRarefaction || riemannStructure_ == RarefactionRarefaction) && eigenValues[2] * middleStateSpeeds_[1] < T(0.0)) { // Transsonic rarefaction, second wave family
          rareBarrier[0] = T(0.2);
        }

        T sqrtGhMiddle = std::sqrt(gravity_ * hMiddle_);
        // Determine the max. rarefaction size (distance between head and tail)
        T rareFactionSize[2]{};
        rareFactionSize[0] = T(3.0) * (sqrtGhLeft_ - sqrtGhMiddle);
        rareFactionSize[1] = T(3.0) * (sqrtGhRight_ - sqrtGhMiddle);

        T maxRareFactionSize = std::max(rareFactionSize[0], rareFactionSize[1]);

        // Set the eigenvalue of the corrector wave in the case of a "strong rarefaction"
        if (maxRareFactionSize > rareBarrier[0] * (eigenValues[2] - eigenValues[0]) && maxRareFactionSize < rareBarrier[1] * (eigenValues[2] - eigenValues[0])) {
          strongRarefaction = true;
          if (rareFactionSize[0] > rareFactionSize[1]) {
            eigenValues[1] = middleStateSpeeds_[0];
          } else {
            eigenValues[1] = middleStateSpeeds_[1];
          }
        }

        // TODO: implemented in clawpack, why?
        // if (hMiddle < std::min(hLeft, hRight) / 5.0) { // Middle state in an HLL solve
        //   strongRarefaction = false;
        // }
      }
#endif

      // Set 2nd eigenvector
      if (strongRarefaction == false) {
        eigenValues[1]     = T(0.5) * (eigenValues[0] + eigenValues[2]);
        eigenVectors[0][1] = 0.0;
        eigenVectors[1][1] = 0.0;
        eigenVectors[2][1] = 1.0;
      } else {
        eigenVectors[0][1] = T(1.0);
        eigenVectors[1][1] = eigenValues[1];
        eigenVectors[2][1] = eigenValues[1] * eigenValues[1];
      }

      // Compute the jump in state
      T rightHandSide[3]{};
      rightHandSide[0] = hRight_ - hLeft_;
      rightHandSide[1] = huRight_ - huLeft_;
      rightHandSide[2] = huRight_ * uRight_ + T(0.5) * gravity_ * hRight_ * hRight_
                         - (huLeft_ * uLeft_ + T(0.5) * gravity_ * hLeft_ * hLeft_);

      // Compute steady state wave
      //   \cite[ch. 4.2.1 \& app. A]{george2008augmented}, \cite[ch. 6.2 \& ch. 4.4]{george2006finite}
      T steadyStateWave[2]{};
      T hBar = (hLeft_ + hRight_) * T(0.5);
#ifdef COMPLEX_STEADY_STATE_WAVE
      T lLambdaBar = (uLeft_ + uRight_) * (uLeft_ + uRight_) * T(0.25) - gravity_ * hBar;

      T lLambdaTilde = std::max(T(0.0), uLeft_ * uRight_) - gravity_ * hBar;

      // Near sonic as defined in geoclaw (TODO: literature?)
      if ((std::fabs(lLambdaBar) < zeroTol_) || (lLambdaBar * lLambdaTilde < zeroTol_) || (lLambdaBar * eigenValues[0] * eigenValues[1] < zeroTol_) || (std::min(std::fabs(eigenValues[0]), std::fabs(eigenValues[2])) < zeroTol_) || (eigenValues[0] < T(0.0) && middleStateSpeeds_[0] > T(0.0)) || (eigenValues[2] > T(0.0) && middleStateSpeeds_[1] < T(0.0)) || ((uLeft_ + sqrt(gravity_ * hLeft_)) * (uRight_ + sqrt(gravity_ * hRight_)) < T(0.0)) || ((uLeft_ - sqrt(gravity_ * hLeft_)) * (uRight_ - sqrt(gravity_ * hRight_)) < T(0.0))) {
#endif
        steadyStateWave[0] = -(bRight_ - bLeft_);
        steadyStateWave[1] = -gravity_ * hBar * (bRight_ - bLeft_);
#ifdef COMPLEX_STEADY_STATE_WAVE
      } else {
        steadyStateWave[0] = gravity_ * (hBar / lLambdaBar) * (bRight_ - bLeft_);
        T hTilde           = hBar * lLambdaTilde / lLambdaBar;

        // Set bounds for problems far from steady state
        hTilde             = std::max(std::min(hLeft_, hRight_), hTilde);
        hTilde             = std::min(std::max(hLeft_, hRight_), hTilde);
        steadyStateWave[1] = -gravity_ * hTilde * (bRight_ - bLeft_);
      }
#endif

      // Preserve depth-positivity
      //   \cite[ch. 6.5.2]{george2006finite}, \cite[ch. 4.2.3]{george2008augmented}
      if (eigenValues[0] < -zeroTol_ && eigenValues[2] > zeroTol_) { // Subsonic
        steadyStateWave[0] = std::max(
          steadyStateWave[0], hLLMiddleHeight * (eigenValues[2] - eigenValues[0]) / eigenValues[0]
        );
        steadyStateWave[0] = std::min(
          steadyStateWave[0], hLLMiddleHeight * (eigenValues[2] - eigenValues[0]) / eigenValues[2]
        );
      } else if (eigenValues[0] > zeroTol_) { // Supersonic right TODO: motivation?
        steadyStateWave[0] = std::max(steadyStateWave[0], -hLeft_);
        steadyStateWave[0] = std::min(
          steadyStateWave[0], hLLMiddleHeight * (eigenValues[2] - eigenValues[0]) / eigenValues[0]
        );
      } else if (eigenValues[2] < -zeroTol_) { // Supersonic left TODO: motivation?
        steadyStateWave[0] = std::max(
          steadyStateWave[0], hLLMiddleHeight * (eigenValues[2] - eigenValues[0]) / eigenValues[2]
        );
        steadyStateWave[0] = std::min(steadyStateWave[0], hRight_);
      }

      // Limit the effect of the source term
      //   \cite[ch. 6.4.2]{george2006finite}
      steadyStateWave[1] = std::min(
        steadyStateWave[1], gravity_ * std::max(-hLeft_ * (bRight_ - bLeft_), -hRight_ * (bRight_ - bLeft_))
      );
      steadyStateWave[1] = std::max(
        steadyStateWave[1], gravity_ * std::min(-hLeft_ * (bRight_ - bLeft_), -hRight_ * (bRight_ - bLeft_))
      );

      // No source term in the case of a wall
      if (wetDryState_ == WavePropagationSolver<T>::WetDryState::WetDryWall || wetDryState_ == WavePropagationSolver<T>::WetDryState::DryWetWall) {
        assert(std::fabs(steadyStateWave[0]) < zeroTol_);
        assert(std::fabs(steadyStateWave[1]) < zeroTol_);
      }

      rightHandSide[0] -= steadyStateWave[0];
      // rightHandSide[1]: no source term
      rightHandSide[2] -= steadyStateWave[1];

      // Everything is ready, solve the equations!
      T beta[3]{};
      solveLinearEquation(eigenVectors, rightHandSide, beta);

      // Compute f-waves and wave-speeds
      if (wetDryState_ == WavePropagationSolver<T>::WetDryState::WetDryWall) { // Zero ghost updates (wall boundary)
        // Care about the left going wave (0) only
        o_fWaves[0][0] = beta[0] * eigenVectors[1][0];
        o_fWaves[0][1] = beta[0] * eigenVectors[2][0];

        // Set the rest to zero
        o_fWaves[1][0] = o_fWaves[1][1] = T(0.0);
        o_fWaves[2][0] = o_fWaves[2][1] = T(0.0);

        o_waveSpeeds[0] = eigenValues[0];
        o_waveSpeeds[1] = o_waveSpeeds[2] = T(0.0);

        assert(eigenValues[0] < zeroTol_);
      } else if (wetDryState_ == WavePropagationSolver<T>::WetDryState::DryWetWall) { // Zero ghost updates (wall
                                                                                      // boundary)
        // Care about the right going wave (2) only
        o_fWaves[2][0] = beta[2] * eigenVectors[1][2];
        o_fWaves[2][1] = beta[2] * eigenVectors[2][2];

        // Set the rest to zero
        o_fWaves[0][0] = o_fWaves[0][1] = T(0.0);
        o_fWaves[1][0] = o_fWaves[1][1] = T(0.0);

        o_waveSpeeds[2] = eigenValues[2];
        o_waveSpeeds[0] = o_waveSpeeds[1] = 0.;

        assert(eigenValues[2] > -zeroTol_);
      } else {
        // Compute f-waves (default)
        for (int waveNumber = 0; waveNumber < 3; waveNumber++) {
          o_fWaves[waveNumber][0] = beta[waveNumber] * eigenVectors[1][waveNumber]; // Select 2nd and
          o_fWaves[waveNumber][1] = beta[waveNumber] * eigenVectors[2][waveNumber]; // 3rd component of the augmented
                                                                                    // decomposition
#if ENABLE_AUGMENTED_RIEMANN_EIGEN_COEFFICIENTS
          // Store eigen-coefficients
          o_eigenCoefficients[waveNumber] = beta[waveNumber];
#endif
        }

        o_waveSpeeds[0] = eigenValues[0];
        o_waveSpeeds[1] = eigenValues[1];
        o_waveSpeeds[2] = eigenValues[2];
      }
    }

    /**
     * Computes the middle state of the homogeneous Riemann-problem.
     *   -> (\cite[ch. 13]{leveque2002finite})
     *
     * @param hLeft height on the left side of the edge.
     * @param hRight height on the right side of the edge.
     * @param uLeft velocity on the left side of the edge.
     * @param uRight velocity on the right side of the edge.
     * @param huLeft momentum on the left side of the edge.
     * @param huRight momentum on the right side of the edge.
     * @param maxNumberOfNewtonIterations maximum number of Newton iterations.
     */
    void computeMiddleState(
      const T&   hLeft,
      const T&   hRight,
      const T&   uLeft,
      const T&   uRight,
      const T&   huLeft,
      const T&   huRight,
      const int& maxNumberOfNewtonIterations = 1
    ) {
      // Set everything to zero
      hMiddle_              = T(0.0);
      middleStateSpeeds_[0] = T(0.0);
      middleStateSpeeds_[1] = T(0.0);

      // Compute local square roots
      // (not necessarily the same ones as in computeNetUpdates!)
      T lsqrtGhRight = std::sqrt(gravity_ * hRight);
      T lsqrtGhLeft  = std::sqrt(gravity_ * hLeft);

      // Single rarefaction in the case of a wet/dry interface
      if (hLeft < dryTol_) {
        middleStateSpeeds_[1] = middleStateSpeeds_[0] = uRight - T(2) * lsqrtGhRight;
        riemannStructure_                             = DrySingleRarefaction;
        return;
      } else if (hRight < dryTol_) {
        middleStateSpeeds_[0] = middleStateSpeeds_[1] = uLeft + T(2) * lsqrtGhLeft;
        riemannStructure_                             = SingleRarefactionDry;
        return;
      }

      // Determine the wave structure of the Riemann-problem
      riemannStructure_ = determineRiemannStructure(hLeft, hRight, uLeft, uRight);

      // Will be computed later
      T sqrtGhMiddle = T(0.0);

      if (riemannStructure_ == ShockShock) {
        /* Compute All-Shock Riemann Solution
         *  \cite[ch. 13.7]{leveque2002finite}
         *
         * Compute middle state h_m
         * => Solve non-linear scalar equation by Newton's method
         *    u_l - (h_m - h_l) \sqrt{\frac{g}{2} \left( \frac{1}{h_m} + \frac{1}{h_l} \right) }
         * =  u_r + (h_m - h_r) \sqrt{\frac{g}{2} \left( \frac{1}{h_m} + \frac{1}{h_r} \right) }
         *
         * Therefore determine the root of phi_{ss}:
         *\begin{equation}
         *  \begin{matrix}
         *    \phi_{ss}(h) =&&  u_r - u_l +
         *                  (h-h_l) \sqrt{
         *                            \frac{g}{2} \left( \frac{1}{h} + \frac{1}{h_l} \right)
         *                          } \\
         *             &&  +(h-h_r) \sqrt{
         *                            \frac{g}{2} \left( \frac{1}{h} + \frac{1}{h_r} \right)
         *                          }
         *  \end{matrix}
         *\end{equation}
         *
         *\begin{equation}
         *  \begin{matrix}
         *    \phi_{ss}'(h) = &&\sqrt{\frac{g}{2} \left( \frac{1}{h} + \frac{1}{h_l} \right) }
         *                      +\sqrt{\frac{g}{2} \left( \frac{1}{h} + \frac{1}{h_r} \right) }\\
         *                    && -\frac{g}{4}
         *                        \frac{h-h_l}
         *                        {
         *                         h^2\sqrt{\frac{g}{2} \left( \frac{1}{h} + \frac{1}{h_l} \right) }
         *                        }-
         *                        \frac{g}{4}
         *                        \frac{h-h_r}
         *                        {
         *                         h^2\sqrt{\frac{g}{2} \left( \frac{1}{h} + \frac{1}{h_r} \right) }
         *                        }
         *  \end{matrix}
         *\end{equation}
         */

        // hMiddle_ = (hLeft + hRight) * T(0.618); // First estimate
        hMiddle_ = std::min(hLeft, hRight);

        T lsqrtTermH[2] = {0, 0};

        for (int i = 0; i < maxNumberOfNewtonIterations; i++) {
          // sqrtTermHLow = std::sqrt(T(0.5) * gravity_ * (T(1/hMiddle_) + (T(1/hLeft)));
          lsqrtTermH[0] = std::sqrt(T(0.5) * gravity_ * ((hMiddle_ + hLeft) / (hMiddle_ * hLeft)));
          // sqrtTermHHigh = std::sqrt(T(0.5) * gravity_ * (T(1/hMiddle_) + T(1/hRight)));
          lsqrtTermH[1] = std::sqrt(T(0.5) * gravity_ * ((hMiddle_ + hRight) / (hMiddle_ * hRight)));

          T phi = uRight - uLeft + (hMiddle_ - hLeft) * lsqrtTermH[0] + (hMiddle_ - hRight) * lsqrtTermH[1];

          if (std::fabs(phi) < newtonTol_) {
            break;
          }

          T derivativePhi = lsqrtTermH[0] + lsqrtTermH[1]
                            - T(0.25) * gravity_ * (hMiddle_ - hLeft) / (lsqrtTermH[0] * hMiddle_ * hMiddle_)
                            - T(0.25) * gravity_ * (hMiddle_ - hRight) / (lsqrtTermH[1] * hMiddle_ * hMiddle_);

          hMiddle_ = hMiddle_ - phi / derivativePhi; // Newton step
          assert(hMiddle >= dryTol_);

          // if (i == maxNumberOfNewtonIterations_ - 1) {
          //   std::cerr << "Newton-Method did not converge" << std::endl;
          //   std::cerr << "std::fabs(phi): " << std::fabs(phi) << std::endl;
          //   std::cerr << "hMiddle: " << hMiddle_ << std::endl;
          //   assert(false);
          // }
        }

        sqrtGhMiddle = std::sqrt(gravity_ * hMiddle_);

        // Compute middle speed u_m
        //\begin{equation}
        //   \label{eq:hmshock1stfamsimp}
        //   u_m = u_l - (h_m - h_l) \sqrt{\frac{g}{2} \left( \frac{1}{h_m} + \frac{1}{h_l} \right) }
        //\end{equation}
        //
        //\begin{equation}
        //   \label{eq:hmshock2ndfamsimp}
        //   u_m = u_r + (h_m - h_r) \sqrt{\frac{g}{2} \left( \frac{1}{h_m} + \frac{1}{h_r} \right) }
        //\end{equation}

        // T uMiddleEstimates[2];
        // uMiddleEstimates[0] = uLeft - (hMiddle_ - hLeft) * lsqrtTermH[0];
        // uMiddleEstimates[1] = uRight + (hMiddle_ - hRight) * lsqrtTermH[1];
        // uMiddle = T(0.5) * (uMiddleEstimates[0] + uMiddleEstimates[1]);

        // Middle state speeds as they are implemented in clawpack, TODO: why?
        //        middleStateSpeeds[0] = uMiddleEstimates[0] - sqrtGhMiddle;
        //        middleStateSpeeds[1] = uMiddleEstimates[1] + sqrtGhMiddle;
      }

      if (riemannStructure_ == RarefactionRarefaction) {
        // Compute All-Rarefaction Riemann Solution
        //   \cite[ch. 13.8.6]{leveque2002finite}

        // Compute middle state height h_m
        hMiddle_ = std::max(
          T(0.0), uLeft - uRight + T(2.0) * (lsqrtGhLeft + lsqrtGhRight)
        ); // std::max -> Text after formula (13.56), page 279
        hMiddle_ = T(1) / (T(16) * gravity_) * hMiddle_ * hMiddle_;

        sqrtGhMiddle = std::sqrt(gravity_ * hMiddle_);

        // Middle state speeds as they are implemented in clawpack, why?
        //        middleStateSpeeds[0] = uLeft + T(2.0) * lsqrtGhLow - T(3.0) * sqrtGhMiddle;
        //        middleStateSpeeds[1] = uRight - T(2.0) * lsqrtGhHigh + T(3.0) * sqrtGhMiddle;
      }

      if (riemannStructure_ == ShockRarefaction || riemannStructure_ == RarefactionShock) {
        // Compute dam-break Riemann-solution
        // TODO: reference
        T hMin, hMax;
        if (riemannStructure_ == ShockRarefaction) {
          hMin = hLeft;
          hMax = hRight;
        } else {
          hMin = hRight;
          hMax = hLeft;
        }

        // hMiddle_ = (hLeft + hRight) * T(0.618); // First estimate
        hMiddle_ = hMin;

        sqrtGhMiddle = std::sqrt(gravity_ * hMiddle_);
        T sqrtGhMax  = std::sqrt(gravity_ * hMax);
        for (int i = 0; i < maxNumberOfNewtonIterations; i++) {
          /*
           * Compute middle state h_m
           * => Solve non-linear scalar equation by Newton's method
           *
           * Therefore determine the root of phi_{sr}/phi_{rs}:
           *\begin{equation}
           *  \begin{matrix}
           *    h_{min} = \min(h_l, h_r)\\
           *    h_{max} = \max(h_l, h_r)
           *  \end{matrix}
           *\end{equation}
           *\begin{equation}
           *  \begin{matrix}
           *    \phi_{sr}(h) = \phi_{rs}(h) =&& u_r - u_l + (h - h_{min}) \sqrt{
           *                                            \frac{g}{2} \left( \frac{1}{h} + \frac{1}{h_{min}} \right)
           *                                          } \\
           *               && + 2 (\sqrt{gh} - \sqrt{gh_{max}})
           *  \end{matrix}
           *\end{equation}
           *\begin{equation}
           *    \phi'_{sr}(h) = \phi'_{rs}(h) = \sqrt{\frac{g}{2} \left( \frac{1}{h} + \frac{1}{h_{min}} \right)  }
           *                     -\frac{g}{4} \frac{h-h_{min}}
           *                        { h^2
           *                          \sqrt{
           *                            \frac{g}{2}  \left(\frac{1}{h} + \frac{1}{h_{min}} \right)
           *                          }
           *                        }
           *                     + \sqrt{\frac{g}{h}}
           *\end{equation}
           */

          T sqrtTermHMin = std::sqrt(T(0.5) * gravity_ * ((hMiddle_ + hMin) / (hMiddle_ * hMin)));

          T phi = uRight - uLeft + (hMiddle_ - hMin) * sqrtTermHMin + T(2) * (sqrtGhMiddle - sqrtGhMax);

          if (std::fabs(phi) < newtonTol_) {
            break;
          }

          T derivativePhi = sqrtTermHMin - T(0.25) * gravity_ * (hMiddle_ - hMin) / (hMiddle_ * hMiddle_ * sqrtTermHMin)
                            + sqrtG_ / sqrtGhMiddle;

          hMiddle_ = hMiddle_ - phi / derivativePhi; // Newton step

#ifndef NDEBUG
          if (hMiddle_ < hMin) {
            std::cout << phi << std::endl;
            std::cout << derivativePhi << std::endl;
            std::cerr << "hMiddle(" << hMiddle_ << ") < hMin(" << hMin << ")" << std::endl;
            assert(false);
          }
#endif

          // if (i == maxNumberOfNewtonIterations - 1) {
          //   std::cerr << "Newton-Method did not converge" << std::endl;
          //   std::cerr << "std::fabs(phi): " << std::fabs(phi) << std::endl;
          //   std::cerr << "hMiddle: " << hMiddle_ << std::endl;
          //   assert(false);
          // }

          sqrtGhMiddle = std::sqrt(gravity_ * hMiddle_);
        }

        // Middle state speeds as they are implemented in clawpack, TODO: why?
        // if (riemannStructure_ == ShockRarefaction) {
        //   uMiddle = uRight - T(2.0) * (lsqrtGhHigh - sqrtGhMiddle);
        //   middleStateSpeeds[0] = uRight - T(2.0) * lsqrtGhHigh + sqrtGhMiddle;
        //   middleStateSpeeds[1] = uRight - T(2.0) * lsqrtGhHigh + T(3.0) * sqrtGhMiddle;
        // } else {
        //   uMiddle = uLeft + T(2.0) * (lsqrtGhLow - sqrtGhMiddle);
        //   middleStateSpeeds[0] = uLeft + T(2.0) * lsqrtGhLow - T(3.0) * sqrtGhMiddle;
        //   middleStateSpeeds[1] = uLeft + T(2.0) * lsqrtGhLow - sqrtGhMiddle;
        // }
      }

      middleStateSpeeds_[0] = uLeft + T(2.0) * lsqrtGhLeft - T(3.0) * sqrtGhMiddle;
      middleStateSpeeds_[1] = uRight - T(2.0) * lsqrtGhRight + T(3.0) * sqrtGhMiddle;

      assert(hMiddle_ >= 0);
    }

    /**
     * Determine the Riemann-structure of a given problem.
     *   -> \cite[theorem 4.2]{george2006finite}, \cite[appendix B]{george2008augmented}
     *
     * @param hLeft height on the left side of the edge.
     * @param hRight height on the right side of the edge.
     * @param uLeft velocity on the left side of the edge.
     * @param uRight velocity on the right side of the edge.
     *
     * @return Riemann-structure of a given problem.
     */
    RiemannStructure determineRiemannStructure(const T& hLeft, const T& hRight, const T& uLeft, const T& uRight) const {
      T hMin = std::min(hLeft, hRight);
      T hMax = std::max(hLeft, hRight);

      T uDif = uRight - uLeft;

      if (0 <= T(2.0) * (std::sqrt(gravity_ * hMin) - std::sqrt(gravity_ * hMax)) + uDif) {
        return RarefactionRarefaction;
      }

      if ((hMax - hMin) * std::sqrt(gravity_ * T(0.5) * (1 / hMax + 1 / hMin)) + uDif <= 0) {
        return ShockShock;
      }

      if (hLeft < hRight) {
        return ShockRarefaction;
      }

      return RarefactionShock;
    }

    /**
     * Solve the linear equation:
     * A * x = b with A \in \mathbb{R}^{3\times3}, x,b \in \mathbb{R}^3
     *
     * @param matrix the matrix
     * @param b right hand side
     * @param o_x solution
     */
    static void solveLinearEquation(const T matrix[3][3], const T b[3], T o_x[3]) {
#if 1
      // Compute inverse of 3x3 matrix
      const T m[3][3] = {
        {(matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]),
         -(matrix[0][1] * matrix[2][2] - matrix[0][2] * matrix[2][1]),
         (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1])},
        {-(matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]),
         (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]),
         -(matrix[0][0] * matrix[1][2] - matrix[0][2] * matrix[1][0])},
        {(matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]),
         -(matrix[0][0] * matrix[2][1] - matrix[0][1] * matrix[2][0]),
         (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) }};
      T d = (matrix[0][0] * m[0][0] + matrix[0][1] * m[1][0] + matrix[0][2] * m[2][0]);

#ifndef NDEBUG
      if (std::fabs(d) < 0.000000001) {
        std::cerr << "Division close to zero!" << std::endl;
        std::cout << "Matrix: " << std::endl;
        std::cout << matrix[0][0] << ", " << matrix[0][1] << ", " << matrix[0][2] << std::endl;
        std::cout << matrix[1][0] << ", " << matrix[1][1] << ", " << matrix[1][2] << std::endl;
        std::cout << matrix[2][0] << ", " << matrix[2][1] << ", " << matrix[2][2] << std::endl;
        std::cout << "b: " << std::endl;
        std::cout << b[0] << ", " << b[1] << ", " << b[2] << std::endl;
        std::cout << "d: " << d << std::endl;
        assert(false);
      }
#endif

      // m stores not really the inverse matrix, but the inverse multiplied by d
      T s = T(1) / d;

      // Compute m*b
      o_x[0] = (m[0][0] * b[0] + m[0][1] * b[1] + m[0][2] * b[2]) * s;
      o_x[1] = (m[1][0] * b[0] + m[1][1] * b[1] + m[1][2] * b[2]) * s;
      o_x[2] = (m[2][0] * b[0] + m[2][1] * b[1] + m[2][2] * b[2]) * s;

#else
      T origDet = computeDeterminant(matrix);

      if (std::fabs(origDet) > zeroTol_) {
        T modifiedMatrix[3][3]{};

        for (int column = 0; column < 3; column++) {
          memcpy(modifiedMatrix, matrix, sizeof(T) * 3 * 3);

          // Set a column of the matrix to b
          modifiedMatrix[0][column] = b[0];
          modifiedMatrix[1][column] = b[1];
          modifiedMatrix[2][column] = b[2];

          o_x[column] = computeDeterminant(modifiedMatrix) / origDet;
        }
      } else {
        std::cerr << "Warning: Linear dependent eigenvectors! (using Jacobi Solver)" << std::endl;
        std::cerr << "Determinant: " << origDet << std::endl;
        std::cerr << matrix[0][0] << "\t" << matrix[0][1] << "\t" << matrix[0][2] << std::endl;
        std::cerr << matrix[1][0] << "\t" << matrix[1][1] << "\t" << matrix[1][2] << std::endl;
        std::cerr << matrix[2][0] << "\t" << matrix[2][1] << "\t" << matrix[2][2] << std::endl;
#if 1
        T xTemp[3]{};
        xTemp[0] = xTemp[1] = xTemp[2] = T(0.0);
        for (int m = 0; m < 20; m++) {
          for (int row = 0; row < 3; row++) {
            o_x[row] = 0.;
            for (int col = 0; col < 3; col++) {
              if (col != row) {
                o_x[row] += matrix[row][col] * xTemp[col];
              }
            }
            o_x[row] = (b[row] - o_x[row]) / matrix[row][row];
          }

          if (fabs(o_x[0] - xTemp[0]) + fabs(o_x[1] - xTemp[1]) + fabs(o_x[2] - xTemp[2]) < zeroTol_ * 10.) {
            break;
          } else {
            std::cout
              << "Error: " << fabs(o_x[0] - xTemp[0]) + fabs(o_x[1] - xTemp[1]) + fabs(o_x[2] - xTemp[2]) << std::endl;
            xTemp[0] = o_x[0];
            xTemp[1] = o_x[1];
            xTemp[2] = o_x[2];
          }
        }
        std::cout << "Solution:" << std::endl;
        std::cout << "\t" << o_x[0] << std::endl;
        std::cout << "x=\t" << o_x[1] << std::endl;
        std::cout << "\t" << o_x[2] << std::endl;
        std::cout << "*********" << std::endl;
        std::cout << "\t" << b[0] << std::endl;
        std::cout << "b=\t" << b[1] << std::endl;
        std::cout << "\t" << b[2] << std::endl;
        std::cout << "***A*x****" << std::endl;
        for (int row = 0; row < 3; row++) {
          std::cout << "\t" << matrix[row][0] * o_x[0] + matrix[row][1] * o_x[1] + matrix[row][2] * o_x[2] << std::endl;
        }
#endif
        assert(false);
      }
#endif
    }

#if 0
    /**
     * Compute the determinant of a 3x3 matrix
     * using formula: |A|=a11 * (a22a33-a32a23) + a12 * (a23a31-a33a21) + a13 * (a21a32-a31a22)
     *                |A|=a00 * (a11a22-a21a12) + a01 * (a12a20-a22a10) + a02 * (a10a21-a20a11)
     */
    T computeDeterminant(T matrix[3][3]) const {
      T determinant = matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[2][1] * matrix[1][2])
                      + matrix[0][1] * (matrix[1][2] * matrix[2][0] - matrix[2][2] * matrix[1][0])
                      + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[2][0] * matrix[1][1]);
      return determinant;
    }
#endif

#undef sqrtGhLeft
#undef sqrtGhRight
#undef sqrthLeft
#undef sqrthRight
  };

} // namespace Solvers

#undef COMPLEX_STEADY_STATE_WAVE
#undef CORRECT_RARE_FACTIONS
