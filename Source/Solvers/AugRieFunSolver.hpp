/**
 * AugRieFun.hpp
 *
 ****
 **** Approximate Augmented Riemann Solver for the Shallow Water Equations
 **** functional implementation (based on AugRieCUDA)
 ****
 *
 *  Created on: May 28, 2013
 *  Last Update: Jan 1, 2014
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
 *  CUDA: Wolfgang Hölzl
 *    E-Mail: hoelzlw AT in.tum.de
 *
 *  Further optimizations: Michael Bader
 *    Homepage: http://www5.in.tum.de/wiki/index.php/Michael_Bader
 *    E-Mail: bader AT in.tum.de
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

#ifndef AUGRIE_FUN_HPP
#define AUGRIE_FUN_HPP

#include <cassert>
#include <cmath>
#include <iostream>
//#include <string>
//#include <vector>

// constants to classify wet-dry-state of a pairs of cells:
const int DryDry               = 0;
const int WetWet               = 1;
const int WetDryInundation     = 2;
const int WetDryWall           = 3;
const int WetDryWallInundation = 4;
const int DryWetInundation     = 5;
const int DryWetWall           = 6;
const int DryWetWallInundation = 7;

// constants to classify Riemann state of a pairs of cells:
const int DrySingleRarefaction   = 0;
const int SingleRarefactionDry   = 1;
const int ShockShock             = 2;
const int ShockRarefaction       = 3;
const int RarefactionShock       = 4;
const int RarefactionRarefaction = 5;

namespace solver {

  /**
   *
   */
  template <typename real>
  class AugRieFun {
  private:
    const real         dryTol;
    const real         g;      // gravity constant
    const real         half_g; // 0.5 * gravity constant
    const real         sqrt_g; // square root of the gravity constant
    const real         zeroTol;
    const real         newtonTol;                   // tolerance for the Newton iterative solver
    const unsigned int maxNumberOfNewtonIterations; // maximum number of performed Newton iterations

  public:
    /**
     * AugRieFun Constructor, takes three problem parameters
     * @param dryTol "dry tolerance": if the water height falls below dryTol, wall boundary conditions are applied
     * (default value is 100)
     * @param gravity takes the value of the gravity constant (default value is 9.81 m/s^2)
     * @param zeroTol computed f-waves with an absolute value < zeroTol are treated as static waves (default value is
     * 10^{-7})
     */
    AugRieFun(
      real i_dryTol        = (real)1.0,
      real i_gravity       = (real)9.81,
      real i_zeroTol       = (real)0.0000001,
      real i_newtonTol     = (real)0.0000001,
      real i_maxNewtonIter = 1
    ):
      dryTol(i_dryTol),
      g(i_gravity),
      half_g(static_cast<real>(.5) * i_gravity),
      sqrt_g(std::sqrt(i_gravity)),
      zeroTol(i_zeroTol),
      newtonTol(i_newtonTol),
      maxNumberOfNewtonIterations(i_maxNewtonIter) {}

    /**
     * Compute net updates for the left/right cell of the edge.
     *
     * maxWaveSpeed will be set to the maximum (linearized) wave speed => CFL
     */
    void computeNetUpdates(
      const real i_hLeft,
      const real i_hRight,
      const real i_huLeft,
      const real i_huRight,
      const real i_bLeft,
      const real i_bRight,

      real& o_hUpdateLeft,
      real& o_hUpdateRight,
      real& o_huUpdateLeft,
      real& o_huUpdateRight,
      real& o_maxWaveSpeed
    ) const {
      real hLeft   = i_hLeft;
      real hRight  = i_hRight;
      real uLeft   = static_cast<real>(0);
      real uRight  = static_cast<real>(0);
      real huLeft  = i_huLeft;
      real huRight = i_huRight;
      real bLeft   = i_bLeft;
      real bRight  = i_bRight;

      // declare variables which are used over and over again
      real sqrt_g_hLeft;
      real sqrt_g_hRight;

      real sqrt_hLeft;
      real sqrt_hRight;

      // reset net updates and the maximum wave speed
      o_hUpdateLeft = o_hUpdateRight = o_huUpdateLeft = o_huUpdateRight = static_cast<real>(0);
      o_maxWaveSpeed                                                    = static_cast<real>(0);

      // variables for computing middle states
      real hMiddle              = static_cast<real>(0);
      real middleStateSpeeds[2] = {static_cast<real>(0)};

      /***************************************************************************************
       * Determine Wet Dry State Begin
       * (determine the wet/dry state and compute local variables correspondingly)
       **************************************************************************************/
      int wetDryState;
      // compute speeds or set them to zero (dry cells)
      if (hLeft >= dryTol) {
        uLeft = huLeft / hLeft;
      } else {
        bLeft += hLeft;
        hLeft = huLeft = uLeft = static_cast<real>(0);
      }

      if (hRight >= dryTol) {
        uRight = huRight / hRight;
      } else {
        bRight += hRight;
        hRight = huRight = uRight = static_cast<real>(0);
      }

      // MB: determine wet/dry-state - try to start with most frequent case
      if (hLeft >= dryTol) {
        if (hRight >= dryTol) {
          // simple wet/wet case - expected as most frequently executed branch
          wetDryState = WetWet;
        } else { // hLeft >= dryTol and hRight < dryTol
          // we have a shoreline: left cell wet, right cell dry
          //=>check for simple inundation problems
          if (hLeft + bLeft > bRight) {
            // => dry cell lies lower than the wet cell
            wetDryState = WetDryInundation;
          } else { // hLeft >= dryTol and hRight < dryTol and hLeft + bLeft <= bRight
            // =>dry cell (right) lies higher than the wet cell (left)
            // lets check if the momentum is able to overcome the difference in height
            //  => solve homogeneous Riemann-problem to determine the middle state height
            //     which would arise if there is a wall (wall-boundary-condition)
            //       \cite[ch. 6.8.2]{george2006finite})
            //       \cite[ch. 5.2]{george2008augmented}
            computeMiddleState(hLeft, hLeft, uLeft, -uLeft, maxNumberOfNewtonIterations, hMiddle, middleStateSpeeds);

            if (hMiddle + bLeft > bRight) {
              // momentum is large enough, continue with the original values
              //           bRight = o_hMiddle + bLeft;
              wetDryState = WetDryWallInundation;
              // limit the effect of the source term if there is a "wall"
              //\cite[end of ch. 5.2?]{george2008augmented}, \cite[rpn2ez_fast_geo.f][levequeclawpack]
              bRight = hLeft + bLeft;
            } else {
              hRight  = hLeft;
              uRight  = -uLeft;
              huRight = -huLeft;
              bRight = bLeft = static_cast<real>(0);
              wetDryState    = WetDryWall;
            }
          }
        }
      } else { // hLeft < dryTol
        if (hRight >= dryTol) {
          // we have a shoreline: left cell dry, right cell wet
          //=>check for simple inundation problems
          if (hRight + bRight > bLeft) {
            // => dry cell lies lower than the wet cell
            wetDryState = DryWetInundation;
          } else { // hLeft < dryTol and hRight >= dryTol and hRight + bRight <= bLeft
            // =>dry cell (left) lies higher than the wet cell (right)
            // lets check if the momentum is able to overcome the difference in height
            //  => solve homogeneous Riemann-problem to determine the middle state height
            //     which would arise if there is a wall (wall-boundary-condition)
            //       \cite[ch. 6.8.2]{george2006finite})
            //       \cite[ch. 5.2]{george2008augmented}
            computeMiddleState(
              hRight, hRight, -uRight, uRight, maxNumberOfNewtonIterations, hMiddle, middleStateSpeeds
            );

            if (hMiddle + bRight > bLeft) {
              // momentum is large enough, continue with the original values
              //           bLeft = o_hMiddle + bRight;
              wetDryState = DryWetWallInundation;
              // limit the effect of the source term if there is a "wall"
              //\cite[end of ch. 5.2?]{george2008augmented}, \cite[rpn2ez_fast_geo.f][levequeclawpack]
              bLeft = hRight + bRight;
            } else {
              // momentum is not large enough, use wall-boundary-values
              hLeft  = hRight;
              uLeft  = -uRight;
              huLeft = -huRight;
              bLeft = bRight = static_cast<real>(0);
              wetDryState    = DryWetWall;
            }
          }
        } else { // hLeft < dryTol and hRight < dryTol
          wetDryState = DryDry;
          // nothing to do for dry/dry case, all netUpdates and maxWaveSpeed are 0
          return;
        }
      };

      /***************************************************************************************
       * Determine Wet Dry State End
       **************************************************************************************/
      // MB: not executed for DryDry => returns right after case is identified
      assert(wetDryState != DryDry);

      // precompute some terms which are fixed during
      // the computation after some specific point
      sqrt_hLeft  = std::sqrt(hLeft);
      sqrt_hRight = std::sqrt(hRight);

      sqrt_g_hLeft  = sqrt_g * sqrt_hLeft;
      sqrt_g_hRight = sqrt_g * sqrt_hRight;

      // compute the augmented decomposition
      //   (thats the place where the computational work is done..)
      /***************************************************************************************
       * Compute Wave Decomposition Begin
       **************************************************************************************/
      // compute eigenvalues of the jacobian matrices in states Q_{i-1} and Q_{i} (char. speeds)
      const real characteristicSpeeds[2] = {uLeft - sqrt_g_hLeft, uRight + sqrt_g_hRight};

      // compute "Roe speeds"
      const real hRoe = static_cast<real>(0.5) * (hRight + hLeft);
      const real uRoe = (uLeft * sqrt_hLeft + uRight * sqrt_hRight) / (sqrt_hLeft + sqrt_hRight);

      // optimization for dumb compilers
      const real sqrt_g_hRoe  = sqrt_g * std::sqrt(hRoe);
      const real roeSpeeds[2] = {uRoe - sqrt_g_hRoe, uRoe + sqrt_g_hRoe};

      // compute the middle state of the homogeneous Riemann-Problem
      //  MB: ???computeMiddleState always called? -> possible to rearrange computation to reduce if-statements???
      //  MB: !!!here always only 1 Newton iteration
      if (wetDryState != WetDryWall && wetDryState != DryWetWall) {
        // case WDW and DWW was computed in determineWetDryState already
        computeMiddleState(hLeft, hRight, uLeft, uRight, 1, hMiddle, middleStateSpeeds);
      }

      // compute extended eindfeldt speeds (einfeldt speeds + middle state speeds)
      //   \cite[ch. 5.2]{george2008augmented}, \cite[ch. 6.8]{george2006finite}
      real extEinfeldtSpeeds[2] = {static_cast<real>(0), static_cast<real>(0)};
      if (wetDryState == WetWet || wetDryState == WetDryWall || wetDryState == DryWetWall) {
        extEinfeldtSpeeds[0] = std::min(characteristicSpeeds[0], roeSpeeds[0]);
        extEinfeldtSpeeds[0] = std::min(extEinfeldtSpeeds[0], middleStateSpeeds[1]);

        extEinfeldtSpeeds[1] = std::max(characteristicSpeeds[1], roeSpeeds[1]);
        extEinfeldtSpeeds[1] = std::max(extEinfeldtSpeeds[1], middleStateSpeeds[0]);
      } else if (hLeft < dryTol) { // MB: !!! cases DryWetInundation, DryWetWallInundation
        // ignore undefined speeds
        extEinfeldtSpeeds[0] = std::min(roeSpeeds[0], middleStateSpeeds[1]);
        extEinfeldtSpeeds[1] = std::max(characteristicSpeeds[1], roeSpeeds[1]);

        assert(middleStateSpeeds[0] < extEinfeldtSpeeds[1]);
      } else if (hRight < dryTol) { // MB: !!! cases WetDryInundation, WetDryWallInundation
        // ignore undefined speeds
        extEinfeldtSpeeds[0] = std::min(characteristicSpeeds[0], roeSpeeds[0]);
        extEinfeldtSpeeds[1] = std::max(roeSpeeds[1], middleStateSpeeds[0]);

        assert(middleStateSpeeds[1] > extEinfeldtSpeeds[0]);
      } else {
        assert(false);
      }

      // HLL middle state
      //   \cite[theorem 3.1]{george2006finite}, \cite[ch. 4.1]{george2008augmented}
      const real hLLMiddleHeight = std::max(
        (huLeft - huRight + extEinfeldtSpeeds[1] * hRight - extEinfeldtSpeeds[0] * hLeft
        ) / (extEinfeldtSpeeds[1] - extEinfeldtSpeeds[0]),
        static_cast<real>(0)
      );

      // define eigenvalues
      const real eigenValues[3] = {
        extEinfeldtSpeeds[0],
        static_cast<real>(0.5) * (extEinfeldtSpeeds[0] + extEinfeldtSpeeds[1]),
        extEinfeldtSpeeds[1]};

      // define eigenvectors
      //  MB: no longer used as system matrix
      //      -> but still used to compute f-waves
      const real eigenVectors[3][3] = {
        {static_cast<real>(1), static_cast<real>(0), static_cast<real>(1)},
        {eigenValues[0], static_cast<real>(0), eigenValues[2]},
        {eigenValues[0] * eigenValues[0], static_cast<real>(1), eigenValues[2] * eigenValues[2]}};

      // compute rarefaction corrector wave
      //   \cite[ch. 6.7.2]{george2006finite}, \cite[ch. 5.1]{george2008augmented}
#pragma message \
  "strongRarefaction set to false. Further investigations needed about senseless initialization of eigenValues[1]"

      // compute the jump in state
      real rightHandSide[3] = {
        hRight - hLeft,
        huRight - huLeft,
        (huRight * uRight + half_g * hRight * hRight) - (huLeft * uLeft + half_g * hLeft * hLeft)};

      // compute steady state wave
      //   \cite[ch. 4.2.1 \& app. A]{george2008augmented}, \cite[ch. 6.2 \& ch. 4.4]{george2006finite}
      real steadyStateWave[2] = {-(bRight - bLeft), -half_g * (hLeft + hRight) * (bRight - bLeft)};

      // preserve depth-positivity
      //   \cite[ch. 6.5.2]{george2006finite}, \cite[ch. 4.2.3]{george2008augmented}
      if (eigenValues[0] < -zeroTol && eigenValues[2] > zeroTol) {
        // subsonic
        steadyStateWave[0] = std::max(
          steadyStateWave[0], hLLMiddleHeight * (eigenValues[2] - eigenValues[0]) / eigenValues[0]
        );
        steadyStateWave[0] = std::min(
          steadyStateWave[0], hLLMiddleHeight * (eigenValues[2] - eigenValues[0]) / eigenValues[2]
        );
      } else if (eigenValues[0] > zeroTol) {
        // supersonic right TODO: motivation?
        steadyStateWave[0] = std::max(steadyStateWave[0], -hLeft);
        steadyStateWave[0] = std::min(
          steadyStateWave[0], hLLMiddleHeight * (eigenValues[2] - eigenValues[0]) / eigenValues[0]
        );
      } else if (eigenValues[2] < -zeroTol) {
        // supersonic left TODO: motivation?
        steadyStateWave[0] = std::max(
          steadyStateWave[0], hLLMiddleHeight * (eigenValues[2] - eigenValues[0]) / eigenValues[2]
        );
        steadyStateWave[0] = std::min(steadyStateWave[0], hRight);
      }

      // Limit the effect of the source term
      //   \cite[ch. 6.4.2]{george2006finite}
      steadyStateWave[1] = std::min(
        steadyStateWave[1], g * std::max(-hLeft * (bRight - bLeft), -hRight * (bRight - bLeft))
      );
      steadyStateWave[1] = std::max(
        steadyStateWave[1], g * std::min(-hLeft * (bRight - bLeft), -hRight * (bRight - bLeft))
      );

      rightHandSide[0] -= steadyStateWave[0];
      // rightHandSide[1]: no source term
      rightHandSide[2] -= steadyStateWave[1];

      // everything is ready, solve the equations!
      /***************************************************************************************
       * Solve linear equation begin
       **************************************************************************************/

      // direct solution of specific 3x3 system:
      // (       1           0       1           ) ( beta[0] ) = ( rightHandSide[0] )
      // ( eigenValues[0]    0  eigenValues[2]   ) ( beta[1] )   ( rightHandSide[1] )
      // ( eigenValues[0]^2  1  eigenValues[2]^2 ) ( beta[2] )   ( rightHandSide[2] )

      // step 1: solve the following 2x2 system (1st and 3rd column):
      // (       1              1           ) ( beta[0] ) = ( rightHandSide[0] )
      // ( eigenValues[0]  eigenValues[2]   ) ( beta[2] )   ( rightHandSide[1] )

      // compute the inverse of the wave speed difference:
      real inverseDiff = static_cast<real>(1.) / (eigenValues[2] - eigenValues[0]); // 2 FLOPs (1 div)
      // compute f-waves:
      real beta[3];
      beta[0] = (eigenValues[2] * rightHandSide[0] - rightHandSide[1]) * inverseDiff;  // 3 FLOPs
      beta[2] = (-eigenValues[0] * rightHandSide[0] + rightHandSide[1]) * inverseDiff; // 3 FLOPs

      // step 2: solve 3rd row for beta[1]:
      beta[1] = rightHandSide[2] - eigenValues[0] * eigenValues[0] * beta[0] - eigenValues[2] * eigenValues[2] * beta[2];

#ifdef false
      // compute inverse of 3x3 matrix
      // MB: !!!AVOID COMPUTATION of inverse => see above for better solution for simple matrix
      // MB: TODO -> different choice of matrix in AugRie.hpp -> requires a separate solution procedure
      const real m[3][3] = {
        {(eigenVectors[1][1] * eigenVectors[2][2] - eigenVectors[1][2] * eigenVectors[2][1]),
         -(eigenVectors[0][1] * eigenVectors[2][2] - eigenVectors[0][2] * eigenVectors[2][1]),
         (eigenVectors[0][1] * eigenVectors[1][2] - eigenVectors[0][2] * eigenVectors[1][1])},
        {-(eigenVectors[1][0] * eigenVectors[2][2] - eigenVectors[1][2] * eigenVectors[2][0]),
         (eigenVectors[0][0] * eigenVectors[2][2] - eigenVectors[0][2] * eigenVectors[2][0]),
         -(eigenVectors[0][0] * eigenVectors[1][2] - eigenVectors[0][2] * eigenVectors[1][0])},
        {(eigenVectors[1][0] * eigenVectors[2][1] - eigenVectors[1][1] * eigenVectors[2][0]),
         -(eigenVectors[0][0] * eigenVectors[2][1] - eigenVectors[0][1] * eigenVectors[2][0]),
         (eigenVectors[0][0] * eigenVectors[1][1] - eigenVectors[0][1] * eigenVectors[1][0])}};
      const real d = (eigenVectors[0][0] * m[0][0] + eigenVectors[0][1] * m[1][0] + eigenVectors[0][2] * m[2][0]);

      // m stores not really the inverse matrix, but the inverse multiplied by d
      const real s = 1 / d;

      // compute m*rightHandSide
      const real beta[3] = {
        (m[0][0] * rightHandSide[0] + m[0][1] * rightHandSide[1] + m[0][2] * rightHandSide[2]) * s,
        (m[1][0] * rightHandSide[0] + m[1][1] * rightHandSide[1] + m[1][2] * rightHandSide[2]) * s,
        (m[2][0] * rightHandSide[0] + m[2][1] * rightHandSide[1] + m[2][2] * rightHandSide[2]) * s};
#endif

      /***************************************************************************************
       * Solve linear equation end
       **************************************************************************************/

      // compute f-waves and wave-speeds
      real fWaves[3][2];  // array to store the three f-wave vectors (2 components each)
      real waveSpeeds[3]; // and vector to their speeds

      if (wetDryState == WetDryWall) {
        // zero ghost updates (wall boundary)
        // care about the left going wave (0) only
        fWaves[0][0] = beta[0] * eigenVectors[1][0];
        fWaves[0][1] = beta[0] * eigenVectors[2][0];

        // set the rest to zero
        fWaves[1][0] = fWaves[1][1] = static_cast<real>(0);
        fWaves[2][0] = fWaves[2][1] = static_cast<real>(0);

        waveSpeeds[0] = eigenValues[0];
        waveSpeeds[1] = waveSpeeds[2] = static_cast<real>(0);

        assert(eigenValues[0] < zeroTol);
      } else if (wetDryState == DryWetWall) {
        // zero ghost updates (wall boundary)
        // care about the right going wave (2) only
        fWaves[2][0] = beta[2] * eigenVectors[1][2];
        fWaves[2][1] = beta[2] * eigenVectors[2][2];

        // set the rest to zero
        fWaves[0][0] = fWaves[0][1] = static_cast<real>(0);
        fWaves[1][0] = fWaves[1][1] = static_cast<real>(0);

        waveSpeeds[2] = eigenValues[2];
        waveSpeeds[0] = waveSpeeds[1] = static_cast<real>(0);

        assert(eigenValues[2] > -zeroTol);
      } else {
        // compute f-waves (default)
        for (int waveNumber = 0; waveNumber < 3; waveNumber++) {
          fWaves[waveNumber][0] = beta[waveNumber] * eigenVectors[1][waveNumber]; // select 2nd and
          fWaves[waveNumber][1] = beta[waveNumber] * eigenVectors[2][waveNumber]; // 3rd component of the augmented
                                                                                  // decomposition
        }

        waveSpeeds[0] = eigenValues[0];
        waveSpeeds[1] = eigenValues[1];
        waveSpeeds[2] = eigenValues[2];
      }
      /***************************************************************************************
       * Compute Wave Decomposition End
       **************************************************************************************/

      // compute the updates from the three propagating waves
      // A^-\delta Q = \sum{s[i]<0} \beta[i] * r[i] = A^-\delta Q = \sum{s[i]<0} Z^i
      // A^+\delta Q = \sum{s[i]>0} \beta[i] * r[i] = A^-\delta Q = \sum{s[i]<0} Z^i
      for (int waveNumber = 0; waveNumber < 3; waveNumber++) {
        if (waveSpeeds[waveNumber] < -zeroTol) {
          // left going
          o_hUpdateLeft += fWaves[waveNumber][0];
          o_huUpdateLeft += fWaves[waveNumber][1];
        } else if (waveSpeeds[waveNumber] > zeroTol) {
          // right going
          o_hUpdateRight += fWaves[waveNumber][0];
          o_huUpdateRight += fWaves[waveNumber][1];
        } else {
          // TODO: this case should not happen mathematically, but it does. Where is the bug? Machine accuracy only?
          //  MB: >>> waveSpeeds set to 0 for cases WetDryWall and DryWetWall
          o_hUpdateLeft += static_cast<real>(0.5) * fWaves[waveNumber][0];
          o_huUpdateLeft += static_cast<real>(0.5) * fWaves[waveNumber][1];

          o_hUpdateRight += static_cast<real>(0.5) * fWaves[waveNumber][0];
          o_huUpdateRight += static_cast<real>(0.5) * fWaves[waveNumber][1];
        }
      }

      // compute maximum wave speed (-> CFL-condition)
      waveSpeeds[0] = std::abs(waveSpeeds[0]);
      waveSpeeds[1] = std::abs(waveSpeeds[1]);
      waveSpeeds[2] = std::abs(waveSpeeds[2]);

      o_maxWaveSpeed = std::max(std::max(waveSpeeds[0], waveSpeeds[1]), waveSpeeds[2]);
    }

    /**
     * Computes the middle state of the homogeneous Riemann-problem.
     *   -> (\cite[ch. 13]{leveque2002finite})
     *
     * @param i_hLeft height on the left side of the edge.
     * @param i_hRight height on the right side of the edge.
     * @param i_uLeft velocity on the left side of the edge.
     * @param i_uRight velocity on the right side of the edge.
     * @param i_maxNumberOfNewtonIterations maximum number of Newton iterations.
     */
    inline void computeMiddleState(
      const real&         i_hLeft,
      const real&         i_hRight,
      const real&         i_uLeft,
      const real&         i_uRight,
      const unsigned int& i_maxNumberOfNewtonIterations,
      real&               o_hMiddle,
      real                o_middleStateSpeeds[2]
    ) const {
      //!!!MB: currently called with
      // i_maxNumberOfNewtonIterations = maxNumberOfNewtonIterations for wall boundaries
      // -> only leads to ShockShock or RarefactionRarefaction -> can be simplified
      // i_maxNumberOfNewtonIterations = 1 (computation of middle states)
      // -> could be changed towards a simplified implementation
      // ==> consider to split method into two or more separate implementations?

      // set everything to zero
      o_hMiddle              = static_cast<real>(0);
      o_middleStateSpeeds[0] = static_cast<real>(0);
      o_middleStateSpeeds[1] = static_cast<real>(0);
      // will be computed later
      real sqrt_g_hMiddle = static_cast<real>(0);

      // compute local square roots
      //(not necessarily the same ones as in computeNetUpdates!)
      real l_sqrt_g_hRight = std::sqrt(g * i_hRight);
      real l_sqrt_g_hLeft  = std::sqrt(g * i_hLeft);

      // single rarefaction in the case of a wet/dry interface
      if (i_hLeft < dryTol) {
        // ===== riemannStructure = DrySingleRarefaction; =====
        o_middleStateSpeeds[1] = o_middleStateSpeeds[0] = i_uRight - static_cast<real>(2) * l_sqrt_g_hRight;
        return;
      } else if (i_hRight < dryTol) {
        // ===== riemannStructure = SingleRarefactionDry; =====
        o_middleStateSpeeds[0] = o_middleStateSpeeds[1] = i_uLeft + static_cast<real>(2) * l_sqrt_g_hLeft;
        return;
      }

      /***************************************************************************************
       * Determine the wave structure of the Riemann-problem (no dry cases)
       **************************************************************************************/
      const real hMin = std::min(i_hLeft, i_hRight);
      const real hMax = std::max(i_hLeft, i_hRight);

      const real uDif = i_uRight - i_uLeft;

      if (0 <= static_cast<real>(2) * (std::sqrt(g * hMin) - std::sqrt(g * hMax)) + uDif) {
        // ===== riemannStructure = RarefactionRarefaction; =====
        o_hMiddle = std::max(
          static_cast<real>(0), i_uLeft - i_uRight + static_cast<real>(2) * (l_sqrt_g_hLeft + l_sqrt_g_hRight)
        );
        o_hMiddle = o_hMiddle * o_hMiddle / (static_cast<real>(16) * g);

        sqrt_g_hMiddle = sqrt_g * std::sqrt(o_hMiddle);
        // =======================================================
      } else if ((hMax - hMin) * std::sqrt(half_g * (1 / hMax + 1 / hMin)) + uDif <= 0) {
        // ===== riemannStructure = ShockShock; ==================
        o_hMiddle = std::min(i_hLeft, i_hRight);

        real l_sqrtTermH[2] = {0, 0};

        for (unsigned int i = 0; i < i_maxNumberOfNewtonIterations; i++) {
          l_sqrtTermH[0] = std::sqrt(half_g * ((o_hMiddle + i_hLeft) / (o_hMiddle * i_hLeft)));
          l_sqrtTermH[1] = std::sqrt(half_g * ((o_hMiddle + i_hRight) / (o_hMiddle * i_hRight)));

          real phi = i_uRight - i_uLeft + (o_hMiddle - i_hLeft) * l_sqrtTermH[0]
                     + (o_hMiddle - i_hRight) * l_sqrtTermH[1];

          if (std::abs(phi) < newtonTol)
            break;

          real derivativePhi
            = l_sqrtTermH[0] + l_sqrtTermH[1]
              - static_cast<real>(0.25) * g
                  * ((o_hMiddle - i_hLeft) / (l_sqrtTermH[0] * o_hMiddle * o_hMiddle) + (o_hMiddle - i_hRight) / (l_sqrtTermH[1] * o_hMiddle * o_hMiddle));
          o_hMiddle = o_hMiddle - phi / derivativePhi; // Newton step
          assert(o_hMiddle >= dryTol);
        }

        sqrt_g_hMiddle = sqrt_g * std::sqrt(o_hMiddle);
        // =======================================================
      } else {
        // ===== riemannStructure = ShockRarefaction; ============
        // ===== riemannStructure = RarefactionShock; ============
        o_hMiddle = hMin;

        sqrt_g_hMiddle   = sqrt_g * std::sqrt(o_hMiddle);
        real sqrt_g_hMax = sqrt_g * std::sqrt(hMax);
        for (unsigned int i = 0; i < i_maxNumberOfNewtonIterations; i++) {
          real sqrtTermHMin = half_g * ((o_hMiddle + hMin) / (o_hMiddle * hMin));

          real phi = i_uRight - i_uLeft + (o_hMiddle - hMin) * sqrtTermHMin
                     + static_cast<real>(2) * (sqrt_g_hMiddle - sqrt_g_hMax);

          if (std::abs(phi) < newtonTol)
            break;

          real derivativePhi = sqrtTermHMin
                               - static_cast<real>(0.25) * g * (o_hMiddle - hMin)
                                   / (o_hMiddle * o_hMiddle * sqrtTermHMin)
                               + sqrt_g / sqrt_g_hMiddle;

          o_hMiddle = o_hMiddle - phi / derivativePhi; // Newton step

          sqrt_g_hMiddle = sqrt_g * std::sqrt(o_hMiddle);
        }
        // =======================================================
      }
      /***************************************************************************************
       * Determine the wave structure end
       **************************************************************************************/

      o_middleStateSpeeds[0] = i_uLeft + static_cast<real>(2) * l_sqrt_g_hLeft - static_cast<real>(3) * sqrt_g_hMiddle;
      o_middleStateSpeeds[1] = i_uRight - static_cast<real>(2) * l_sqrt_g_hRight + static_cast<real>(3) * sqrt_g_hMiddle;

      assert(o_hMiddle >= 0);
    }

  }; // end of class AugRieFun

} // end of namespace solver

#endif /* AUGRIE_FUN_HPP */
