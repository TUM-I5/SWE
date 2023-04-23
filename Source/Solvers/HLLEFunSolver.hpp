/**
 * HLLEFun.hpp
 * @file
 * This file is part of Pond.
 *
 ****
 **** HLLE Solver for the Shallow Water Equations
 **** vectorizable functional implementation (based on AugRieFun)
 ****
 ****
 *
 *  Author: Nicolai Schaffroth
 *
 *  AugRiefun: Alexander Breuer
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
 ***
 *
 * @section LICENSE
 *
 * Pond is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Pond is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Pond.  If not, see <http://www.gnu.org/licenses/>.
 *
 ****
 */

#ifndef HLLE_FUN_HPP
#define HLLE_FUN_HPP

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
  template <typename T>
  class HLLEFun {
  private:
    const T dryTol;
    const T g;      // gravity constant
    const T half_g; // 0.5 * gravity constant
    const T sqrt_g; // square root of the gravity constant
    const T zeroTol;

  public:
    /**
     * AugRieFun Constructor, takes three problem parameters
     * @param dryTol "dry tolerance": if the water height falls below dryTol, wall boundary conditions are applied
     * (default value is 0.01)
     * @param gravity takes the value of the gravity constant (default value is 9.81 m/s^2)
     * @param zeroTol computed f-waves with an absolute value < zeroTol are treated as static waves (default value is
     * 10^{-7})
     */
    HLLEFun(T i_dryTol = (T)0.01, T i_gravity = (T)9.81, T i_zeroTol = (T)0.0000001):
      dryTol(i_dryTol),
      g(i_gravity),
      half_g((T).5 * i_gravity),
      sqrt_g(std::sqrt(i_gravity)),
      zeroTol(i_zeroTol) {}

// Define MIN and MAX as Macros, std::min and std::max seem to crash icpc
// For some reason, the Intel C++ compiler crashes during when Vectorization is enabled.
#define MIN(x, y) (((x) <= (y)) ? (x) : (y))
#define MAX(x, y) (((x) >= (y)) ? (x) : (y))

/**
 * Compute net updates for the left/right cell of the edge.
 *
 * maxWaveSpeed will be set to the maximum (linearized) wave speed => CFL
 */
#ifdef VECTORIZE
#pragma omp declare simd
#endif
    void computeNetUpdates(
      const T i_hLeft,
      const T i_hRight,
      const T i_huLeft,
      const T i_huRight,
      const T i_bLeft,
      const T i_bRight,

      T& o_hUpdateLeft,
      T& o_hUpdateRight,
      T& o_huUpdateLeft,
      T& o_huUpdateRight,
      T& o_maxWaveSpeed
    ) const {
      T hLeft   = i_hLeft;
      T hRight  = i_hRight;
      T uLeft   = (T)0;
      T uRight  = (T)0;
      T huLeft  = i_huLeft;
      T huRight = i_huRight;
      T bLeft   = i_bLeft;
      T bRight  = i_bRight;

      // declare variables which are used over and over again
      T sqrt_g_hLeft;
      T sqrt_g_hRight;

      T sqrt_hLeft;
      T sqrt_hRight;

      // reset net updates and the maximum wave speed
      o_hUpdateLeft = o_hUpdateRight = o_huUpdateLeft = o_huUpdateRight = T(0);
      o_maxWaveSpeed                                                    = (T)0;

      // HLLE middle height (Unused?)
      // T hMiddle = (T)0;

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
        hLeft = huLeft = uLeft = (T)0;
      }

      if (hRight >= dryTol) {
        uRight = huRight / hRight;
      } else {
        bRight += hRight;
        hRight = huRight = uRight = (T)0;
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
            // Removed middle state computation -> assuming momentum is never high enough to overcome boundary on its
            // own
            hRight  = hLeft;
            uRight  = -uLeft;
            huRight = -huLeft;
            bRight = bLeft = (T)0;
            wetDryState    = WetDryWall;
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
            // Removed middle state computation -> assuming momentum is never high enough to overcome boundary on its
            // own
            hLeft  = hRight;
            uLeft  = -uRight;
            huLeft = -huRight;
            bLeft = bRight = (T)0;
            wetDryState    = DryWetWall;
          }
        } else { // hLeft < dryTol and hRight < dryTol
          wetDryState = DryDry;
// nothing to do for dry/dry case, all netUpdates and maxWaveSpeed are 0
#ifdef VECTORIZE // Vectorization fails with branching return statement. Let computation resume with dummy values
                 // instead.
          hLeft = hRight = dryTol;
          uLeft = uRight = (T)0;
          huLeft = huRight = (T)0;
          bLeft = bRight = (T)0;
#else
          return;
#endif
        }
      };

      /***************************************************************************************
       * Determine Wet Dry State End
       **************************************************************************************/
      // MB: not executed for DryDry => returns right after case is identified
      // assert(wetDryState != DryDry);

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
      const T characteristicSpeeds[2] = {uLeft - sqrt_g_hLeft, uRight + sqrt_g_hRight};

      // compute "Roe speeds"
      const T hRoe = (T)0.5 * (hRight + hLeft);
      const T uRoe = (uLeft * sqrt_hLeft + uRight * sqrt_hRight) / (sqrt_hLeft + sqrt_hRight);

      // optimization for dumb compilers
      const T sqrt_g_hRoe  = sqrt_g * std::sqrt(hRoe);
      const T roeSpeeds[2] = {uRoe - sqrt_g_hRoe, uRoe + sqrt_g_hRoe};

      // middle state computation removed, using basic Einfeldt speeds

      T extEinfeldtSpeeds[2] = {(T)0, (T)0};
      if (wetDryState == WetWet || wetDryState == WetDryWall || wetDryState == DryWetWall) {
        extEinfeldtSpeeds[0] = MIN(characteristicSpeeds[0], roeSpeeds[0]);

        extEinfeldtSpeeds[1] = MAX(characteristicSpeeds[1], roeSpeeds[1]);
      } else if (hLeft < dryTol) { // MB: !!! cases DryWetInundation, DryWetWallInundation
        // ignore undefined speeds
        extEinfeldtSpeeds[0] = roeSpeeds[0];
        extEinfeldtSpeeds[1] = MAX(characteristicSpeeds[1], roeSpeeds[1]);

      } else if (hRight < dryTol) { // MB: !!! cases WetDryInundation, WetDryWallInundation
        // ignore undefined speeds
        extEinfeldtSpeeds[0] = MIN(characteristicSpeeds[0], roeSpeeds[0]);
        extEinfeldtSpeeds[1] = roeSpeeds[1];

      } else {
// This case corresponds to the early return statement above. Should only be executed with vectorization enabled.
#ifdef VECTORIZE
        extEinfeldtSpeeds[0] = (T)0;
        extEinfeldtSpeeds[1] = (T)1;
#else
        assert(false);
#endif
      }

      // HLL middle state
      //   \cite[theorem 3.1]{george2006finite}, \cite[ch. 4.1]{george2008augmented}
      const T hLLMiddleHeight = MAX(
        (huLeft - huRight + extEinfeldtSpeeds[1] * hRight - extEinfeldtSpeeds[0] * hLeft
        ) / (extEinfeldtSpeeds[1] - extEinfeldtSpeeds[0]),
        (T)0
      );

      // define eigenvalues
      const T eigenValues[3] = {
        extEinfeldtSpeeds[0], (T)0.5 * (extEinfeldtSpeeds[0] + extEinfeldtSpeeds[1]), extEinfeldtSpeeds[1]};

      // define eigenvectors
      //  MB: no longer used as system matrix
      //      -> but still used to compute f-waves
      const T eigenVectors[3][3] = {
        {(T)1, (T)0, (T)1},
        {eigenValues[0], (T)0, eigenValues[2]},
        {eigenValues[0] * eigenValues[0], (T)1, eigenValues[2] * eigenValues[2]}};

      // compute the jump in state
      T rightHandSide[3] = {
        hRight - hLeft,
        huRight - huLeft,
        (huRight * uRight + half_g * hRight * hRight) - (huLeft * uLeft + half_g * hLeft * hLeft)};

      // compute steady state wave
      //   \cite[ch. 4.2.1 \& app. A]{george2008augmented}, \cite[ch. 6.2 \& ch. 4.4]{george2006finite}
      T steadyStateWave[2] = {-(bRight - bLeft), -half_g * (hLeft + hRight) * (bRight - bLeft)};

      // preserve depth-positivity
      //   \cite[ch. 6.5.2]{george2006finite}, \cite[ch. 4.2.3]{george2008augmented}
      if (eigenValues[0] < -zeroTol && eigenValues[2] > zeroTol) {
        // subsonic
        steadyStateWave[0] = MAX(
          steadyStateWave[0], hLLMiddleHeight * (eigenValues[2] - eigenValues[0]) / eigenValues[0]
        );
        steadyStateWave[0] = MIN(
          steadyStateWave[0], hLLMiddleHeight * (eigenValues[2] - eigenValues[0]) / eigenValues[2]
        );
      } else if (eigenValues[0] > zeroTol) {
        // supersonic right TODO: motivation?
        steadyStateWave[0] = MAX(steadyStateWave[0], -hLeft);
        steadyStateWave[0] = MIN(
          steadyStateWave[0], hLLMiddleHeight * (eigenValues[2] - eigenValues[0]) / eigenValues[0]
        );
      } else if (eigenValues[2] < -zeroTol) {
        // supersonic left TODO: motivation?
        steadyStateWave[0] = MAX(
          steadyStateWave[0], hLLMiddleHeight * (eigenValues[2] - eigenValues[0]) / eigenValues[2]
        );
        steadyStateWave[0] = MIN(steadyStateWave[0], hRight);
      }
      // Limit the effect of the source term
      //   \cite[ch. 6.4.2]{george2006finite}
      steadyStateWave[1] = MIN(steadyStateWave[1], g * MAX(-hLeft * (bRight - bLeft), -hRight * (bRight - bLeft)));
      steadyStateWave[1] = MAX(steadyStateWave[1], g * MIN(-hLeft * (bRight - bLeft), -hRight * (bRight - bLeft)));

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
      T inverseDiff = (T)1 / (eigenValues[2] - eigenValues[0]); // 2 FLOPs (1 div)
      // compute f-waves:
      T beta[3];
      beta[0] = (eigenValues[2] * rightHandSide[0] - rightHandSide[1]) * inverseDiff;  // 3 FLOPs
      beta[2] = (-eigenValues[0] * rightHandSide[0] + rightHandSide[1]) * inverseDiff; // 3 FLOPs

      // step 2: solve 3rd row for beta[1]:
      beta[1] = rightHandSide[2] - eigenValues[0] * eigenValues[0] * beta[0] - eigenValues[2] * eigenValues[2] * beta[2];

#ifdef false
      // compute inverse of 3x3 matrix
      // MB: !!!AVOID COMPUTATION of inverse => see above for better solution for simple matrix
      // MB: TODO -> different choice of matrix in AugRie.hpp -> requires a separate solution procedure
      const T m[3][3] = {
        {(eigenVectors[1][1] * eigenVectors[2][2] - eigenVectors[1][2] * eigenVectors[2][1]),
         -(eigenVectors[0][1] * eigenVectors[2][2] - eigenVectors[0][2] * eigenVectors[2][1]),
         (eigenVectors[0][1] * eigenVectors[1][2] - eigenVectors[0][2] * eigenVectors[1][1])},
        {-(eigenVectors[1][0] * eigenVectors[2][2] - eigenVectors[1][2] * eigenVectors[2][0]),
         (eigenVectors[0][0] * eigenVectors[2][2] - eigenVectors[0][2] * eigenVectors[2][0]),
         -(eigenVectors[0][0] * eigenVectors[1][2] - eigenVectors[0][2] * eigenVectors[1][0])},
        {(eigenVectors[1][0] * eigenVectors[2][1] - eigenVectors[1][1] * eigenVectors[2][0]),
         -(eigenVectors[0][0] * eigenVectors[2][1] - eigenVectors[0][1] * eigenVectors[2][0]),
         (eigenVectors[0][0] * eigenVectors[1][1] - eigenVectors[0][1] * eigenVectors[1][0])}};
      const T d = (eigenVectors[0][0] * m[0][0] + eigenVectors[0][1] * m[1][0] + eigenVectors[0][2] * m[2][0]);

      // m stores not really the inverse matrix, but the inverse multiplied by d
      const T s = 1 / d;

      // compute m*rightHandSide
      const T beta[3] = {
        (m[0][0] * rightHandSide[0] + m[0][1] * rightHandSide[1] + m[0][2] * rightHandSide[2]) * s,
        (m[1][0] * rightHandSide[0] + m[1][1] * rightHandSide[1] + m[1][2] * rightHandSide[2]) * s,
        (m[2][0] * rightHandSide[0] + m[2][1] * rightHandSide[1] + m[2][2] * rightHandSide[2]) * s};
#endif

      /***************************************************************************************
       * Solve linear equation end
       **************************************************************************************/

      // compute f-waves and wave-speeds
      T fWaves[3][2];  // array to store the three f-wave vectors (2 components each)
      T waveSpeeds[3]; // and vector to their speeds

      if (wetDryState == WetDryWall) {
        // zero ghost updates (wall boundary)
        // care about the left going wave (0) only
        fWaves[0][0] = beta[0] * eigenVectors[1][0];
        fWaves[0][1] = beta[0] * eigenVectors[2][0];

        // set the rest to zero
        fWaves[1][0] = fWaves[1][1] = (T)0;
        fWaves[2][0] = fWaves[2][1] = (T)0;

        waveSpeeds[0] = eigenValues[0];
        waveSpeeds[1] = waveSpeeds[2] = (T)0;

        assert(eigenValues[0] < zeroTol);
      } else if (wetDryState == DryWetWall) {
        // zero ghost updates (wall boundary)
        // care about the right going wave (2) only
        fWaves[2][0] = beta[2] * eigenVectors[1][2];
        fWaves[2][1] = beta[2] * eigenVectors[2][2];

        // set the rest to zero
        fWaves[0][0] = fWaves[0][1] = (T)0;
        fWaves[1][0] = fWaves[1][1] = (T)0;

        waveSpeeds[2] = eigenValues[2];
        waveSpeeds[0] = waveSpeeds[1] = (T)0;

        assert(eigenValues[2] > -zeroTol);
      } else {
        // compute f-waves (default)
        // loop manually unrolled for vectorization testing
        /*
        for (int waveNumber = 0; waveNumber < 3; waveNumber++) {
            fWaves[waveNumber][0] = beta[waveNumber] * eigenVectors[1][waveNumber]; //select 2nd and
            fWaves[waveNumber][1] = beta[waveNumber] * eigenVectors[2][waveNumber]; //3rd component of the augmented
        decomposition
        }
        */

        fWaves[0][0] = beta[0] * eigenVectors[1][0];
        fWaves[0][1] = beta[0] * eigenVectors[2][0];

        fWaves[1][0] = beta[1] * eigenVectors[1][1];
        fWaves[1][1] = beta[1] * eigenVectors[2][1];

        fWaves[2][0] = beta[2] * eigenVectors[1][2];
        fWaves[2][1] = beta[2] * eigenVectors[2][2];

        waveSpeeds[0] = eigenValues[0];
        waveSpeeds[1] = eigenValues[1];
        waveSpeeds[2] = eigenValues[2];
      }
      /***************************************************************************************
       * Compute Wave Decomposition End
       **************************************************************************************/

#ifdef VECTORIZE
      // DryDry state with dummy data only reaches this point with vectorization enabled.
      // Return values are not updated in that case, since nothing happens in those cells.
      if (wetDryState != DryDry) {
#endif
        // compute the updates from the three propagating waves
        // A^-\delta Q = \sum{s[i]<0} \beta[i] * r[i] = A^-\delta Q = \sum{s[i]<0} Z^i
        // A^+\delta Q = \sum{s[i]>0} \beta[i] * r[i] = A^-\delta Q = \sum{s[i]<0} Z^i

        // loop manually unrolled for vectorization testing
        /*
        for (int waveNumber = 0; waveNumber < 3; waveNumber++) {
            if (waveSpeeds[waveNumber] < -zeroTol) {
                //left going
                o_hUpdateLeft  += fWaves[waveNumber][0];
                o_huUpdateLeft += fWaves[waveNumber][1];
            } else if (waveSpeeds[waveNumber] > zeroTol) {
                //right going
                o_hUpdateRight  += fWaves[waveNumber][0];
                o_huUpdateRight += fWaves[waveNumber][1];
            } else {
                //TODO: this case should not happen mathematically, but it does. Where is the bug? Machine accuracy
        only?
                // MB: >>> waveSpeeds set to 0 for cases WetDryWall and DryWetWall
                o_hUpdateLeft  += (T)0.5 * fWaves[waveNumber][0];
                o_huUpdateLeft += (T)0.5 * fWaves[waveNumber][1];

                o_hUpdateRight  += (T)0.5 * fWaves[waveNumber][0];
                o_huUpdateRight += (T)0.5 * fWaves[waveNumber][1];
            }
        }
        */

        if (waveSpeeds[0] < -zeroTol) {
          // left going
          o_hUpdateLeft += fWaves[0][0];
          o_huUpdateLeft += fWaves[0][1];
        } else if (waveSpeeds[0] > zeroTol) {
          // right going
          o_hUpdateRight += fWaves[0][0];
          o_huUpdateRight += fWaves[0][1];
        } else {
          // TODO: this case should not happen mathematically, but it does. Where is the bug? Machine accuracy only?
          //  MB: >>> waveSpeeds set to 0 for cases WetDryWall and DryWetWall
          o_hUpdateLeft += (T)0.5 * fWaves[0][0];
          o_huUpdateLeft += (T)0.5 * fWaves[0][1];

          o_hUpdateRight += (T)0.5 * fWaves[0][0];
          o_huUpdateRight += (T)0.5 * fWaves[0][1];
        }
        if (waveSpeeds[1] < -zeroTol) {
          // left going
          o_hUpdateLeft += fWaves[1][0];
          o_huUpdateLeft += fWaves[1][1];
        } else if (waveSpeeds[1] > zeroTol) {
          // right going
          o_hUpdateRight += fWaves[1][0];
          o_huUpdateRight += fWaves[1][1];
        } else {
          // TODO: this case should not happen mathematically, but it does. Where is the bug? Machine accuracy only?
          //  MB: >>> waveSpeeds set to 0 for cases WetDryWall and DryWetWall
          o_hUpdateLeft += (T)0.5 * fWaves[1][0];
          o_huUpdateLeft += (T)0.5 * fWaves[1][1];

          o_hUpdateRight += (T)0.5 * fWaves[1][0];
          o_huUpdateRight += (T)0.5 * fWaves[1][1];
        }
        if (waveSpeeds[2] < -zeroTol) {
          // left going
          o_hUpdateLeft += fWaves[2][0];
          o_huUpdateLeft += fWaves[2][1];
        } else if (waveSpeeds[2] > zeroTol) {
          // right going
          o_hUpdateRight += fWaves[2][0];
          o_huUpdateRight += fWaves[2][1];
        } else {
          // TODO: this case should not happen mathematically, but it does. Where is the bug? Machine accuracy only?
          //  MB: >>> waveSpeeds set to 0 for cases WetDryWall and DryWetWall
          o_hUpdateLeft += (T)0.5 * fWaves[2][0];
          o_huUpdateLeft += (T)0.5 * fWaves[2][1];

          o_hUpdateRight += (T)0.5 * fWaves[2][0];
          o_huUpdateRight += (T)0.5 * fWaves[2][1];
        }

        // compute maximum wave speed (-> CFL-condition)
        waveSpeeds[0] = std::abs(waveSpeeds[0]);
        waveSpeeds[1] = std::abs(waveSpeeds[1]);
        waveSpeeds[2] = std::abs(waveSpeeds[2]);

        o_maxWaveSpeed = MAX(MAX(waveSpeeds[0], waveSpeeds[1]), waveSpeeds[2]);
#ifdef VECTORIZE
      }
#endif
    }

  }; // end of class HLLEFun

} // end of namespace solver

#endif /* HLLE_FUN_HPP */
