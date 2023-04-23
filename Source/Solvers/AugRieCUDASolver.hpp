/**
 * AugRie_CUDA.h
 *
 ****
 **** Approximate Augmented Riemann Solver for the Shallow Water Equations
 **** vectorized using SIMD-Extensions
 ****
 *
 *  Created on: May 28, 2013
 *  Last Update: Jul 14, 2013
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

#ifndef AUGRIE_CUDA_HPP_
#define AUGRIE_CUDA_HPP_

#define FLOAT32

#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "SIMD_TYPES.hpp"

const integer DryDry               = 0;
const integer WetWet               = 1;
const integer WetDryInundation     = 2;
const integer WetDryWall           = 3;
const integer WetDryWallInundation = 4;
const integer DryWetInundation     = 5;
const integer DryWetWall           = 6;
const integer DryWetWallInundation = 7;

const integer DrySingleRarefaction   = 0;
const integer SingleRarefactionDry   = 1;
const integer ShockShock             = 2;
const integer ShockRarefaction       = 3;
const integer RarefactionShock       = 4;
const integer RarefactionRarefaction = 5;

__device__ inline void computeMiddleState(
  const real&         i_hLeft,
  const real&         i_hRight,
  const real&         i_uLeft,
  const real&         i_uRight,
  const real          dryTol,
  const real          newtonTol,
  const real          g,
  const real          sqrt_g,
  const unsigned int& i_maxNumberOfNewtonIterations,
  real&               o_hMiddle,
  real                o_middleStateSpeeds[2]
);

/**
 * Compute net updates for the left/right cell of the edge.
 *
 * maxWaveSpeed will be set to the maximum (linearized) wave speed => CFL
 */
__device__ void augRieComputeNetUpdates(
  const real i_hLeft,
  const real i_hRight,
  const real i_huLeft,
  const real i_huRight,
  const real i_bLeft,
  const real i_bRight,

  const real         g,
  const real         dryTol,
  const real         newtonTol,
  const real         zeroTol,
  const unsigned int maxNumberOfNewtonIterations,

  real o_netUpdates[5]
) {
  real hLeft   = i_hLeft;
  real hRight  = i_hRight;
  real uLeft   = static_cast<real>(0);
  real uRight  = static_cast<real>(0);
  real huLeft  = i_huLeft;
  real huRight = i_huRight;
  real bLeft   = i_bLeft;
  real bRight  = i_bRight;

  // declare variables which are used over and over again
  const real sqrt_g = sqrtf(g);
  real       sqrt_g_hLeft;
  real       sqrt_g_hRight;

  real sqrt_hLeft;
  real sqrt_hRight;

  // set speeds to zero (will be determined later)
  uLeft = uRight = 0.;

  // reset net updates and the maximum wave speed
  o_netUpdates[0] = o_netUpdates[1] = o_netUpdates[2] = o_netUpdates[3] = static_cast<real>(0);
  o_netUpdates[4]                                                       = static_cast<real>(0);

  real hMiddle              = static_cast<real>(0);
  real middleStateSpeeds[2] = {static_cast<real>(0)};

  // determine the wet/dry state and compute local variables correspondingly

  /***************************************************************************************
   * Determine Wet Dry State Begin
   **************************************************************************************/
  integer wetDryState;
  // compute speeds or set them to zero (dry cells)
  if (hLeft > dryTol) {
    uLeft = huLeft / hLeft;
  } else {
    bLeft += hLeft;
    hLeft = huLeft = uLeft = 0;
  }

  if (hRight > dryTol) {
    uRight = huRight / hRight;
  } else {
    bRight += hRight;
    hRight = huRight = uRight = 0;
  }

  if (hLeft >= dryTol and hRight >= dryTol) {
    // test for simple wet/wet case since this is most probably the
    // most frequently executed branch
    wetDryState = WetWet;
  } else if (hLeft < dryTol and hRight < dryTol) {
    // check for the dry/dry-case
    wetDryState = DryDry;
  } else if (hLeft < dryTol and hRight + bRight > bLeft) {
    // we have a shoreline: one cell dry, one cell wet

    // check for simple inundation problems
    //  (=> dry cell lies lower than the wet cell)
    wetDryState = DryWetInundation;
  } else if (hRight < dryTol and hLeft + bLeft > bRight) {
    wetDryState = WetDryInundation;
  } else if (hLeft < dryTol) {
    // dry cell lies higher than the wet cell
    // lets check if the momentum is able to overcome the difference in height
    //   => solve homogeneous Riemann-problem to determine the middle state height
    //      which would arise if there is a wall (wall-boundary-condition)
    //        \cite[ch. 6.8.2]{george2006finite})
    //        \cite[ch. 5.2]{george2008augmented}
    computeMiddleState(
      hRight,
      hRight,
      -uRight,
      uRight,
      dryTol,
      newtonTol,
      g,
      sqrt_g,
      maxNumberOfNewtonIterations,
      hMiddle,
      middleStateSpeeds
    );

    if (hMiddle + bRight > bLeft) {
      // momentum is large enough, continue with the original values
      //           bLeft = o_hMiddle + bRight;
      wetDryState = DryWetWallInundation;
    } else {
      // momentum is not large enough, use wall-boundary-values
      hLeft  = hRight;
      uLeft  = -uRight;
      huLeft = -huRight;
      bLeft = bRight = static_cast<real>(0);
      wetDryState    = DryWetWall;
    }
  } else if (hRight < dryTol) {
    // lets check if the momentum is able to overcome the difference in height
    //   => solve homogeneous Riemann-problem to determine the middle state height
    //      which would arise if there is a wall (wall-boundary-condition)
    //        \cite[ch. 6.8.2]{george2006finite})
    //        \cite[ch. 5.2]{george2008augmented}
    computeMiddleState(
      hLeft, hLeft, uLeft, -uLeft, dryTol, newtonTol, g, sqrt_g, maxNumberOfNewtonIterations, hMiddle, middleStateSpeeds
    );

    if (hMiddle + bLeft > bRight) {
      // momentum is large enough, continue with the original values
      //           bRight = o_hMiddle + bLeft;
      wetDryState = WetDryWallInundation;
    } else {
      hRight  = hLeft;
      uRight  = -uLeft;
      huRight = -huLeft;
      bRight = bLeft = static_cast<real>(0);
      wetDryState    = WetDryWall;
    }
  } else {
    // done with all cases
    assert(false);
  }

  // limit the effect of the source term if there is a "wall"
  //\cite[end of ch. 5.2?]{george2008augmented}
  //\cite[rpn2ez_fast_geo.f][levequeclawpack]
  if (wetDryState == DryWetWallInundation) {
    bLeft = hRight + bRight;
  } else if (wetDryState == WetDryWallInundation) {
    bRight = hLeft + bLeft;
  }
  /***************************************************************************************
   * Determine Wet Dry State End
   **************************************************************************************/

  if (wetDryState != DryDry) {
    // precompute some terms which are fixed during
    // the computation after some specific point
    sqrt_hLeft  = std::sqrt(hLeft);
    sqrt_hRight = std::sqrt(hRight);

    sqrt_g_hLeft  = sqrt_g * sqrt_hLeft;
    sqrt_g_hRight = sqrt_g * sqrt_hRight;

    // where to store the three waves
    real fWaves[3][2];
    // and their speeds
    real waveSpeeds[3];

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
    const real sqrt_g_hRoe  = std::sqrt(g * hRoe);
    const real roeSpeeds[2] = {uRoe - sqrt_g_hRoe, uRoe + sqrt_g_hRoe};

    // compute the middle state of the homogeneous Riemann-Problem
    if (wetDryState != WetDryWall and wetDryState != DryWetWall) {
      // case WDW and DWW was computed in determineWetDryState already
      computeMiddleState(hLeft, hRight, uLeft, uRight, dryTol, newtonTol, g, sqrt_g, 1, hMiddle, middleStateSpeeds);
    }

    // compute extended eindfeldt speeds (einfeldt speeds + middle state speeds)
    //   \cite[ch. 5.2]{george2008augmented}, \cite[ch. 6.8]{george2006finite}
    real extEinfeldtSpeeds[2] = {static_cast<real>(0), static_cast<real>(0)};
    if (wetDryState == WetWet or wetDryState == WetDryWall or wetDryState == DryWetWall) {
      extEinfeldtSpeeds[0] = fmin(characteristicSpeeds[0], roeSpeeds[0]);
      extEinfeldtSpeeds[0] = fmin(extEinfeldtSpeeds[0], middleStateSpeeds[1]);

      extEinfeldtSpeeds[1] = fmax(characteristicSpeeds[1], roeSpeeds[1]);
      extEinfeldtSpeeds[1] = fmax(extEinfeldtSpeeds[1], middleStateSpeeds[0]);
    } else if (hLeft < dryTol) {
      // ignore undefined speeds
      extEinfeldtSpeeds[0] = fmin(roeSpeeds[0], middleStateSpeeds[1]);
      extEinfeldtSpeeds[1] = fmax(characteristicSpeeds[1], roeSpeeds[1]);

      assert(middleStateSpeeds[0] < extEinfeldtSpeeds[1]);
    } else if (hRight < dryTol) {
      // ignore undefined speeds
      extEinfeldtSpeeds[0] = fmin(characteristicSpeeds[0], roeSpeeds[0]);
      extEinfeldtSpeeds[1] = fmax(roeSpeeds[1], middleStateSpeeds[0]);

      assert(middleStateSpeeds[1] > extEinfeldtSpeeds[0]);
    } else {
      assert(false);
    }

    // HLL middle state
    //   \cite[theorem 3.1]{george2006finite}, \cite[ch. 4.1]{george2008augmented}
    const real hLLMiddleHeight = fmax(
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
      huRight * uRight + static_cast<real>(0.5) * g * hRight * hRight
        - (huLeft * uLeft + static_cast<real>(0.5) * g * hLeft * hLeft)};

    // compute steady state wave
    //   \cite[ch. 4.2.1 \& app. A]{george2008augmented}, \cite[ch. 6.2 \& ch. 4.4]{george2006finite}
    const real hBar = (hLeft + hRight) * static_cast<real>(0.5);

    real steadyStateWave[2] = {-(bRight - bLeft), -g * hBar * (bRight - bLeft)};

    // preserve depth-positivity
    //   \cite[ch. 6.5.2]{george2006finite}, \cite[ch. 4.2.3]{george2008augmented}
    if (eigenValues[0] < -zeroTol and eigenValues[2] > zeroTol) {
      // subsonic
      steadyStateWave[0] = fmax(
        steadyStateWave[0], hLLMiddleHeight * (eigenValues[2] - eigenValues[0]) / eigenValues[0]
      );
      steadyStateWave[0] = fmin(
        steadyStateWave[0], hLLMiddleHeight * (eigenValues[2] - eigenValues[0]) / eigenValues[2]
      );
    } else if (eigenValues[0] > zeroTol) {
      // supersonic right TODO: motivation?
      steadyStateWave[0] = fmax(steadyStateWave[0], -hLeft);
      steadyStateWave[0] = fmin(
        steadyStateWave[0], hLLMiddleHeight * (eigenValues[2] - eigenValues[0]) / eigenValues[0]
      );
    } else if (eigenValues[2] < -zeroTol) {
      // supersonic left TODO: motivation?
      steadyStateWave[0] = fmax(
        steadyStateWave[0], hLLMiddleHeight * (eigenValues[2] - eigenValues[0]) / eigenValues[2]
      );
      steadyStateWave[0] = fmin(steadyStateWave[0], hRight);
    }

    // Limit the effect of the source term
    //   \cite[ch. 6.4.2]{george2006finite}
    steadyStateWave[1] = fmin(steadyStateWave[1], g * fmax(-hLeft * (bRight - bLeft), -hRight * (bRight - bLeft)));
    steadyStateWave[1] = fmax(steadyStateWave[1], g * fmin(-hLeft * (bRight - bLeft), -hRight * (bRight - bLeft)));

    rightHandSide[0] -= steadyStateWave[0];
    // rightHandSide[1]: no source term
    rightHandSide[2] -= steadyStateWave[1];

    // everything is ready, solve the equations!
    /***************************************************************************************
     * Solve linear equation begin
     **************************************************************************************/
    // compute inverse of 3x3 matrix
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
    /***************************************************************************************
     * Solve linear equation end
     **************************************************************************************/

    // compute f-waves and wave-speeds
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
        o_netUpdates[0] += fWaves[waveNumber][0];
        o_netUpdates[2] += fWaves[waveNumber][1];
      } else if (waveSpeeds[waveNumber] > zeroTol) {
        // right going
        o_netUpdates[1] += fWaves[waveNumber][0];
        o_netUpdates[3] += fWaves[waveNumber][1];
      } else {
        // TODO: this case should not happen mathematically, but it does. Where is the bug? Machine accuracy only?
        o_netUpdates[0] += static_cast<real>(0.5) * fWaves[waveNumber][0];
        o_netUpdates[2] += static_cast<real>(0.5) * fWaves[waveNumber][1];

        o_netUpdates[1] += static_cast<real>(0.5) * fWaves[waveNumber][0];
        o_netUpdates[3] += static_cast<real>(0.5) * fWaves[waveNumber][1];
      }
    }

    // compute maximum wave speed (-> CFL-condition)
    waveSpeeds[0] = fabs(waveSpeeds[0]);
    waveSpeeds[1] = fabs(waveSpeeds[1]);
    waveSpeeds[2] = fabs(waveSpeeds[2]);

    o_netUpdates[4] = fmax(waveSpeeds[0], waveSpeeds[1]);
    o_netUpdates[4] = fmax(o_netUpdates[4], waveSpeeds[2]);
  }
}

/**
 * Computes the middle state of the homogeneous Riemann-problem.
 *   -> (\cite[ch. 13]{leveque2002finite})
 *
 * @param i_hLeft height on the left side of the edge.
 * @param i_hRight height on the right side of the edge.
 * @param i_uLeft velocity on the left side of the edge.
 * @param i_uRight velocity on the right side of the edge.
 * @param i_huLeft momentum on the left side of the edge.
 * @param i_huRight momentum on the right side of the edge.
 * @param i_maxNumberOfNewtonIterations maximum number of Newton iterations.
 */
__device__ inline void computeMiddleState(
  const real&         i_hLeft,
  const real&         i_hRight,
  const real&         i_uLeft,
  const real&         i_uRight,
  const real          dryTol,
  const real          newtonTol,
  const real          g,
  const real          sqrt_g,
  const unsigned int& i_maxNumberOfNewtonIterations,
  real&               o_hMiddle,
  real                o_middleStateSpeeds[2]
) {
  // set everything to zero
  o_hMiddle              = static_cast<real>(0);
  o_middleStateSpeeds[0] = static_cast<real>(0);
  o_middleStateSpeeds[1] = static_cast<real>(0);

  // compute local square roots
  //(not necessarily the same ones as in computeNetUpdates!)
  real l_sqrt_g_hRight = std::sqrt(g * i_hRight);
  real l_sqrt_g_hLeft  = std::sqrt(g * i_hLeft);

  // single rarefaction in the case of a wet/dry interface
  integer riemannStructure;
  if (i_hLeft < dryTol) {
    o_middleStateSpeeds[1] = o_middleStateSpeeds[0] = i_uRight - static_cast<real>(2) * l_sqrt_g_hRight;
    riemannStructure                                = DrySingleRarefaction;
    return;
  } else if (i_hRight < dryTol) {
    o_middleStateSpeeds[0] = o_middleStateSpeeds[1] = i_uLeft + static_cast<real>(2) * l_sqrt_g_hLeft;
    riemannStructure                                = SingleRarefactionDry;
    return;
  }

  // determine the wave structure of the Riemann-problem
  /***************************************************************************************
   * Determine riemann structure begin
   **************************************************************************************/
  const real hMin = fmin(i_hLeft, i_hRight);
  const real hMax = fmax(i_hLeft, i_hRight);

  const real uDif = i_uRight - i_uLeft;

  if (0 <= static_cast<real>(2) * (std::sqrt(g * hMin) - std::sqrt(g * hMax)) + uDif) {
    riemannStructure = RarefactionRarefaction;
  } else if ((hMax - hMin) * std::sqrt(g * static_cast<real>(0.5) * (1 / hMax + 1 / hMin)) + uDif <= 0) {
    riemannStructure = ShockShock;
  } else if (i_hLeft < i_hRight) {
    riemannStructure = ShockRarefaction;
  } else {
    riemannStructure = RarefactionShock;
  }
  /***************************************************************************************
   * Determine riemann structure end
   **************************************************************************************/

  // will be computed later
  real sqrt_g_hMiddle = static_cast<real>(0);

  if (riemannStructure == ShockShock) {
    o_hMiddle = fmin(i_hLeft, i_hRight);

    real l_sqrtTermH[2] = {0, 0};

    for (unsigned int i = 0; i < i_maxNumberOfNewtonIterations; i++) {
      l_sqrtTermH[0] = std::sqrt(static_cast<real>(0.5) * g * ((o_hMiddle + i_hLeft) / (o_hMiddle * i_hLeft)));
      l_sqrtTermH[1] = std::sqrt(static_cast<real>(0.5) * g * ((o_hMiddle + i_hRight) / (o_hMiddle * i_hRight)));

      real phi = i_uRight - i_uLeft + (o_hMiddle - i_hLeft) * l_sqrtTermH[0] + (o_hMiddle - i_hRight) * l_sqrtTermH[1];

      if (std::fabs(phi) < newtonTol) {
        break;
      }

      real derivativePhi = l_sqrtTermH[0] + l_sqrtTermH[1]
                           - static_cast<real>(0.25) * g * (o_hMiddle - i_hLeft)
                               / (l_sqrtTermH[0] * o_hMiddle * o_hMiddle)
                           - static_cast<real>(0.25) * g * (o_hMiddle - i_hRight)
                               / (l_sqrtTermH[1] * o_hMiddle * o_hMiddle);

      o_hMiddle = o_hMiddle - phi / derivativePhi; // Newton step
      assert(o_hMiddle >= dryTol);
    }

    sqrt_g_hMiddle = std::sqrt(g * o_hMiddle);
  }

  if (riemannStructure == RarefactionRarefaction) {
    o_hMiddle = fmax(
      static_cast<real>(0), i_uLeft - i_uRight + static_cast<real>(2) * (l_sqrt_g_hLeft + l_sqrt_g_hRight)
    );
    o_hMiddle = static_cast<real>(1) / (static_cast<real>(16) * g) * o_hMiddle * o_hMiddle;

    sqrt_g_hMiddle = std::sqrt(g * o_hMiddle);
  }

  if (riemannStructure == ShockRarefaction or riemannStructure == RarefactionShock) {
    real hMin, hMax;
    if (riemannStructure == ShockRarefaction) {
      hMin = i_hLeft;
      hMax = i_hRight;
    } else {
      hMin = i_hRight;
      hMax = i_hLeft;
    }

    o_hMiddle = hMin;

    sqrt_g_hMiddle   = std::sqrt(g * o_hMiddle);
    real sqrt_g_hMax = std::sqrt(g * hMax);
    for (unsigned int i = 0; i < i_maxNumberOfNewtonIterations; i++) {
      real sqrtTermHMin = std::sqrt(static_cast<real>(0.5) * g * ((o_hMiddle + hMin) / (o_hMiddle * hMin)));

      real phi = i_uRight - i_uLeft + (o_hMiddle - hMin) * sqrtTermHMin
                 + static_cast<real>(2) * (sqrt_g_hMiddle - sqrt_g_hMax);

      if (std::fabs(phi) < newtonTol) {
        break;
      }

      real derivativePhi = sqrtTermHMin
                           - static_cast<real>(0.25) * g * (o_hMiddle - hMin) / (o_hMiddle * o_hMiddle * sqrtTermHMin)
                           + sqrt_g / sqrt_g_hMiddle;

      o_hMiddle = o_hMiddle - phi / derivativePhi; // Newton step

      sqrt_g_hMiddle = std::sqrt(g * o_hMiddle);
    }
  }

  o_middleStateSpeeds[0] = i_uLeft + static_cast<real>(2) * l_sqrt_g_hLeft - static_cast<real>(3) * sqrt_g_hMiddle;
  o_middleStateSpeeds[1] = i_uRight - static_cast<real>(2) * l_sqrt_g_hRight + static_cast<real>(3) * sqrt_g_hMiddle;

  assert(o_hMiddle >= 0);
}

#endif /* AUGRIE_CUDA_HPP_ */
