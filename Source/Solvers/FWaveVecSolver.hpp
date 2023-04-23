/**
 * FWaveVecSolver.hpp
 *
 ****
 **** This is a vectorizable C++ implementation of the F-Wave solver (FWaveSolver.hpp).
 ****
 *
 * Created on: Nov 13, 2012
 * Last Update: Dec 28, 2013
 *
 ****
 *
 *  Author: Sebastian Rettenberger
 *    Homepage: http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.
 *    E-Mail: rettenbs AT in.tum.de
 *  Some optimzations: Michael Bader
 *    Homepage: http://www5.in.tum.de/wiki/index.php/Michael_Bader
 *    E-Mail: bader AT in.tum.de
 *
 ****
 *
 * (Main) Literature:
 *
 *   @article{bale2002wave,
 *            title={A wave propagation method for conservation laws and balance laws with spatially varying flux
 *functions}, author={Bale, D.S. and LeVeque, R.J. and Mitran, S. and Rossmanith, J.A.}, journal={SIAM Journal on
 *Scientific Computing}, volume={24}, number={3}, pages={955--978}, year={2002}}
 *
 *   @book{leveque2002finite,
 *         Author = {LeVeque, R. J.},
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
 */

#pragma once

#include <cmath>

namespace Solvers
{

    template <class T>
    class FWaveVecSolver
    {
      private:
        const T dryTol_;
        const T zeroTol_;
        const T halfGravity_;
        const T sqrtGravity_;

      public:
        /**
         * FWaveVec Constructor, takes three problem parameters
         * @param dryTol "dry tolerance": if the water height falls below dryTol, wall boundary conditions are applied
         * (default value is 100)
         * @param gravity takes the value of the gravity constant (default value is 9.81 m/s^2)
         * @param zeroTol computed f-waves with an absolute value < zeroTol are treated as static waves (default value
         * is 10^{-7})
         */
        FWaveVecSolver(T dryTol = T(1.0), T gravity = T(9.81), T zeroTol = T(0.0000001)):
            dryTol_(dryTol),
            zeroTol_(zeroTol),
            halfGravity_(T(0.5) * gravity),
            sqrtGravity_(std::sqrt(gravity))
        {
        }

        /**
         * Takes the water height, discharge and bathymatry in the left and right cell
         * and computes net updates (left and right going waves) according to the f-wave approach.
         * It also returns the maximum wave speed.
         */
#ifdef ENABLE_VECTORIZATION
#pragma omp declare simd
#endif
        void computeNetUpdates(
            T  hLeft,
            T  hRight,
            T  huLeft,
            T  huRight,
            T  bLeft,
            T  bRight,
            T& o_hUpdateLeft,
            T& o_hUpdateRight,
            T& o_huUpdateLeft,
            T& o_huUpdateRight,
            T& o_maxWaveSpeed
        ) const
        {
            if (hLeft >= dryTol_)
            {
                if (hRight < dryTol_)
                {
                    // Dry/Wet case
                    // Set values according to wall boundary condition
                    hRight  = hLeft;
                    huRight = -huLeft;
                    bRight  = bLeft;
                }
            }
            else if (hRight >= dryTol_)
            {
                // Wet/Dry case
                // Set values according to wall boundary condition
                hLeft  = hRight;
                huLeft = -huRight;
                bLeft  = bRight;
            }
            else
            {
                // Dry/Dry case
                // Set dummy values such that the result is zero
                hLeft   = dryTol_;
                huLeft  = T(0.0);
                bLeft   = T(0.0);
                hRight  = dryTol_;
                huRight = T(0.0);
                bRight  = T(0.0);
            }

            // Velocity on the left side of the edge
            T uLeft = huLeft / hLeft; // 1 FLOP (div)
            // Velocity on the right side of the edge
            T uRight = huRight / hRight; // 1 FLOP (div)

            /// Wave speeds of the f-waves
            T waveSpeeds0 = T(0.0);
            T waveSpeeds1 = T(0.0);

            fWaveComputeWaveSpeeds(
                hLeft, hRight, huLeft, huRight, uLeft, uRight, bLeft, bRight, waveSpeeds0, waveSpeeds1
            ); // 20 FLOPs (incl. 3 sqrt, 1 div, 2 min/max)

            // Variables to store the two f-waves
            T fWaves0 = T(0.0);
            T fWaves1 = T(0.0);

            // Compute the decomposition into f-waves
            fWaveComputeWaveDecomposition(
                hLeft, hRight, huLeft, huRight, uLeft, uRight, bLeft, bRight, waveSpeeds0, waveSpeeds1, fWaves0, fWaves1
            ); // 23 FLOPs (incl. 1 div)

            // Compute the net-updates
            o_hUpdateLeft   = T(0.0);
            o_hUpdateRight  = T(0.0);
            o_huUpdateLeft  = T(0.0);
            o_huUpdateRight = T(0.0);

            // 1st wave family
            if (waveSpeeds0 < -zeroTol_)
            { // Left going
                o_hUpdateLeft += fWaves0;
                o_huUpdateLeft += fWaves0 * waveSpeeds0; // 3 FLOPs (assume left going wave ...)
            }
            else if (waveSpeeds0 > zeroTol_)
            { // Right going
                o_hUpdateRight += fWaves0;
                o_huUpdateRight += fWaves0 * waveSpeeds0;
            }
            else
            { // Split waves, if waveSpeeds0 close to 0
                o_hUpdateLeft += T(0.5) * fWaves0;
                o_huUpdateLeft += T(0.5) * fWaves0 * waveSpeeds0;
                o_hUpdateRight += T(0.5) * fWaves0;
                o_huUpdateRight += T(0.5) * fWaves0 * waveSpeeds0;
            }

            // 2nd wave family
            if (waveSpeeds1 > zeroTol_)
            { // Right going
                o_hUpdateRight += fWaves1;
                o_huUpdateRight += fWaves1 * waveSpeeds1; // 3 FLOPs (assume right going wave ...)
            }
            else if (waveSpeeds1 < -zeroTol_)
            { // Left going
                o_hUpdateLeft += fWaves1;
                o_huUpdateLeft += fWaves1 * waveSpeeds1;
            }
            else
            { // Split waves
                o_hUpdateLeft += T(0.5) * fWaves1;
                o_huUpdateLeft += T(0.5) * fWaves1 * waveSpeeds1;
                o_hUpdateRight += T(0.5) * fWaves1;
                o_huUpdateRight += T(0.5) * fWaves1 * waveSpeeds1;
            }

            // Compute maximum wave speed (-> CFL-condition)
            o_maxWaveSpeed = std::max(std::abs(waveSpeeds0), std::abs(waveSpeeds1)); // 3 FLOPs (2 abs, 1 max)

            // ========================
            // 54 FLOPs (3 sqrt, 4 div, 2 abs, 3 min/max)
        }

#ifdef ENABLE_VECTORIZATION
#pragma omp declare simd
#endif
        void fWaveComputeWaveSpeeds(
            const T hLeft,
            const T hRight,
            const T huLeft,
            const T huRight,
            const T uLeft,
            const T uRight,
            const T bLeft,
            const T bRight,
            T&      o_waveSpeed0,
            T&      o_waveSpeed1
        ) const
        {
            // Helper variables for sqrt of h:
            T sqrtHLeft  = std::sqrt(hLeft);  // 1 FLOP (sqrt)
            T sqrtHRight = std::sqrt(hRight); // 1 FLOP (sqrt)

            // Compute eigenvalues of the jacobian matrices in states Q_{i-1} and Q_{i}
            T characteristicSpeed0 = uLeft - sqrtGravity_ * sqrtHLeft;   // 2 FLOPs
            T characteristicSpeed1 = uRight + sqrtGravity_ * sqrtHRight; // 2 FLOPs

            // Compute "Roe averages"
            T hRoe     = T(0.5) * (hRight + hLeft);              // 2 FLOPs
            T sqrtHRoe = std::sqrt(hRoe);                        // 1 FLOP (sqrt)
            T uRoe     = uLeft * sqrtHRoe + uRight * sqrtHRight; // 3 FLOPs
            uRoe /= sqrtHLeft + sqrtHRight;                      // 2 FLOPs (1 div)

            // Compute "Roe speeds" from Roe averages
            T roeSpeed0 = uRoe - sqrtGravity_ * sqrtHRoe; // 2 FLOPs
            T roeSpeed1 = uRoe + sqrtGravity_ * sqrtHRoe; // 2 FLOPs

            // Compute Eindfeldt speeds (returned as output parameters)
            o_waveSpeed0 = std::min(characteristicSpeed0, roeSpeed0); // 1 FLOP (min)
            o_waveSpeed1 = std::max(characteristicSpeed1, roeSpeed1); // 1 FLOP (max)

            // ==============
            // 20 FLOPs (incl. 3 sqrt, 1 div, 2 min/max)
        }

#ifdef ENABLE_VECTORIZATION
#pragma omp declare simd
#endif
        void fWaveComputeWaveDecomposition(
            const T hLeft,
            const T hRight,
            const T huLeft,
            const T huRight,
            const T uLeft,
            const T uRight,
            const T bLeft,
            const T bRight,
            const T waveSpeed0,
            const T waveSpeed1,
            T&      o_fWave0,
            T&      o_fWave1
        ) const
        {
            // Calculate modified (bathymetry) flux difference
            // f(Q_i) - f(Q_{i-1}) -> serve as right hand sides
            T fDif0 = huRight - huLeft; // 1 FLOP
            T fDif1 = huRight * uRight + halfGravity_ * hRight * hRight
                      - (huLeft * uLeft + halfGravity_ * hLeft * hLeft); // 9 FLOPs

            // \delta x \Psi[2]
            fDif1 += halfGravity_ * (hRight + hLeft) * (bRight - bLeft); // 5 FLOPs

            // Solve linear system of equations to obtain f-waves:
            // (       1            1      ) ( o_fWave0 ) = ( fDif0 )
            // ( waveSpeed0     waveSpeed1 ) ( o_fWave1 )   ( fDif1 )

            // Compute the inverse of the wave speed difference:
            T inverseSpeedDiff = T(1.0) / (waveSpeed1 - waveSpeed0); // 2 FLOPs (1 div)
            // Compute f-waves:
            o_fWave0 = (waveSpeed1 * fDif0 - fDif1) * inverseSpeedDiff;  // 3 FLOPs
            o_fWave1 = (-waveSpeed0 * fDif0 + fDif1) * inverseSpeedDiff; // 3 FLOPs

            // =========
            // 23 FLOPs in total (incl. 1 div)
        }
    };

} // namespace Solvers
