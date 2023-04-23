/**
 * FWaveSolver.hpp
 *
 ****
 **** F-Wave Riemann Solver for the Shallow Water Equation
 ****
 *
 *  Created on: Aug 25, 2011
 *  Last Update: Feb 18, 2012
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
 *   @article{bale2002wave,
 *            title={A wave propagation method for conservation laws and balance laws with spatially varying flux
 *functions}, author={Bale, D.S. and LeVeque, R.J. and Mitran, S. and Rossmanith, J.A.}, journal={SIAM Journal on
 *Scientific Computing}, volume={24}, number={3}, pages={955--978}, year={2002}, publisher={Citeseer}}
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

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

#include "WavePropagationSolver.hpp"

#include "Tools/HelperFunctions.hpp"

namespace Solvers
{

    /**
     * FWave Riemann Solver for the Shallow Water Equations.
     *
     * T should be double or float.
     */
    template <class T>
    class FWaveSolver: public WavePropagationSolver<T>
    {
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

        void computeWaveSpeeds(T o_waveSpeeds[2 * STRIDE]) const
        {
            for (size_t i = 0; i < STRIDE; i++)
            {
                // Compute eigenvalues of the Jacobian matrices in states Q_{i-1} and Q_{i}
                T characteristicSpeeds[2]{};
                characteristicSpeeds[0] = uLeft_[i] - std::sqrt(gravity_[i] * hLeft_[i]);
                characteristicSpeeds[1] = uRight_[i] + std::sqrt(gravity_[i] * hRight_[i]);

                // Compute "Roe speeds"
                T hRoe = T(0.5) * (hRight_[i] + hLeft_[i]);
                T uRoe = uLeft_[i] * std::sqrt(hLeft_[i]) + uRight_[i] * std::sqrt(hRight_[i]);
                uRoe /= std::sqrt(hLeft_[i]) + std::sqrt(hRight_[i]);

                T roeSpeeds[2]{};
                roeSpeeds[0] = uRoe - std::sqrt(gravity_[i] * hRoe);
                roeSpeeds[1] = uRoe + std::sqrt(gravity_[i] * hRoe);

                // Compute eindfeldt speeds
                T einfeldtSpeeds[2]{};
                einfeldtSpeeds[0] = std::min(characteristicSpeeds[0], roeSpeeds[0]);
                einfeldtSpeeds[1] = std::max(characteristicSpeeds[1], roeSpeeds[1]);

                // Set wave speeds                           // i = 0 , 1, 2, 3
                o_waveSpeeds[0 + i * 2] = einfeldtSpeeds[0]; //     0 , 2, 4, 6
                o_waveSpeeds[1 + i * 2] = einfeldtSpeeds[1]; //     1 , 3, 5, 7
            }
        }

        /**
         * Compute the decomposition into f-Waves.
         *
         * @param waveSpeeds speeds of the linearized waves (eigenvalues).
         * @param o_fWaves will be set to: Decomposition into f-Waves.
         */

        void computeWaveDecomposition(const T waveSpeeds[2 * STRIDE], T o_fWaves[2][2], const int idx) const
        {
            // Eigenvalues***********************************************************************************************
            // Computed somewhere before.
            // An option would be to use the char. Speeds:
            //
            // lambda^1 = u_{i-1} - sqrt(g * h_{i-1})
            // lambda^2 = u_i     + sqrt(g * h_i)
            // Matrix of right
            // eigenvectors******************************************************************************
            //     1                              1
            // R =
            //     u_{i-1} - sqrt(g * h_{i-1})    u_i + sqrt(g * h_i)
            // **********************************************************************************************************
            //                                                                      u_i + sqrt(g * h_i)              -1
            // R^{-1} = 1 / (u_i - sqrt(g * h_i) - u_{i-1} + sqrt(g * h_{i-1}) *
            //                                                                   -( u_{i-1} - sqrt(g * h_{i-1}) )     1
            // **********************************************************************************************************
            //                hu
            // f(q) =
            //         hu^2 + 1/2 g * h^2
            // **********************************************************************************************************
            //                                    0
            // \delta x \Psi =
            //                  -g * 1/2 * (h_i + h_{i-1}) * (b_i - b_{i+1})
            // **********************************************************************************************************
            // beta = R^{-1} * (f(Q_i) - f(Q_{i-1}) - \delta x \Psi)
            // **********************************************************************************************************

            // assert: wave speed of the 1st wave family should be less than the speed of the 2nd wave family.
            assert(waveSpeeds[idx * 2 + 0] < waveSpeeds[idx * 2 + 1]);

            T lambdaDif = waveSpeeds[idx * 2 + 1] - waveSpeeds[idx * 2 + 0];

            // assert: no division by zero
            assert(std::abs(lambdaDif) > zeroTol_[idx]);

            // Compute the inverse matrix R^{-1}
            T Rinv[2][2]{};

            T oneDivLambdaDif = T(1.0) / lambdaDif;
            Rinv[0][0]        = oneDivLambdaDif * waveSpeeds[idx * 2 + 1];
            Rinv[0][1]        = -oneDivLambdaDif;

            Rinv[1][0] = oneDivLambdaDif * -waveSpeeds[idx * 2 + 0];
            Rinv[1][1] = oneDivLambdaDif;

            // Right hand side
            T fDif[2]{};

            // Calculate modified (bathymetry!) flux difference
            // f(Q_i) - f(Q_{i-1})
            fDif[0] = huRight_[idx] - huLeft_[idx];
            fDif[1] = huRight_[idx] * uRight_[idx] + T(0.5) * gravity_[idx] * hRight_[idx] * hRight_[idx]
                      - (huLeft_[idx] * uLeft_[idx] + T(0.5) * gravity_[idx] * hLeft_[idx] * hLeft_[idx]);

            // \delta x \Psi[2]
            T psi = -gravity_[idx] * T(0.5) * (hRight_[idx] + hLeft_[idx]) * (bRight_[idx] - bLeft_[idx]);
            fDif[1] -= psi;

            // Solve linear equations
            T beta[2]{};
            beta[0] = Rinv[0][0] * fDif[0] + Rinv[0][1] * fDif[1];
            beta[1] = Rinv[1][0] * fDif[0] + Rinv[1][1] * fDif[1];

            // Return f-waves
            o_fWaves[0][0] = beta[0];
            o_fWaves[0][1] = beta[0] * waveSpeeds[idx * 2 + 0];

            o_fWaves[1][0] = beta[1];
            o_fWaves[1][1] = beta[1] * waveSpeeds[idx * 2 + 1];
        }

        /**
         * Compute net updates for the cell on the left/right side of the edge.
         * Its assumed that the member variables are set already.
         *
         * @param waveSpeeds speeds of the linearized waves (eigenvalues).
         *
         * @param o_hUpdateLeft will be set to: Net-update for the height of the cell on the left side of the edge.
         * @param o_hUpdateRight will be set to: Net-update for the height of the cell on the right side of the edge.
         * @param o_huUpdateLeft will be set to: Net-update for the momentum of the cell on the left side of the edge.
         * @param o_huUpdateRight will be set to: Net-update for the momentum of the cell on the right side of the edge.
         * @param o_maxWaveSpeed will be set to: Maximum (linearized) wave speed -> Should be used in the CFL-condition.
         */

        /**
         * Compute net updates for the cell on the left/right side of the edge.
         * Its assumed that the member variables are set already.
         *
         * @param waveSpeeds speeds of the linearized waves (eigenvalues).
         *
         * @param o_hUpdateLeft will be set to: Net-update for the height of the cell on the left side of the edge.
         * @param o_hUpdateRight will be set to: Net-update for the height of the cell on the right side of the edge.
         * @param o_huUpdateLeft will be set to: Net-update for the momentum of the cell on the left side of the edge.
         * @param o_huUpdateRight will be set to: Net-update for the momentum of the cell on the right side of the edge.
         * @param o_maxWaveSpeed will be set to: Maximum (linearized) wave speed -> Should be used in the CFL-condition.
         */

        void computeNetUpdatesWithWaveSpeeds(
            const T   waveSpeeds[2 * STRIDE],
            T&        o_hUpdateLeft,
            T&        o_hUpdateRight,
            T&        o_huUpdateLeft,
            T&        o_huUpdateRight,
            T&        o_maxWaveSpeed,
            const int idx
        )
        {
            // Reset net updates
            o_hUpdateLeft   = T(0.0);
            o_hUpdateRight  = T(0.0);
            o_huUpdateLeft  = T(0.0);
            o_huUpdateRight = T(0.0);

            //! Where to store the two f-waves
            T fWaves[2][2];

            // Compute the decomposition into f-waves
            computeWaveDecomposition(waveSpeeds, fWaves, idx);

            // Compute the net-updates
            // 1st wave family
            if (waveSpeeds[idx * 2] < -zeroTol_[idx])
            { // Left going
                o_hUpdateLeft += fWaves[0][0];
                o_huUpdateLeft += fWaves[0][1];
            }
            else if (waveSpeeds[idx * 2] > zeroTol_[idx])
            { // Right going
                o_hUpdateRight += fWaves[0][0];
                o_huUpdateRight += fWaves[0][1];
            }
            else
            { // Split waves
                o_hUpdateLeft += T(0.5) * fWaves[0][0];
                o_huUpdateLeft += T(0.5) * fWaves[0][1];
                o_hUpdateRight += T(0.5) * fWaves[0][0];
                o_huUpdateRight += T(0.5) * fWaves[0][1];
            }

            // 2nd wave family
            if (waveSpeeds[idx * 2 + 1] < -zeroTol_[idx])
            { // Left going
                o_hUpdateLeft += fWaves[1][0];
                o_huUpdateLeft += fWaves[1][1];
            }
            else if (waveSpeeds[idx * 2 + 1] > zeroTol_[idx])
            { // Right going
                o_hUpdateRight += fWaves[1][0];
                o_huUpdateRight += fWaves[1][1];
            }
            else
            { // Split waves
                o_hUpdateLeft += T(0.5) * fWaves[1][0];
                o_huUpdateLeft += T(0.5) * fWaves[1][1];
                o_hUpdateRight += T(0.5) * fWaves[1][0];
                o_huUpdateRight += T(0.5) * fWaves[1][1];
            }

            // Compute maximum wave speed (-> CFL-condition)
            o_maxWaveSpeed = std::max(std::fabs(waveSpeeds[idx * 2]), std::fabs(waveSpeeds[idx * 2 + 1]));
        }

      protected:
        void determineWetDryState() override
        {
            for (size_t i = 0; i < STRIDE; i++)
            {
                // Determine the wet/dry state
                if (hLeft_[i] < dryTol_[i] && hRight_[i] < dryTol_[i])
                { // Both cells are dry
                    wetDryState_[i] = WavePropagationSolver<T>::WetDryState::DryDry;
                }
                else if (hLeft_[i] < dryTol_[i])
                { // Left cell dry, right cell wet
                    uRight_[i] = huRight_[i] / hRight_[i];

                    // Set wall boundary conditions.
                    // This is not correct in the case of inundation problems.
                    hLeft_[i]       = hRight_[i];
                    bLeft_[i]       = bRight_[i];
                    huLeft_[i]      = -huRight_[i];
                    uLeft_[i]       = -uRight_[i];
                    wetDryState_[i] = WavePropagationSolver<T>::WetDryState::DryWetWall;
                }
                else if (hRight_[i] < dryTol_[i])
                { // Left cell wet, right cell dry
                    uLeft_[i] = huLeft_[i] / hLeft_[i];

                    // Set wall boundary conditions.
                    // This is not correct in the case of inundation problems.
                    hRight_[i]      = hLeft_[i];
                    bRight_[i]      = bLeft_[i];
                    huRight_[i]     = -huLeft_[i];
                    uLeft_[i]       = -uRight_[i];
                    wetDryState_[i] = WavePropagationSolver<T>::WetDryState::WetDryWall;
                }
                else
                { // Both cells wet
                    uLeft_[i]  = huLeft_[i] / hLeft_[i];
                    uRight_[i] = huRight_[i] / hRight_[i];

                    wetDryState_[i] = WavePropagationSolver<T>::WetDryState::WetWet;
                }
            }
        }

      public:
        /**
         * Constructor of the f-Wave solver with optional parameters.
         *
         * @param dryTolerance numerical definition of "dry".
         * @param gravity gravity constant.
         * @param zeroTolerance numerical definition of zero.
         */
        FWaveSolver(
            T dryTolerance  = static_cast<T>(0.01),
            T gravity       = static_cast<T>(9.81),
            T zeroTolerance = static_cast<T>(0.000000001)
        ):
            WavePropagationSolver<T>(dryTolerance, gravity, zeroTolerance)
        {
        }

        ~FWaveSolver() override = default;

        void computeNetUpdates(
            const T& hLeft,
            const T& hRight,
            const T& huLeft,
            const T& huRight,
            const T& bLeft,
            const T& bRight,
            T&       o_hUpdateLeft,
            T&       o_hUpdateRight,
            T&       o_huUpdateLeft,
            T&       o_huUpdateRight,
            T&       o_maxWaveSpeed
        ) override
        {

            for (size_t i = 0; i < STRIDE; i++)
            {
                // Set speeds to zero (will be determined later)
                uLeft_[i]  = 0;
                uRight_[i] = 0;

                // Reset the maximum wave speed
                PTR(o_maxWaveSpeed)[i] = 0;
            }

            // Store parameters to member variables
            // NOTE serialized inside
            WavePropagationSolver<T>::storeParameters(hLeft, hRight, huLeft, huRight, bLeft, bRight);

            // Determine the wet/dry state and compute local variables correspondingly
            // NOTE serialized inside <-------------------------------------------------
            determineWetDryState();

            for (size_t i = 0; i < STRIDE; i++)
            {
                // Zero updates and return in the case of dry cells
                if (wetDryState_[i] == WavePropagationSolver<T>::WetDryState::DryDry)
                {
                    PTR(o_hUpdateLeft)[i]   = T(0.0);
                    PTR(o_hUpdateRight)[i]  = T(0.0);
                    PTR(o_huUpdateLeft)[i]  = T(0.0);
                    PTR(o_huUpdateRight)[i] = T(0.0);
                    return;
                }
            }

            //! Wave speeds of the f-waves
            T waveSpeeds[2 * STRIDE];
            // Compute the wave speeds
            // NOTE serialized inside <-------------------------------------------------
            computeWaveSpeeds(waveSpeeds);

            // Use the wave speeds to compute the net-updates
            for (size_t i = 0; i < STRIDE; i++)
            {
                computeNetUpdatesWithWaveSpeeds(
                    PTR(waveSpeeds)[0],
                    PTR(o_hUpdateLeft)[i],
                    PTR(o_hUpdateRight)[i],
                    PTR(o_huUpdateLeft)[i],
                    PTR(o_huUpdateRight)[i],
                    PTR(o_maxWaveSpeed)[i],
                    i
                );

                // Zero ghost updates (wall boundary)
                if (wetDryState_[i] == WavePropagationSolver<T>::WetDryState::WetDryWall)
                {
                    PTR(o_hUpdateRight)[i]  = 0;
                    PTR(o_huUpdateRight)[i] = 0;
                }
                else if (wetDryState_[i] == WavePropagationSolver<T>::WetDryState::DryWetWall)
                {
                    PTR(o_hUpdateLeft)[i]  = 0;
                    PTR(o_huUpdateLeft)[i] = 0;
                }
            }
        }
    };

} // namespace Solvers
