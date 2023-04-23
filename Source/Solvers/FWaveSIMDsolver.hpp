/**
 * FWaveSIMDsolver.hpp
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

#if defined(ENABLE_VECTORIZATION) && defined(ENABLE_VECTORIZATION_WITH_SIMD)

// #pragma ICPC target("avx2")
// #define __AVX__
// #define __AVX2__
#include <algorithm>
#include <cassert>
#include <cmath>
#include <immintrin.h>
#include <iostream>

#include "SIMD_defs.hpp"
#include "WavePropagationSolver.hpp"

#include "Tools/HelperFunctions.hpp"

// #include "SIMDDefinitions.hpp" it says never include this directly  ¯\_(ツ)_/¯
// #include "SIMDTypes.hpp"

namespace Solvers
{

    /**
     * FWave Riemann Solver for the Shallow Water Equations.
     *
     * T should be double or float.
     */
    template <class T>
    class FWaveSIMDsolver: public WavePropagationSolver<T>
    {
      private:
        using WavePropagationSolver<T>::dryTol_v;
        using WavePropagationSolver<T>::gravity_v;
        using WavePropagationSolver<T>::zeroTol_v;

        using WavePropagationSolver<T>::hLeft_vv;
        using WavePropagationSolver<T>::hRight_vv;
        using WavePropagationSolver<T>::huLeft_vv;
        using WavePropagationSolver<T>::huRight_vv;
        using WavePropagationSolver<T>::bLeft_vv;
        using WavePropagationSolver<T>::bRight_vv;
        using WavePropagationSolver<T>::uLeft_vv;
        using WavePropagationSolver<T>::uRight_vv;

        using WavePropagationSolver<T>::wetDryState_v;

        using WavePropagationSolver<T>::WetWet_v;
        using WavePropagationSolver<T>::DryWetWall_v;
        using WavePropagationSolver<T>::WetDryWall_v;
        using WavePropagationSolver<T>::DryDry_v;

        void computeWaveSpeeds(__m256d& o_waveSpeeds_xv, __m256d& o_waveSpeeds_yv) const
        {
            __m256d characteristicSpeeds_x, characteristicSpeeds_y;
            __m256d mulop          = _mm256_mul_pd(gravity_v, hLeft_vv);
            characteristicSpeeds_x = uLeft_vv - _mm256_sqrt_pd(mulop);
            mulop                  = _mm256_mul_pd(gravity_v, hRight_vv);
            characteristicSpeeds_y = uRight_vv + _mm256_sqrt_pd(mulop);

            __m256d addop                  = _mm256_add_pd(hRight_vv, hLeft_vv);
            __m256d hRoe_v                 = _mm256_mul_pd(_mm256_set1_pd(0.5), addop);
            __m256d sqrtL                  = _mm256_sqrt_pd(hLeft_vv);
            __m256d sqrtR                  = _mm256_sqrt_pd(hRight_vv);
            __m256d uRoe_v                 = _mm256_mul_pd(uLeft_vv, sqrtL) + _mm256_mul_pd(uRight_vv, sqrtR);
            __m256d is_sqrtL_or_sqrtR_zero = _mm256_cmp_pd((sqrtL + sqrtR), _mm256_set1_pd(0.0), _CMP_EQ_OS);
            uRoe_v                         = _mm256_div_pd(
                uRoe_v, _mm256_add_pd((sqrtL + sqrtR), _mm256_and_pd(is_sqrtL_or_sqrtR_zero, _mm256_set1_pd(1.0)))
            );
            uRoe_v = _mm256_andnot_pd(is_sqrtL_or_sqrtR_zero, uRoe_v); // 0   1.5    1.2    0
            __m256d roeSpeeds_x, roeSpeeds_y;
            mulop       = _mm256_mul_pd(gravity_v, hRoe_v); // 0   4.5       4.5    0
            roeSpeeds_x = uRoe_v - _mm256_sqrt_pd(mulop);   // 0   1.5       1.2    0
            roeSpeeds_y = uRoe_v + _mm256_sqrt_pd(mulop);

            // REVIEW[epic=SWE,seq=55] optimize using permute
            // TODO[epic=SWE,seq=58] move to aligned memory for perf improve
            // Compute eindfeldt speeds
            // Set wave speeds
            o_waveSpeeds_xv = _mm256_min_pd(characteristicSpeeds_x, roeSpeeds_x);
            o_waveSpeeds_yv = _mm256_max_pd(characteristicSpeeds_y, roeSpeeds_y);
        }

        /**
         * Compute the decomposition into f-Waves.
         *
         * @param waveSpeeds speeds of the linearized waves (eigenvalues).
         * @param o_fWaves will be set to: Decomposition into f-Waves.
         */
        void computeWaveDecomposition(
            const __m256d& o_waveSpeeds_xv,
            const __m256d& o_waveSpeeds_yv,
            T              o_fWaves00[4],
            T              o_fWaves01[4],
            T              o_fWaves10[4],
            T              o_fWaves11[4]
        ) const
        {

            // assert: no division by zero
            // FIXME[epic=SWE,seq=55] enable assert in debug mode
            //   #ifndef NDEBUG
            //    auto asrt_v = _mm256_cmplt_epi64_mask(_mm256_abs_epi64(_mm256_cvtpd_epu64(lambdaDif)),
            //    _mm256_set1_epi64x(0));
            //   #endif

            // Compute the inverse matrix R^{-1}
            __m256d Rinv00, Rinv01, Rinv10, Rinv11;
            __m256d lambdaDif       = o_waveSpeeds_yv - o_waveSpeeds_xv;
            __m256d is_lambda_zero  = _mm256_cmp_pd(lambdaDif, _mm256_set1_pd(0.0), _CMP_EQ_OS);
            __m256d oneDivLambdaDif = _mm256_div_pd(
                _mm256_set1_pd(1.0), _mm256_add_pd(lambdaDif, _mm256_and_pd(is_lambda_zero, _mm256_set1_pd(1.0)))
            );
            oneDivLambdaDif = _mm256_andnot_pd(is_lambda_zero, oneDivLambdaDif);
            Rinv00          = _mm256_mul_pd(oneDivLambdaDif, o_waveSpeeds_yv);
            Rinv01          = -oneDivLambdaDif;

            Rinv10 = _mm256_mul_pd(oneDivLambdaDif, -o_waveSpeeds_xv);
            Rinv11 = oneDivLambdaDif;

            // Right hand side
            __m256d fDif_x = _mm256_set1_pd(0.0);
            __m256d fDif_y = _mm256_set1_pd(0.0);

            // Calculate modified (bathymetry!) flux difference
            // f(Q_i) - f(Q_{i-1})
            fDif_x        = huRight_vv - huLeft_vv;
            __m256d mulop = _mm256_mul_pd(_mm256_set1_pd(0.5), gravity_v);
            fDif_y        = (huRight_vv * uRight_vv) + (mulop * hRight_vv * hRight_vv) - (huLeft_vv * uLeft_vv)
                     - (mulop * hLeft_vv * hLeft_vv);

            __m256d psi = (-mulop) * (hRight_vv + hLeft_vv) * (bRight_vv - bLeft_vv);
            fDif_y      = fDif_y - psi;

            // Solve linear equations
            // T beta[2]{};
            __m256d beta_x = _mm256_set1_pd(0.0);
            __m256d beta_y = _mm256_set1_pd(0.0);
            beta_x         = Rinv00 * fDif_x + Rinv01 * fDif_y;
            beta_y         = Rinv10 * fDif_x + Rinv11 * fDif_y;

            // Return f-waves
            //  o_fWaves00= beta_x;
            _mm256_storeu_pd(o_fWaves00, beta_x);
            mulop = _mm256_mul_pd(beta_x, o_waveSpeeds_xv);
            _mm256_storeu_pd(o_fWaves01, mulop);

            _mm256_storeu_pd(o_fWaves10, beta_y);
            mulop = _mm256_mul_pd(beta_y, o_waveSpeeds_yv);
            _mm256_storeu_pd(o_fWaves11, mulop);
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
            const __m256d& o_waveSpeeds_xv,
            const __m256d& o_waveSpeeds_yv,
            __m256d&       o_hUpdateLeft_v,
            __m256d&       o_hUpdateRight_v,
            __m256d&       o_huUpdateLeft_v,
            __m256d&       o_huUpdateRight_v,
            __m256d&       o_maxWaveSpeed_v
        )
        {
            o_hUpdateLeft_v   = _mm256_setzero_pd();
            o_hUpdateRight_v  = _mm256_setzero_pd();
            o_huUpdateLeft_v  = _mm256_setzero_pd();
            o_huUpdateRight_v = _mm256_setzero_pd();

            //! Where to store the two f-waves
            T fWaves00[4]{};
            T fWaves01[4]{};
            T fWaves11[4]{};
            T fWaves10[4]{};

            // Compute the decomposition into f-waves
            computeWaveDecomposition(o_waveSpeeds_xv, o_waveSpeeds_yv, fWaves00, fWaves01, fWaves10, fWaves11);

            // Compute the net-updates
            // 1st wave family
            // if (waveSpeeds[0] < -zeroTol_) { // Left going

            __m256d cmp_v, cmp2_v, z;
            z = _mm256_set1_pd(-0.0);

            __m256d fWaves00v = _mm256_loadu_pd(fWaves00);
            __m256d fWaves01v = _mm256_loadu_pd(fWaves01);
            __m256d fWaves11v = _mm256_loadu_pd(fWaves11);
            __m256d fWaves10v = _mm256_loadu_pd(fWaves10);
            cmp_v             = _mm256_cmp_pd(o_waveSpeeds_xv, -zeroTol_v, _CMP_LT_OS);

            __m256d andop = _mm256_and_pd(cmp_v, fWaves00v);

            o_hUpdateLeft_v  = _mm256_add_pd(o_hUpdateLeft_v, andop); // fWaves[0v][0];
            andop            = _mm256_and_pd(cmp_v, fWaves01v);
            o_huUpdateLeft_v = _mm256_add_pd(o_huUpdateLeft_v, andop);

            //} else if (waveSpeeds[0] > zeroTol_) { // Right going
            cmp2_v = _mm256_cmp_pd(o_waveSpeeds_xv, zeroTol_v, _CMP_GT_OS);

            andop             = _mm256_and_pd(cmp2_v, fWaves00v);
            o_hUpdateRight_v  = _mm256_add_pd(o_hUpdateRight_v, andop);
            andop             = _mm256_and_pd(cmp2_v, fWaves01v);
            o_huUpdateRight_v = _mm256_add_pd(o_huUpdateRight_v, andop);

            //} else { // Split waves
            cmp_v  = _mm256_cmp_pd(o_waveSpeeds_xv, -zeroTol_v, _CMP_GE_OS);
            cmp2_v = _mm256_cmp_pd(o_waveSpeeds_xv, zeroTol_v, _CMP_LE_OS);
            cmp2_v = _mm256_and_pd(cmp_v, cmp2_v);

            __m256d mulop   = _mm256_mul_pd(_mm256_set1_pd(0.5), fWaves00v);
            andop           = _mm256_and_pd(cmp2_v, mulop);
            o_hUpdateLeft_v = _mm256_add_pd(o_hUpdateLeft_v, andop); //[0][0]

            mulop            = _mm256_mul_pd(_mm256_set1_pd(0.5), fWaves01v);
            andop            = _mm256_and_pd(cmp2_v, mulop);
            o_huUpdateLeft_v = _mm256_add_pd(o_huUpdateLeft_v, andop); //[0][1]

            mulop            = _mm256_mul_pd(_mm256_set1_pd(0.5), fWaves00v);
            andop            = _mm256_and_pd(cmp2_v, mulop);
            o_hUpdateRight_v = _mm256_add_pd(o_hUpdateRight_v, andop); //[0][0]

            mulop             = _mm256_mul_pd(_mm256_set1_pd(0.5), fWaves01v);
            andop             = _mm256_and_pd(cmp2_v, mulop);
            o_huUpdateRight_v = _mm256_add_pd(o_huUpdateRight_v, andop); //[0][1]

            //}//////////////////////////////////////////////

            cmp_v             = _mm256_cmp_pd(o_waveSpeeds_yv, -zeroTol_v, _CMP_LT_OS);
            andop             = _mm256_and_pd(cmp_v, fWaves10v);
            o_hUpdateLeft_v   = _mm256_add_pd(o_hUpdateLeft_v, andop); // fWaves[0v][0];
            andop             = _mm256_and_pd(cmp_v, fWaves11v);
            o_huUpdateLeft_v  = _mm256_add_pd(o_huUpdateLeft_v, andop);
            cmp2_v            = _mm256_cmp_pd(o_waveSpeeds_yv, zeroTol_v, _CMP_GT_OS);
            andop             = _mm256_and_pd(cmp2_v, fWaves10v);
            o_hUpdateRight_v  = _mm256_add_pd(o_hUpdateRight_v, andop);
            andop             = _mm256_and_pd(cmp2_v, fWaves11v);
            o_huUpdateRight_v = _mm256_add_pd(o_huUpdateRight_v, andop);
            cmp_v             = _mm256_cmp_pd(o_waveSpeeds_yv, -zeroTol_v, _CMP_GE_OS);
            cmp2_v            = _mm256_cmp_pd(o_waveSpeeds_yv, zeroTol_v, _CMP_LE_OS);
            cmp2_v            = _mm256_and_pd(cmp_v, cmp2_v);

            mulop             = _mm256_mul_pd(_mm256_set1_pd(0.5), fWaves10v);
            andop             = _mm256_and_pd(cmp2_v, mulop);
            o_hUpdateLeft_v   = _mm256_add_pd(o_hUpdateLeft_v, andop); //[0][0]
            mulop             = _mm256_mul_pd(_mm256_set1_pd(0.5), fWaves11v);
            andop             = _mm256_and_pd(cmp2_v, mulop);
            o_huUpdateLeft_v  = _mm256_add_pd(o_huUpdateLeft_v, andop); //[0][1]
            mulop             = _mm256_mul_pd(_mm256_set1_pd(0.5), fWaves10v);
            andop             = _mm256_and_pd(cmp2_v, mulop);
            o_hUpdateRight_v  = _mm256_add_pd(o_hUpdateRight_v, andop); //[0][0]
            mulop             = _mm256_mul_pd(_mm256_set1_pd(0.5), fWaves11v);
            andop             = _mm256_and_pd(cmp2_v, mulop);
            o_huUpdateRight_v = _mm256_add_pd(o_huUpdateRight_v, andop); //[0][1]

            // Compute maximum wave speed (-> CFL-condition)
            __m256d andnot_x = _mm256_andnot_pd(z, o_waveSpeeds_xv);
            __m256d andnot_y = _mm256_andnot_pd(z, o_waveSpeeds_yv);
            o_maxWaveSpeed_v = _mm256_max_pd(andnot_x, andnot_y);

            __m256d y  = _mm256_permute2f128_pd(o_maxWaveSpeed_v, o_maxWaveSpeed_v, 1); // permute 128-bit values
            __m256d m1 = _mm256_max_pd(o_maxWaveSpeed_v, y); // m1[0] = max(x[0], x[2]), m1[1] = max(x[1], x[3]), etc.
            __m256d m2 = _mm256_permute_pd(m1, 5);           // set m2[0] = m1[1], m2[1] = m1[0], etc.
            o_maxWaveSpeed_v = _mm256_max_pd(m1, m2);
        }

      protected:
        void determineWetDryState() override
        {
            __m256d cmp_v       = _mm256_cmp_pd(hLeft_vv, dryTol_v, _CMP_LT_OS);
            __m256d cmp2_v      = _mm256_cmp_pd(hRight_vv, dryTol_v, _CMP_LT_OS);
            __m256d fnd0_v      = _mm256_cmp_pd(hRight_vv, _mm256_set1_pd(0.0), _CMP_EQ_OS);
            __m256d fnd02_v     = _mm256_cmp_pd(hLeft_vv, _mm256_set1_pd(0.0), _CMP_EQ_OS);
            __m256d cmp3_v      = _mm256_and_pd(cmp2_v, cmp_v);
            __m256d andop       = _mm256_and_pd(fnd0_v, _mm256_set1_pd(1.0));
            __m256d uRightDiv_v = _mm256_div_pd(huRight_vv, (hRight_vv + andop));
            andop               = _mm256_and_pd(fnd02_v, _mm256_set1_pd(1.0));
            __m256d uLeftDiv_v  = _mm256_div_pd(huLeft_vv, (hLeft_vv + andop));

            __m256d xor_v  = _mm256_xor_pd(cmp2_v, cmp_v);
            __m256d cmp4_v = _mm256_and_pd(xor_v, cmp_v);

            uRight_vv = _mm256_andnot_pd(cmp4_v, uRight_vv) + _mm256_and_pd(cmp4_v, uRightDiv_v);

            // Set wall boundary conditions.
            // This is not correct in the case of inundation problems.
            hLeft_vv  = _mm256_andnot_pd(cmp4_v, hLeft_vv) + _mm256_and_pd(cmp4_v, hRight_vv);
            bLeft_vv  = _mm256_andnot_pd(cmp4_v, bLeft_vv) + _mm256_and_pd(cmp4_v, bRight_vv);
            huLeft_vv = _mm256_andnot_pd(cmp4_v, huLeft_vv) + _mm256_and_pd(cmp4_v, -huRight_vv);
            uLeft_vv  = _mm256_andnot_pd(cmp4_v, uLeft_vv) + _mm256_and_pd(cmp4_v, -uRight_vv);

            __m256d cmp5_v = _mm256_and_pd(xor_v, cmp_v);

            uLeft_vv = _mm256_andnot_pd(cmp5_v, uLeft_vv) + _mm256_and_pd(cmp5_v, uLeftDiv_v);

            // Set wall boundary conditions.
            // This is not correct in the case of inundation problems.
            hRight_vv      = _mm256_andnot_pd(cmp5_v, hRight_vv) + _mm256_and_pd(cmp5_v, hLeft_vv);
            bRight_vv      = _mm256_andnot_pd(cmp5_v, bRight_vv) + _mm256_and_pd(cmp5_v, bLeft_vv);
            huRight_vv     = _mm256_andnot_pd(cmp5_v, huRight_vv) + _mm256_and_pd(cmp5_v, -huLeft_vv);
            uLeft_vv       = _mm256_andnot_pd(cmp5_v, uLeft_vv) + _mm256_and_pd(cmp5_v, -uRight_vv);
            __m256d cmp6_v = _mm256_cmp_pd(hLeft_vv, dryTol_v, _CMP_GE_OS);
            __m256d cmp7_v = _mm256_cmp_pd(hRight_vv, dryTol_v, _CMP_GE_OS);
            __m256d cmp8_v = _mm256_or_pd(cmp7_v, cmp6_v);

            uLeft_vv  = _mm256_andnot_pd(cmp8_v, uLeft_vv) + _mm256_and_pd(cmp8_v, uLeftDiv_v);
            uRight_vv = _mm256_andnot_pd(cmp8_v, uRight_vv) + _mm256_and_pd(cmp8_v, uRightDiv_v);

            wetDryState_v = _mm256_and_pd(cmp3_v, DryDry_v) + _mm256_and_pd(cmp4_v, DryWetWall_v)
                            + _mm256_and_pd(cmp5_v, WetDryWall_v) + _mm256_and_pd(cmp8_v, WetWet_v);
        }

      public:
        /**
         * Constructor of the f-Wave solver with optional parameters.
         *
         * @param dryTolerance numerical definition of "dry".
         * @param gravity gravity constant.
         * @param zeroTolerance numerical definition of zero.
         */
        FWaveSIMDsolver(
            T dryTolerance  = static_cast<T>(0.01),
            T gravity       = static_cast<T>(9.81),
            T zeroTolerance = static_cast<T>(0.000000001)
        ):
            WavePropagationSolver<T>(dryTolerance, gravity, zeroTolerance)
        {
        }

        ~FWaveSIMDsolver() override = default;

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
            // Set speeds to zero (will be determined later)
            uLeft_vv  = _mm256_setzero_pd();
            uRight_vv = _mm256_setzero_pd();

            // Store parameters to member variables
            WavePropagationSolver<T>::storeParameters(hLeft, hRight, huLeft, huRight, bLeft, bRight);

            // Determine the wet/dry state and compute local variables correspondingly
            determineWetDryState();

            // Zero ghost updates (wall boundary)
            // FIXME[epic=SWE,seq=57] optimize the following code
            // FIXME[epic=SWE,seq=51] stop calculations and return if wetDryState_v != DryDry_v
            // FIXME[epic=SWE,seq=54] vectorize the whole pipeline! from outside
            __m256d o_hUpdateLeft_v;
            __m256d o_hUpdateRight_v;
            __m256d o_huUpdateLeft_v;
            __m256d o_huUpdateRight_v;
            __m256d o_maxWaveSpeed_v = _mm256_setzero_pd();

            // Dry Dry Wall
            __m256d cmp_DryDryWall_OR_WetDryWall_OR_DryWetWall_v = _mm256_cmp_pd(wetDryState_v, DryDry_v, _CMP_NEQ_OS);
            o_hUpdateLeft_v   = _mm256_and_pd(cmp_DryDryWall_OR_WetDryWall_OR_DryWetWall_v, o_hUpdateLeft_v);
            o_hUpdateRight_v  = _mm256_and_pd(cmp_DryDryWall_OR_WetDryWall_OR_DryWetWall_v, o_hUpdateRight_v);
            o_huUpdateLeft_v  = _mm256_and_pd(cmp_DryDryWall_OR_WetDryWall_OR_DryWetWall_v, o_huUpdateLeft_v);
            o_huUpdateRight_v = _mm256_and_pd(cmp_DryDryWall_OR_WetDryWall_OR_DryWetWall_v, o_huUpdateRight_v);

            //! Wave speeds of the f-waves
            __m256d waveSpeeds_xv;
            __m256d waveSpeeds_yv;

            // Compute the wave speeds
            computeWaveSpeeds(waveSpeeds_xv, waveSpeeds_yv);

            // Use the wave speeds to compute the net-updates
            computeNetUpdatesWithWaveSpeeds(
                waveSpeeds_xv,
                waveSpeeds_yv,
                o_hUpdateLeft_v,
                o_hUpdateRight_v,
                o_huUpdateLeft_v,
                o_huUpdateRight_v,
                o_maxWaveSpeed_v
            );

            // Wet Dry Wall
            cmp_DryDryWall_OR_WetDryWall_OR_DryWetWall_v = _mm256_cmp_pd(wetDryState_v, WetDryWall_v, _CMP_NEQ_OS);
            o_hUpdateRight_v  = _mm256_and_pd(cmp_DryDryWall_OR_WetDryWall_OR_DryWetWall_v, o_hUpdateRight_v);
            o_huUpdateRight_v = _mm256_and_pd(cmp_DryDryWall_OR_WetDryWall_OR_DryWetWall_v, o_huUpdateRight_v);

            // Dry Wet Wall
            cmp_DryDryWall_OR_WetDryWall_OR_DryWetWall_v = _mm256_cmp_pd(wetDryState_v, DryWetWall_v, _CMP_NEQ_OS);
            o_hUpdateLeft_v  = _mm256_and_pd(cmp_DryDryWall_OR_WetDryWall_OR_DryWetWall_v, o_hUpdateLeft_v);
            o_huUpdateLeft_v = _mm256_and_pd(cmp_DryDryWall_OR_WetDryWall_OR_DryWetWall_v, o_huUpdateLeft_v);

            _mm256_storeu_pd(&o_maxWaveSpeed, o_maxWaveSpeed_v);
            _mm256_storeu_pd(&o_hUpdateLeft, o_hUpdateLeft_v);
            _mm256_storeu_pd(&o_huUpdateLeft, o_huUpdateLeft_v);
            _mm256_storeu_pd(&o_hUpdateRight, o_hUpdateRight_v);
            _mm256_storeu_pd(&o_huUpdateRight, o_huUpdateRight_v);
        }
    };

} // namespace Solvers

#endif