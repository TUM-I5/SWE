#pragma once

#if defined(ENABLE_VECTORIZATION) && defined(ENABLE_VECTORIZATION_WITH_SIMD)
#include <immintrin.h>
#endif

#include "SIMD_defs.hpp"

namespace Solvers
{

    /**
     * Abstract wave propagation solver for the Shallow Water Equations.
     *
     * T should be double or float.
     */
    template <class T>
    class WavePropagationSolver
    {
        // protected:
      public:
        /**
         * The wet/dry state of the Riemann-problem.
         */
        enum WetDryState : long long
        {
            DryDry,           /**< Both cells are dry. */
            WetWet,           /**< Both cells are wet. */
            WetDryInundation, /**< 1st cell: wet, 2nd cell: dry. 1st cell lies higher than the 2nd one. */
            WetDryWall, /**< 1st cell: wet, 2nd cell: dry. 1st cell lies lower than the 2nd one. Momentum is not large
                           enough to overcome the difference. */
            WetDryWallInundation, /**< 1st cell: wet, 2nd cell: dry. 1st cell lies lower than the 2nd one. Momentum is
                                     large enough to overcome the difference. */
            DryWetInundation,     /**< 1st cell: dry, 2nd cell: wet. 1st cell lies lower than the 2nd one. */
            DryWetWall, /**< 1st cell: dry, 2nd cell: wet. 1st cell lies higher than the 2nd one. Momentum is not large
                           enough to overcome the difference. */
            DryWetWallInundation /**< 1st cell: dry, 2nd cell: wet. 1st cell lies higher than the 2nd one. Momentum is
                                    large enough to overcome the difference. */
        };
        WetDryState wetDryState_[STRIDE]; //! wet/dry state of our Riemann-problem (determined by determineWetDryState).

        //! Edge-local variables.
        T dryTol_[STRIDE];  //! Numerical definition of "dry".
        T gravity_[STRIDE]; //! Gravity constant.
        T zeroTol_[STRIDE]; //! Numerical definition of zero.
        T hLeft_[STRIDE];   //! Height on the left side of the edge (could change during execution).
        T hRight_[STRIDE];  //! Height on the right side of the edge (could change during execution).
        T huLeft_[STRIDE];  //! Momentum on the left side of the edge (could change during execution).
        T huRight_[STRIDE]; //! Momentum on the right side of the edge (could change during execution).
        T bLeft_[STRIDE];   //! Bathymetry on the left side of the edge (could change during execution).
        T bRight_[STRIDE];  //! Bathymetry on the right side of the edge (could change during execution).
        T uLeft_[STRIDE];   //! Velocity on the left side of the edge (computed by determineWetDryState).
        T uRight_[STRIDE];  //! Velocity on the right side of the edge (computed by determineWetDryState).

#if defined(ENABLE_VECTORIZATION) && defined(ENABLE_VECTORIZATION_WITH_SIMD)
        // VECTORIZED SIMD
        __m256d dryTol_v;
        __m256d gravity_v;
        __m256d zeroTol_v;
        __m256d wetDryState_v;
        __m256d hLeft_vv;
        __m256d hRight_vv;
        __m256d huLeft_vv;
        __m256d huRight_vv;
        __m256d bLeft_vv;
        __m256d bRight_vv;
        __m256d uLeft_vv;
        __m256d uRight_vv;

        __m256d WetWet_v     = _mm256_set1_pd(static_cast<double>(WetDryState::WetWet));
        __m256d DryWetWall_v = _mm256_set1_pd(static_cast<double>(WetDryState::DryWetWall));
        __m256d WetDryWall_v = _mm256_set1_pd(static_cast<double>(WetDryState::WetDryWall));
        __m256d DryDry_v     = _mm256_set1_pd(static_cast<double>(WetDryState::DryDry));
#endif

        //! Determine the wet/dry-state and set local values if we have to.
        virtual void determineWetDryState() = 0;

        /**
         * Constructor of a wave propagation solver.
         *
         * @param gravity gravity constant.
         * @param dryTolerance numerical definition of "dry".
         * @param zeroTolerance numerical definition of zero.
         */
        WavePropagationSolver(T dryTolerance, T gravity, T zeroTolerance)
        {
            // FIXME[epic=SWE,seq=54] vectorize
            for (size_t i = 0; i < STRIDE; i++)
            {
                dryTol_[i]  = dryTolerance;
                gravity_[i] = gravity;
                zeroTol_[i] = zeroTolerance;
            }

#if defined(ENABLE_VECTORIZATION) && defined(ENABLE_VECTORIZATION_WITH_SIMD)
            dryTol_v  = _mm256_set1_pd(dryTolerance);
            gravity_v = _mm256_set1_pd(gravity);
            zeroTol_v = _mm256_set1_pd(zeroTolerance);
#endif
        }

        /**
         * Store parameters to member variables.
         *
         * @param hLeft height on the left side of the edge.
         * @param hRight height on the right side of the edge.
         * @param huLeft momentum on the left side of the edge.
         * @param huRight momentum on the right side of the edge.
         * @param bLeft bathymetry on the left side of the edge.
         * @param bRight bathymetry on the right side of the edge.
         */
        void storeParameters(
            const T& hLeft, const T& hRight, const T& huLeft, const T& huRight, const T& bLeft, const T& bRight
        )
        {
#if defined(ENABLE_VECTORIZATION) && defined(ENABLE_VECTORIZATION_WITH_SIMD)
            hLeft_vv   = _mm256_loadu_pd(&hLeft);
            hRight_vv  = _mm256_loadu_pd(&hRight);
            huLeft_vv  = _mm256_loadu_pd(&huLeft);
            huRight_vv = _mm256_loadu_pd(&huRight);
            bLeft_vv   = _mm256_loadu_pd(&bLeft);
            bRight_vv  = _mm256_loadu_pd(&bRight);

            // FIXME[epic=SWE,seq=53] delete when vectorization done
            _mm256_storeu_pd(hLeft_, hLeft_vv);
            _mm256_storeu_pd(hRight_, hRight_vv);
            _mm256_storeu_pd(huLeft_, huLeft_vv);
            _mm256_storeu_pd(huRight_, huRight_vv);
            _mm256_storeu_pd(bLeft_, bLeft_vv);
            _mm256_storeu_pd(bRight_, bRight_vv);
#else
            for (size_t i = 0; i < STRIDE; i++)
            {
                hLeft_[i]   = PTR(hLeft)[i];
                hRight_[i]  = PTR(hRight)[i];
                huLeft_[i]  = PTR(huLeft)[i];
                huRight_[i] = PTR(huRight)[i];
                bLeft_[i]   = PTR(bLeft)[i];
                bRight_[i]  = PTR(bRight)[i];
            }
#endif
        }

        /**
         * Store parameters to member variables.
         *
         * @param hLeft height on the left side of the edge.
         * @param hRight height on the right side of the edge.
         * @param huLeft momentum on the left side of the edge.
         * @param huRight momentum on the right side of the edge.
         * @param bLeft bathymetry on the left side of the edge.
         * @param bRight bathymetry on the right side of the edge.
         * @param uLeft velocity on the left side of the edge.
         * @param uRight velocity on the right side of the edge.
         */
        void storeParameters(
            const T& hLeft,
            const T& hRight,
            const T& huLeft,
            const T& huRight,
            const T& bLeft,
            const T& bRight,
            const T& uLeft,
            const T& uRight
        )
        {
            storeParameters(hLeft, hRight, huLeft, huRight, bLeft, bRight);

            uLeft_  = uLeft;
            uRight_ = uRight;
        }

      public:
        virtual ~WavePropagationSolver() = default;

        /**
         * Compute net updates for the cell on the left/right side of the edge.
         * This is the default method every standalone wave propagation solver should provide.
         *
         * @param hLeft height on the left side of the edge.
         * @param hRight height on the right side of the edge.
         * @param huLeft momentum on the left side of the edge.
         * @param huRight momentum on the right side of the edge.
         * @param bLeft bathymetry on the left side of the edge.
         * @param bRight bathymetry on the right side of the edge.
         *
         * @param o_hUpdateLeft will be set to: Net-update for the height of the cell on the left side of the edge.
         * @param o_hUpdateRight will be set to: Net-update for the height of the cell on the right side of the edge.
         * @param o_huUpdateLeft will be set to: Net-update for the momentum of the cell on the left side of the edge.
         * @param o_huUpdateRight will be set to: Net-update for the momentum of the cell on the right side of the edge.
         * @param o_maxWaveSpeed will be set to: Maximum (linearized) wave speed -> Should be used in the CFL-condition.
         */

        virtual void computeNetUpdates(
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
#ifdef ENABLE_AUGMENTED_RIEMANN_EIGEN_COEFFICIENTS
            ,
            T o_eigenCoefficients[3]
#endif
        ) = 0;
    };

} // namespace Solvers
