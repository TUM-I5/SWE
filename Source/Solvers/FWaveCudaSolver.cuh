/**
 * FWaveCuda.h
 *
 ****
 **** This is a pure C implementation of the standard C++ f-Wave solver (FWave.hpp).
 ****
 *
 * Created on: Jan 12, 2012
 * Last Update: Jan 12, 2012
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

#include <cmath>
#include <cstdio>

#include "Tools/RealType.hpp"

typedef RealType T;

//! numerical definition of zero.
const T zeroTol = (T)0.0000001;

//! numerical definition of a dry cell
const T dryTol = (T)1.;

namespace cuda{

inline __device__ void fWaveComputeNetUpdates(
    const T i_gravity,
    T       i_hLeft,
    T       i_hRight,
    T       i_huLeft,
    T       i_huRight,
    T       i_bLeft,
    T       i_bRight,

    T o_netUpdates[5]
);

inline __device__ void fWaveComputeWaveSpeeds(
    const T i_gravity,
    const T i_hLeft,
    const T i_hRight,
    const T i_huLeft,
    const T i_huRight,
    const T i_uLeft,
    const T i_uRight,
    const T i_bLeft,
    const T i_bRight,

    T o_waveSpeeds[2]
);

inline __device__ void fWaveComputeWaveDecomposition(
    const T i_gravity,
    const T i_hLeft,
    const T i_hRight,
    const T i_huLeft,
    const T i_huRight,
    const T i_uLeft,
    const T i_uRight,
    const T i_bLeft,
    const T i_bRight,
    const T i_waveSpeeds[2],

    T o_fWaves[2][2]
);

/**
 * Compute net updates for the cell on the left/right side of the edge
 *
 * The order of o_netUpdates is given by:
 *   0: hUpdateLeft   - Net-update for the height of the cell on the left side of the edge.
 *   1: hUpdateRight  - Net-update for the height of the cell on the right side of the edge.
 *   2: huUpdateLeft  - Net-update for the momentum of the cell on the left side of the edge
 *   3: huUpdateRight - Net-update for the momentum of the cell on the right side of the edge.
 *   4: maxWaveSpeed  - Maximum (linearized) wave speed -> Should be used in the CFL-condition.
 *
 * TODO: A wet/wet state is assumend / no boundaries are implemented within the solver to
 *       allow an execution with with a few jumps only.
 *
 * @param i_gravity gravity constant.
 * @param i_hLeft height on the left side of the edge.
 * @param i_hRight height on the right side of the edge.
 * @param i_huLeft momentum on the left side of the edge.
 * @param i_huRight momentum on the right side of the edge.
 * @param i_bLeft bathymetry on the left side of the edge.
 * @param i_bRight bathymetry on the right side of the edge.
 *
 * @param o_netUpdates will be set to the net updates.
 */
inline __device__ void fWaveComputeNetUpdates(
    const T i_gravity,
    T       i_hLeft,
    T       i_hRight,
    T       i_huLeft,
    T       i_huRight,
    T       i_bLeft,
    T       i_bRight,

    T o_netUpdates[5]
)
{
    // reset net updates
    o_netUpdates[0] = (T)0.; // hUpdateLeft
    o_netUpdates[1] = (T)0.; // hUpdateRight
    o_netUpdates[2] = (T)0.; // huUpdateLeft
    o_netUpdates[3] = (T)0.; // huUpdateRight

    // reset the maximum wave speed
    o_netUpdates[4] = (T)0.; // maxWaveSpeed

    // determine the wet dry state and corr. values, if necessary.
    if (!(i_hLeft > dryTol && i_hRight > dryTol))
    {
        if (i_hLeft < dryTol && i_hRight < dryTol)
            return;
        else if (i_hLeft < dryTol)
        {
            i_hLeft  = i_hRight;
            i_huLeft = -i_huRight;
            i_bLeft  = i_bRight;
        }
        else
        { //( i_hRight < dryTol )
            i_hRight  = i_hLeft;
            i_huRight = -i_huLeft;
            i_bRight  = i_bLeft;
        }
    }

    //! velocity on the left side of the edge
    T uLeft = i_huLeft / i_hLeft;
    //! velocity on the right side of the edge
    T uRight = i_huRight / i_hRight;

    //! wave speeds of the f-waves
    T waveSpeeds[2];

    // compute the wave speeds
    fWaveComputeWaveSpeeds(
        i_gravity,
        i_hLeft,
        i_hRight,
        i_huLeft,
        i_huRight,
        uLeft,
        uRight,
        i_bLeft,
        i_bRight,

        waveSpeeds
    );

    //! where to store the two f-waves
    T fWaves[2][2];

    // compute the decomposition into f-waves
    fWaveComputeWaveDecomposition(
        i_gravity,
        i_hLeft,
        i_hRight,
        i_huLeft,
        i_huRight,
        uLeft,
        uRight,
        i_bLeft,
        i_bRight,

        waveSpeeds,
        fWaves
    );

    // compute the net-updates
    // 1st wave family
    if (waveSpeeds[0] < -zeroTol)
    {                                    // left going
        o_netUpdates[0] += fWaves[0][0]; // hUpdateLeft
        o_netUpdates[2] += fWaves[0][1]; // huUpdateLeft
    }
    else if (waveSpeeds[0] > zeroTol)
    {                                    // right going
        o_netUpdates[1] += fWaves[0][0]; // hUpdateRight
        o_netUpdates[3] += fWaves[0][1]; // huUpdateRight
    }
    else
    {                                            // split waves
        o_netUpdates[0] += (T).5 * fWaves[0][0]; // hUpdateLeft
        o_netUpdates[2] += (T).5 * fWaves[0][1]; // huUpdateLeft
        o_netUpdates[1] += (T).5 * fWaves[0][0]; // hUpdateRight
        o_netUpdates[3] += (T).5 * fWaves[0][1]; // huUpdateRight
    }

    // 2nd wave family
    if (waveSpeeds[1] < -zeroTol)
    {                                    // left going
        o_netUpdates[0] += fWaves[1][0]; // hUpdateLeft
        o_netUpdates[2] += fWaves[1][1]; // huUpdateLeft
    }
    else if (waveSpeeds[1] > zeroTol)
    {                                    // right going
        o_netUpdates[1] += fWaves[1][0]; // hUpdateRight
        o_netUpdates[3] += fWaves[1][1]; // huUpdateRight
    }
    else
    {                                            // split waves
        o_netUpdates[0] += (T).5 * fWaves[1][0]; // hUpdateLeft
        o_netUpdates[2] += (T).5 * fWaves[1][1]; // huUpdateLeft
        o_netUpdates[1] += (T).5 * fWaves[1][0]; // hUpdateRight
        o_netUpdates[3] += (T).5 * fWaves[1][1]; // huUpdateRight
    }

    // compute maximum wave speed (-> CFL-condition)
    o_netUpdates[4] = fmax(fabs(waveSpeeds[0]), fabs(waveSpeeds[1]));
}

/**
 * Compute the edge local eigenvalues.
 *
 * @param i_gravity gravity constant.
 * @param i_hLeft height on the left side of the edge.
 * @param i_hRight height on the right side of the edge.
 * @param i_huLeft momentum on the left side of the edge.
 * @param i_huRight momentum on the right side of the edge.
 * @param i_uLeft velocity on the left side of the edge.
 * @param i_uRight velocity on the right side of the edge.
 * @param i_bLeft bathymetry on the left side of the edge.
 * @param i_bRight bathymetry on the right side of the edge.
 *
 * @param o_waveSpeeds will be set to: speeds of the linearized waves (eigenvalues).
 */
inline __device__ void fWaveComputeWaveSpeeds(
    const T i_gravity,
    const T i_hLeft,
    const T i_hRight,
    const T i_huLeft,
    const T i_huRight,
    const T i_uLeft,
    const T i_uRight,
    const T i_bLeft,
    const T i_bRight,

    T o_waveSpeeds[2]
)
{
    // compute eigenvalues of the jacobian matrices in states Q_{i-1} and Q_{i}
    T characteristicSpeeds[2];
    characteristicSpeeds[0] = i_uLeft - sqrt(i_gravity * i_hLeft);
    characteristicSpeeds[1] = i_uRight + sqrt(i_gravity * i_hRight);

    // compute "Roe speeds"
    T hRoe = (T).5 * (i_hRight + i_hLeft);
    T uRoe = i_uLeft * sqrt(i_hLeft) + i_uRight * sqrt(i_hRight);
    uRoe /= sqrt(i_hLeft) + sqrt(i_hRight);

    T roeSpeeds[2];
    roeSpeeds[0] = uRoe - sqrt(i_gravity * hRoe);
    roeSpeeds[1] = uRoe + sqrt(i_gravity * hRoe);

    // computer eindfeldt speeds
    T einfeldtSpeeds[2];
    einfeldtSpeeds[0] = fmin(characteristicSpeeds[0], roeSpeeds[0]);
    einfeldtSpeeds[1] = fmax(characteristicSpeeds[1], roeSpeeds[1]);

    // set wave speeds
    o_waveSpeeds[0] = einfeldtSpeeds[0];
    o_waveSpeeds[1] = einfeldtSpeeds[1];
}

/**
 * Compute the decomposition into f-Waves.
 *
 * @param i_gravity gravity constant.
 * @param i_hLeft height on the left side of the edge.
 * @param i_hRight height on the right side of the edge.
 * @param i_huLeft momentum on the left side of the edge.
 * @param i_huRight momentum on the right side of the edge.
 * @param i_uLeft velocity on the left side of the edge.
 * @param i_uRight velocity on the right side of the edge.
 * @param i_bLeft bathymetry on the left side of the edge.
 * @param i_bRight bathymetry on the right side of the edge.
 * @param i_waveSpeeds speeds of the linearized waves (eigenvalues).
 * @param o_fWaves  will be set to: Decomposition into f-Waves.
 */
inline __device__ void fWaveComputeWaveDecomposition(
    const T i_gravity,
    const T i_hLeft,
    const T i_hRight,
    const T i_huLeft,
    const T i_huRight,
    const T i_uLeft,
    const T i_uRight,
    const T i_bLeft,
    const T i_bRight,
    const T i_waveSpeeds[2],

    T o_fWaves[2][2]
)
{
    T lambdaDif = i_waveSpeeds[1] - i_waveSpeeds[0];

    // compute the inverse matrix R^{-1}
    T Rinv[2][2];

    T oneDivLambdaDif = (T)1. / lambdaDif;
    Rinv[0][0]        = oneDivLambdaDif * i_waveSpeeds[1];
    Rinv[0][1]        = -oneDivLambdaDif;

    Rinv[1][0] = oneDivLambdaDif * -i_waveSpeeds[0];
    Rinv[1][1] = oneDivLambdaDif;

    // right hand side
    T fDif[2];

    // calculate modified (bathymetry!) flux difference
    //  f(Q_i) - f(Q_{i-1})
    fDif[0] = i_huRight - i_huLeft;
    fDif[1] = i_huRight * i_uRight + (T).5 * i_gravity * i_hRight * i_hRight
              - (i_huLeft * i_uLeft + (T).5 * i_gravity * i_hLeft * i_hLeft);

    // \delta x \Psi[2]
    T psi = -i_gravity * (T).5 * (i_hRight + i_hLeft) * (i_bRight - i_bLeft);
    fDif[1] -= psi;

    // solve linear equations
    T beta[2];
    beta[0] = Rinv[0][0] * fDif[0] + Rinv[0][1] * fDif[1];
    beta[1] = Rinv[1][0] * fDif[0] + Rinv[1][1] * fDif[1];

    // return f-waves
    o_fWaves[0][0] = beta[0];
    o_fWaves[0][1] = beta[0] * i_waveSpeeds[0];

    o_fWaves[1][0] = beta[1];
    o_fWaves[1][1] = beta[1] * i_waveSpeeds[1];

    /*
     * in the case of a "zero strength" wave set the wave speeds to zero
     *   => non present waves don't affect the CFL-condition.
     * TODO: Seems not to affect the time step size
     */
    //  if( fmax(fabs(beta[0]), fabs(beta[1])) < zeroTol ) {
    //    io_waveSpeeds[0] = (T) 0.;
    //    io_waveSpeeds[1] = (T) 0.;
    //  }
}
} // namespace cuda
