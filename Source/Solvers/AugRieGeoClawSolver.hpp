/** @file This file is part of the swe_solvers repository: https://github.com/TUM-I5/swe_solvers
 *
 * @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *
 * @section LICENSE
 * This software was developed at Technische Universitaet Muenchen, who is the owner of the software.
 *
 * According to good scientific practice, publications on results achieved in whole or in part due to this software
 * should cite at least one paper or referring to an URL presenting the this software software.
 *
 * The owner wishes to make the software available to all users to use, reproduce, modify, distribute and redistribute
 * also for commercial purposes under the following conditions of the original BSD license. Linking this software module
 * statically or dynamically with other modules is making a combined work based on this software. Thus, the terms and
 * conditions of this license cover the whole combination. As a special exception, the copyright holders of this
 * software give you permission to link it with independent modules or to instantiate templates and macros from this
 * software's source files to produce an executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice, provided that you also meet, for each
 * linked independent module, the terms and conditions of this license of that module.
 *
 * Copyright (c) 2012
 * Technische Universitaet Muenchen
 * Department of Informatics
 * Chair of Scientific Computing
 * http://www5.in.tum.de/
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the
 * following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this list of conditions and the following
 * disclaimer. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the
 * following disclaimer in the documentation and/or other materials provided with the distribution. All advertising
 * materials mentioning features or use of this software must display the following acknowledgement: This product
 * includes software developed by the Technische Universitaet Muenchen (TUM), Germany, and its contributors. Neither the
 * name of the Technische Universitaet Muenchen, Munich, Germany nor the names of its contributors may be used to
 * endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Binds GeoClaw's augmented solver to C.
 *
 * @section ACKNOWLEDGMENTS
 * Special thanks go to R.J. LeVeque and D.L. George for publishing their code.
 */

// number of f-waves in the Problem (2 shock linearization = 2, augmented = 3)
#ifndef NUMBER_OF_FWAVES
#define NUMBER_OF_FWAVES 3
#endif

/**
 * Extern declaration of the c_bing_geoclaw_riemann_aug_JCP routine (fixed arrays).
 *
 * @param i_maxNumberOfRiemannIterations maximum number Riemann iterations (solver might iterate over the Riemann
 * problem one day, currently fixed to 1)
 * @param i_variablesLeft array containing variables of the left cell: $(h_l, hu_l, b_l)^T$
 * @param i_variablesRight array containing variables of the right cell: $(h_l, hu_l, b_l)^T$
 * @param i_dryTol dry tolerance (definition of wet/dry cells).
 * @param i_g gravity constant (typically 9.81).
 * @param o_netUpdatesLeft will be set to: Net-update for the height(0)/normal momentum(1)/parallel momentum(2) of the
 * cell on the left side of the edge.
 * @param o_netUpdatesRight will be set to: Net-update for the height(0)/normal momentum(1)/parallel momentum(2) of the
 * cell on the right side of the edge.
 * @param o_waveSpeeds will be set to: (linearized) wave speeds -> Should be used in the CFL-condition.
 */
extern "C" void c_bind_geoclaw_riemann_aug_JCP(
  const int&    i_maxNumberOfRiemannIterations,
  const double  i_variablesLeft[3],
  const double  i_variablesRight[3],
  const double& i_dryTol,
  const double& i_g,
  double        o_netUpdatesLeft[3],
  double        o_netUpdatesRight[3],
  double        o_waveSpeeds[NUMBER_OF_FWAVES]
#if AUGMENTED_RIEMANN_EIGEN_COEFFICIENTS
  ,
  double o_eigenCoefficients[NUMBER_OF_FWAVES]
#endif
);

/**
 * Extern declaration of the c_bing_geoclaw_riemann_aug_JCP routine (pointers).
 *
 * @param i_maxNumberOfRiemannIterations maximum number Riemann iterations (solver might iterate over the Riemann
 * problem one day, currently fixed to 1)
 * @param i_variablesLeft array containing variables of the left cell: $(h_l, hu_l, b_l)^T$
 * @param i_variablesRight array containing variables of the right cell: $(h_l, hu_l, b_l)^T$
 * @param i_dryTol dry tolerance (definition of wet/dry cells).
 * @param i_g gravity constant (typically 9.81).
 * @param o_netUpdatesLeft will be set to: Net-update for the height(0)/normal momentum(1)/parallel momentum(2) of the
 * cell on the left side of the edge.
 * @param o_netUpdatesRight will be set to: Net-update for the height(0)/normal momentum(1)/parallel momentum(2) of the
 * cell on the right side of the edge.
 * @param o_waveSpeeds will be set to: (linearized) wave speeds -> Should be used in the CFL-condition.
 */
extern "C" void c_bind_geoclaw_riemann_aug_JCP(
  const int&    i_maxNumberOfRiemannIterations,
  const double* i_variablesLeft,
  const double* i_variablesRight,
  const double& i_dryTol,
  const double& i_g,
  double*       o_netUpdatesLeft,
  double*       o_netUpdatesRight,
  double*       o_waveSpeeds
#if AUGMENTED_RIEMANN_EIGEN_COEFFICIENTS
  ,
  double* o_eigenCoefficients
#endif
);
