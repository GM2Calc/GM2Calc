/* ====================================================================
 * This file is part of GM2Calc.
 *
 * GM2Calc is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * GM2Calc is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GM2Calc.  If not, see
 * <http://www.gnu.org/licenses/>.
 * ==================================================================== */

#ifndef GM2_THDM_H
#define GM2_THDM_H

/**
 * @file THDM.h
 * @brief contains declarations of C interface functions for the model
 *
 * This file contains the declarations for the C interface functions
 * used to set and retrieve the model parameters and masses.
 */

#include "gm2calc/gm2_error.h"

#ifdef __cplusplus
extern "C" {
#endif

/** general THDM handle */
struct THDM;
typedef struct THDM THDM;

/** Yukawa schemes */
typedef enum {
   THDM_type_1,
   THDM_type_2,
   THDM_type_X,
   THDM_type_Y,
   THDM_general
} THDM_yukawa_scheme;

/** general THDM general basis input */
struct THDM_general_basis {
   THDM_yukawa_scheme yukawa_scheme;
   double lambda1;
   double lambda2;
   double lambda3;
   double lambda4;
   double lambda5;
   double lambda6;
   double lambda7;
   double tan_beta;
   double m122;
   double Xu[3][3];
   double Xd[3][3];
   double Xl[3][3];
};
typedef struct THDM_general_basis THDM_general_basis;

/** general THDM physical basis input */
struct THDM_physical_basis {
   THDM_yukawa_scheme yukawa_scheme;
   double mh;
   double mH;
   double mA;
   double mHp;
   double sin_beta_minus_alpha;
   double lambda6;
   double lambda7;
   double tan_beta;
   double m122;
   double Xu[3][3];
   double Xd[3][3];
   double Xl[3][3];
};
typedef struct THDM_physical_basis THDM_physical_basis;

struct SM;
typedef struct SM SM;

/** allocate new general THDM model with general basis input */
gm2calc_error gm2calc_thdm_new_with_general_basis(THDM**, THDM_general_basis*, SM*);

/** allocate new general THDM model with physical basis input */
gm2calc_error gm2calc_thdm_new_with_physical_basis(THDM**, THDM_physical_basis*, SM*);

/** delete general THDM model */
void gm2calc_thdm_free(THDM*);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
