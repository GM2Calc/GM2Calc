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

#ifndef GM2_GeneralTHDM_H
#define GM2_GeneralTHDM_H

/**
 * @file GeneralTHDM.h
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
struct GeneralTHDM;
typedef struct GeneralTHDM GeneralTHDM;

/** Yukawa schemes */
typedef enum {
   GeneralTHDM_type_1,
   GeneralTHDM_type_2,
   GeneralTHDM_type_X,
   GeneralTHDM_type_Y,
   GeneralTHDM_general
} GeneralTHDM_yukawa_scheme;

/** general THDM general basis input */
struct GeneralTHDM_general_basis {
   GeneralTHDM_yukawa_scheme yukawa_scheme;
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
typedef struct GeneralTHDM_general_basis GeneralTHDM_general_basis;

/** general THDM physical basis input */
struct GeneralTHDM_physical_basis {
   GeneralTHDM_yukawa_scheme yukawa_scheme;
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
typedef struct GeneralTHDM_physical_basis GeneralTHDM_physical_basis;

struct gm2calc_SM;
typedef struct gm2calc_SM gm2calc_SM;

/** allocate new general THDM model with general basis input */
gm2calc_error gm2calc_generalthdm_new_with_general_basis(GeneralTHDM**, GeneralTHDM_general_basis*, gm2calc_SM*);

/** allocate new general THDM model with physical basis input */
gm2calc_error gm2calc_generalthdm_new_with_physical_basis(GeneralTHDM**, GeneralTHDM_physical_basis*, gm2calc_SM*);

/** delete general THDM model */
void gm2calc_generalthdm_free(GeneralTHDM*);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
