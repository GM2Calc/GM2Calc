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

/** general THDM general basis input */
struct GeneralTHDM_general_basis;
typedef struct GeneralTHDM_general_basis GeneralTHDM_general_basis;

/** general THDM physical basis input */
struct GeneralTHDM_physical_basis;
typedef struct GeneralTHDM_physical_basis GeneralTHDM_physical_basis;

struct gm2calc_SM;
typedef struct gm2calc_SM gm2calc_SM;

/** Yukawa schemes */
typedef enum {
   GeneralTHDM_type_1,
   GeneralTHDM_type_2,
   GeneralTHDM_type_X,
   GeneralTHDM_type_Y,
   GeneralTHDM_general
} GeneralTHDM_yukawa_scheme;

/* ********** general basis ********** */

/** allocate new general basis input for the general THDM */
GeneralTHDM_general_basis* gm2calc_generalthdm_general_basis_new();

/** delete general basis */
void gm2calc_generalthdm_general_basis_free(GeneralTHDM_general_basis*);

void gm2calc_generalthdm_general_basis_set_yukawa_scheme(GeneralTHDM_general_basis*, GeneralTHDM_yukawa_scheme);
void gm2calc_generalthdm_general_basis_set_lambda_1(GeneralTHDM_general_basis*, double);
void gm2calc_generalthdm_general_basis_set_lambda_2(GeneralTHDM_general_basis*, double);
void gm2calc_generalthdm_general_basis_set_lambda_3(GeneralTHDM_general_basis*, double);
void gm2calc_generalthdm_general_basis_set_lambda_4(GeneralTHDM_general_basis*, double);
void gm2calc_generalthdm_general_basis_set_lambda_5(GeneralTHDM_general_basis*, double);
void gm2calc_generalthdm_general_basis_set_lambda_6(GeneralTHDM_general_basis*, double);
void gm2calc_generalthdm_general_basis_set_lambda_7(GeneralTHDM_general_basis*, double);
void gm2calc_generalthdm_general_basis_set_tan_beta(GeneralTHDM_general_basis*, double);
void gm2calc_generalthdm_general_basis_set_m122(GeneralTHDM_general_basis*, double);
void gm2calc_generalthdm_general_basis_set_Xu(GeneralTHDM_general_basis*, unsigned, unsigned, double);
void gm2calc_generalthdm_general_basis_set_Xd(GeneralTHDM_general_basis*, unsigned, unsigned, double);
void gm2calc_generalthdm_general_basis_set_Xe(GeneralTHDM_general_basis*, unsigned, unsigned, double);

/* ********** physical basis ********** */

/** allocate new physical basis input for the general THDM */
GeneralTHDM_physical_basis* gm2calc_generalthdm_physical_basis_new();

/** delete physical basis */
void gm2calc_generalthdm_physical_basis_free(GeneralTHDM_physical_basis*);

void gm2calc_generalthdm_physical_basis_set_yukawa_scheme(GeneralTHDM_physical_basis*, GeneralTHDM_yukawa_scheme);
void gm2calc_generalthdm_physical_basis_set_mh(GeneralTHDM_physical_basis*, double);
void gm2calc_generalthdm_physical_basis_set_mH(GeneralTHDM_physical_basis*, double);
void gm2calc_generalthdm_physical_basis_set_mA(GeneralTHDM_physical_basis*, double);
void gm2calc_generalthdm_physical_basis_set_mHp(GeneralTHDM_physical_basis*, double);
void gm2calc_generalthdm_physical_basis_set_sin_beta_minus_alpha(GeneralTHDM_physical_basis*, double);
void gm2calc_generalthdm_physical_basis_set_lambda_6(GeneralTHDM_physical_basis*, double);
void gm2calc_generalthdm_physical_basis_set_lambda_7(GeneralTHDM_physical_basis*, double);
void gm2calc_generalthdm_physical_basis_set_tan_beta(GeneralTHDM_physical_basis*, double);
void gm2calc_generalthdm_physical_basis_set_m122(GeneralTHDM_physical_basis*, double);
void gm2calc_generalthdm_physical_basis_set_Xu(GeneralTHDM_physical_basis*, unsigned, unsigned, double);
void gm2calc_generalthdm_physical_basis_set_Xd(GeneralTHDM_physical_basis*, unsigned, unsigned, double);
void gm2calc_generalthdm_physical_basis_set_Xe(GeneralTHDM_physical_basis*, unsigned, unsigned, double);

/* ********** THDM model ********** */

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
