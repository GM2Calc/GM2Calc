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

#ifndef GM2_MSSMNoFV_onshell_C_H
#define GM2_MSSMNoFV_onshell_C_H

#ifdef __cplusplus
extern "C" {
#endif

/** handle */
struct MSSMNoFV_onshell;
typedef struct MSSMNoFV_onshell MSSMNoFV_onshell;

/** allocate new MSSM model */
MSSMNoFV_onshell* gm2calc_mssmnofv_new(void);
/** delete MSSM model */
void gm2calc_mssmnofv_free(MSSMNoFV_onshell*);

/** set tan(beta) */
void gm2calc_mssmnofv_set_TB(MSSMNoFV_onshell*, double);
/** set Ae */
void gm2calc_mssmnofv_set_Ae(MSSMNoFV_onshell*, unsigned, unsigned, double);
/** set Au */
void gm2calc_mssmnofv_set_Au(MSSMNoFV_onshell*, unsigned, unsigned, double);
/** set Ad */
void gm2calc_mssmnofv_set_Ad(MSSMNoFV_onshell*, unsigned, unsigned, double);
/** set Mu parameter */
void gm2calc_mssmnofv_set_Mu(MSSMNoFV_onshell*, double);
/** set bino mass */
void gm2calc_mssmnofv_set_MassB(MSSMNoFV_onshell*, double);
/** set wino mass */
void gm2calc_mssmnofv_set_MassWB(MSSMNoFV_onshell*, double);
/** set gluino mass */
void gm2calc_mssmnofv_set_MassG(MSSMNoFV_onshell*, double);
/** set MA0 */
void gm2calc_mssmnofv_set_MA0(MSSMNoFV_onshell*, double);
/** set mq2 */
void gm2calc_mssmnofv_set_mq2(MSSMNoFV_onshell*, unsigned, unsigned, double);
/** set mu2 */
void gm2calc_mssmnofv_set_mu2(MSSMNoFV_onshell*, unsigned, unsigned, double);
/** set md2 */
void gm2calc_mssmnofv_set_md2(MSSMNoFV_onshell*, unsigned, unsigned, double);
/** set ml2 */
void gm2calc_mssmnofv_set_ml2(MSSMNoFV_onshell*, unsigned, unsigned, double);
/** set me2 */
void gm2calc_mssmnofv_set_me2(MSSMNoFV_onshell*, unsigned, unsigned, double);
/** set renormalization scale */
void gm2calc_mssmnofv_set_scale(MSSMNoFV_onshell*, double);

/** set calculate mass spectrum */
void gm2calc_mssmnofv_calculate_masses(MSSMNoFV_onshell*);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
