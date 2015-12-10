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
struct C_MSSMNoFV_onshell;
typedef struct C_MSSMNoFV_onshell C_MSSMNoFV_onshell;

/** allocate new MSSM model */
C_MSSMNoFV_onshell* c_mssm_onshell_new(void);
/** delete MSSM model */
void c_mssm_onshell_free(C_MSSMNoFV_onshell*);

/** set tan(beta) */
void c_mssm_onshell_set_TB(C_MSSMNoFV_onshell*, double);
/** set Ae */
void c_mssm_onshell_set_Ae(C_MSSMNoFV_onshell*, unsigned, unsigned, double);
/** set Au */
void c_mssm_onshell_set_Au(C_MSSMNoFV_onshell*, unsigned, unsigned, double);
/** set Ad */
void c_mssm_onshell_set_Ad(C_MSSMNoFV_onshell*, unsigned, unsigned, double);
/** set Mu parameter */
void c_mssm_onshell_set_Mu(C_MSSMNoFV_onshell*, double);
/** set bino mass */
void c_mssm_onshell_set_MassB(C_MSSMNoFV_onshell*, double);
/** set wino mass */
void c_mssm_onshell_set_MassWB(C_MSSMNoFV_onshell*, double);
/** set gluino mass */
void c_mssm_onshell_set_MassG(C_MSSMNoFV_onshell*, double);
/** set MA0 */
void c_mssm_onshell_set_MA0(C_MSSMNoFV_onshell*, double);
/** set mq2 */
void c_mssm_onshell_set_mq2(C_MSSMNoFV_onshell*, unsigned, unsigned, double);
/** set mu2 */
void c_mssm_onshell_set_mu2(C_MSSMNoFV_onshell*, unsigned, unsigned, double);
/** set md2 */
void c_mssm_onshell_set_md2(C_MSSMNoFV_onshell*, unsigned, unsigned, double);
/** set ml2 */
void c_mssm_onshell_set_ml2(C_MSSMNoFV_onshell*, unsigned, unsigned, double);
/** set me2 */
void c_mssm_onshell_set_me2(C_MSSMNoFV_onshell*, unsigned, unsigned, double);
/** set renormalization scale */
void c_mssm_onshell_set_scale(C_MSSMNoFV_onshell*, double);

/** set calculate mass spectrum */
void c_mssm_onshell_calculate_masses(C_MSSMNoFV_onshell*);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
