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

/** error codes */
enum EError { NoError = 0, InvalidInput, PhysicalProblem, UnknownError };

/** allocate new MSSM model */
MSSMNoFV_onshell* gm2calc_mssmnofv_new(void);
/** delete MSSM model */
void gm2calc_mssmnofv_free(MSSMNoFV_onshell*);

/** set alpha_em(MZ) */
void gm2calc_mssmnofv_set_alpha_MZ(MSSMNoFV_onshell*, double);

/** set alpha_em(0) */
void gm2calc_mssmnofv_set_alpha_thompson(MSSMNoFV_onshell*, double);

/** set Ae */
void gm2calc_mssmnofv_set_Ae(MSSMNoFV_onshell*, unsigned, unsigned, double);
/** set Au */
void gm2calc_mssmnofv_set_Au(MSSMNoFV_onshell*, unsigned, unsigned, double);
/** set Ad */
void gm2calc_mssmnofv_set_Ad(MSSMNoFV_onshell*, unsigned, unsigned, double);
/** set gauge coupling g3 */
void gm2calc_mssmnofv_set_g3(MSSMNoFV_onshell*, double);
/** set bino mass */
void gm2calc_mssmnofv_set_MassB(MSSMNoFV_onshell*, double);
/** set wino mass */
void gm2calc_mssmnofv_set_MassWB(MSSMNoFV_onshell*, double);
/** set gluino mass */
void gm2calc_mssmnofv_set_MassG(MSSMNoFV_onshell*, double);
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
/** set Mu parameter */
void gm2calc_mssmnofv_set_Mu(MSSMNoFV_onshell*, double);
/** set tan(beta) */
void gm2calc_mssmnofv_set_TB(MSSMNoFV_onshell*, double);

/** set renormalization scale */
void gm2calc_mssmnofv_set_scale(MSSMNoFV_onshell*, double);

/** set CP-odd Higgs pole mass */
void gm2calc_mssmnofv_set_MAh_pole(MSSMNoFV_onshell*, double);

/** set Z boson pole mass */
void gm2calc_mssmnofv_set_MZ_pole(MSSMNoFV_onshell*, double);

/** set W boson pole mass */
void gm2calc_mssmnofv_set_MW_pole(MSSMNoFV_onshell*, double);

/** set top-Quark pole mass */
void gm2calc_mssmnofv_set_MFt_pole(MSSMNoFV_onshell*, double);

/** set running bottom-Quark mass mb(mb) */
void gm2calc_mssmnofv_set_MFb_running(MSSMNoFV_onshell*, double);

/** set tau-Lepton pole mass */
void gm2calc_mssmnofv_set_MFtau_pole(MSSMNoFV_onshell*, double);

/** set muon pole mass */
void gm2calc_mssmnofv_set_MFm_pole(MSSMNoFV_onshell*, double);

/** set smuon pole masses */
void gm2calc_mssmnofv_set_MSm_pole(MSSMNoFV_onshell*, unsigned, double);

/** set sneutrino pole mass */
void gm2calc_mssmnofv_set_MSvmL_pole(MSSMNoFV_onshell*, double);

/** set chargino pole mass */
void gm2calc_mssmnofv_set_MCha_pole(MSSMNoFV_onshell*, unsigned, double);

/** set neutralino pole mass */
void gm2calc_mssmnofv_set_MChi_pole(MSSMNoFV_onshell*, unsigned, double);

/** set verbose output */
void gm2calc_mssmnofv_set_verbose_output(MSSMNoFV_onshell*, int);

/** convert parameters to mixed on-shell/DR-bar scheme */
int gm2calc_mssmnofv_convert_to_onshell(MSSMNoFV_onshell*);

int gm2calc_mssmnofv_convert_to_onshell_params(
   MSSMNoFV_onshell*, double precision, unsigned max_iterations);

/** set calculate mass spectrum */
int gm2calc_mssmnofv_calculate_masses(MSSMNoFV_onshell*);

/** print model */
void print_mssmnofv(const MSSMNoFV_onshell*);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
