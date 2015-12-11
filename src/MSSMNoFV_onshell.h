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

/**
 * @file MSSMNoFV_onshell.h
 * @brief contains declarations of C interface functions for the model
 *
 * This file contains the declarations for the C interface functions
 * used to set and retrieve the model parameters and masses.
 */

#ifdef __cplusplus
extern "C" {
#endif

/** handle */
struct MSSMNoFV_onshell;
typedef struct MSSMNoFV_onshell MSSMNoFV_onshell;

/** error codes */
enum EError { NoError = 0, InvalidInput, PhysicalProblem, UnknownError };

/** translate error codes into a string */
const char* gm2calc_error_str(enum EError);

/** allocate new MSSMNoFV model */
MSSMNoFV_onshell* gm2calc_mssmnofv_new(void);

/** delete MSSMNoFV model */
void gm2calc_mssmnofv_free(MSSMNoFV_onshell*);

/** set alpha_em(MZ) */
void gm2calc_mssmnofv_set_alpha_MZ(MSSMNoFV_onshell*, double);

/** set alpha_em(0) */
void gm2calc_mssmnofv_set_alpha_thompson(MSSMNoFV_onshell*, double);

/** set trilinear coupling Ae(i,k) */
void gm2calc_mssmnofv_set_Ae(MSSMNoFV_onshell*, unsigned, unsigned, double);
/** set trilinear coupling Au(i,k) */
void gm2calc_mssmnofv_set_Au(MSSMNoFV_onshell*, unsigned, unsigned, double);
/** set trilinear coupling Ad(i,k) */
void gm2calc_mssmnofv_set_Ad(MSSMNoFV_onshell*, unsigned, unsigned, double);
/** set gauge coupling g3 */
void gm2calc_mssmnofv_set_g3(MSSMNoFV_onshell*, double);
/** set bino mass */
void gm2calc_mssmnofv_set_MassB(MSSMNoFV_onshell*, double);
/** set wino mass */
void gm2calc_mssmnofv_set_MassWB(MSSMNoFV_onshell*, double);
/** set gluino mass */
void gm2calc_mssmnofv_set_MassG(MSSMNoFV_onshell*, double);
/** set soft-breaking squared mass parameter mq2(i,k) */
void gm2calc_mssmnofv_set_mq2(MSSMNoFV_onshell*, unsigned, unsigned, double);
/** set soft-breaking squared mass parameter mu2(i,k) */
void gm2calc_mssmnofv_set_mu2(MSSMNoFV_onshell*, unsigned, unsigned, double);
/** set soft-breaking squared mass parameter md2(i,k) */
void gm2calc_mssmnofv_set_md2(MSSMNoFV_onshell*, unsigned, unsigned, double);
/** set soft-breaking squared mass parameter ml2(i,k) */
void gm2calc_mssmnofv_set_ml2(MSSMNoFV_onshell*, unsigned, unsigned, double);
/** set soft-breaking squared mass parameter me2(i,k) */
void gm2calc_mssmnofv_set_me2(MSSMNoFV_onshell*, unsigned, unsigned, double);
/** set soft-breaking squared mass parameter Mu parameter */
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
void gm2calc_mssmnofv_set_MT_pole(MSSMNoFV_onshell*, double);

/** set running bottom-Quark mass mb(mb) */
void gm2calc_mssmnofv_set_MB_running(MSSMNoFV_onshell*, double);

/** set tau-Lepton pole mass */
void gm2calc_mssmnofv_set_ML_pole(MSSMNoFV_onshell*, double);

/** set muon pole mass */
void gm2calc_mssmnofv_set_MM_pole(MSSMNoFV_onshell*, double);

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

/** get Ae(i,k) */
double gm2calc_mssmnofv_get_Ae(const MSSMNoFV_onshell*, unsigned, unsigned);

/** get Ad(i,k) */
double gm2calc_mssmnofv_get_Ad(const MSSMNoFV_onshell*, unsigned, unsigned);

/** get Au(i,k) */
double gm2calc_mssmnofv_get_Au(const MSSMNoFV_onshell*, unsigned, unsigned);

/** get electromagnetic gauge coupling at MZ w/o hadronic corrections */
double gm2calc_mssmnofv_get_EL(const MSSMNoFV_onshell*);

/** get electromagnetic gauge coupling in Thomson limit */
double gm2calc_mssmnofv_get_EL0(const MSSMNoFV_onshell*);

/** get Hypercharge gauge coupling (not GUT normalized) */
double gm2calc_mssmnofv_get_gY(const MSSMNoFV_onshell*);

/** get Hypercharge gauge coupling (GUT normalized) */
double gm2calc_mssmnofv_get_g1(const MSSMNoFV_onshell*);

/** get left gauge coupling */
double gm2calc_mssmnofv_get_g2(const MSSMNoFV_onshell*);

/** get strong gauge coupling */
double gm2calc_mssmnofv_get_g3(const MSSMNoFV_onshell*);

/** get tan(beta) */
double gm2calc_mssmnofv_get_TB(const MSSMNoFV_onshell*);

/** get bino mass */
double gm2calc_mssmnofv_get_MassB(const MSSMNoFV_onshell*);

/** get wino mass */
double gm2calc_mssmnofv_get_MassWB(const MSSMNoFV_onshell*);

/** get gluino mass */
double gm2calc_mssmnofv_get_MassG(const MSSMNoFV_onshell*);

/** get Mu parameter */
double gm2calc_mssmnofv_get_Mu(const MSSMNoFV_onshell*);

/** get left-handed up-Squark soft-breaking squared mass */
double gm2calc_mssmnofv_get_mq2(const MSSMNoFV_onshell*, unsigned, unsigned);

/** get right-handed down-Squark soft-breaking squared mass */
double gm2calc_mssmnofv_get_md2(const MSSMNoFV_onshell*, unsigned, unsigned);

/** get right-handed up-Squark soft-breaking squared mass */
double gm2calc_mssmnofv_get_mu2(const MSSMNoFV_onshell*, unsigned, unsigned);

/** get left-handed up-Lepton soft-breaking squared mass */
double gm2calc_mssmnofv_get_ml2(const MSSMNoFV_onshell*, unsigned, unsigned);

/** get right-handed down-Lepton soft-breaking squared mass */
double gm2calc_mssmnofv_get_me2(const MSSMNoFV_onshell*, unsigned, unsigned);

/** get vacuum expectation value */
double gm2calc_mssmnofv_get_vev(const MSSMNoFV_onshell*);

/** get renormalization scale */
double gm2calc_mssmnofv_get_scale(const MSSMNoFV_onshell*);

/** get MW */
double gm2calc_mssmnofv_get_MW(const MSSMNoFV_onshell*);

/** get MZ */
double gm2calc_mssmnofv_get_MZ(const MSSMNoFV_onshell*);

/** get electron mass */
double gm2calc_mssmnofv_get_ME(const MSSMNoFV_onshell*);

/** get muon pole mass */
double gm2calc_mssmnofv_get_MM(const MSSMNoFV_onshell*);

/** get tau mass */
double gm2calc_mssmnofv_get_ML(const MSSMNoFV_onshell*);

/** get up-Quark mass */
double gm2calc_mssmnofv_get_MU(const MSSMNoFV_onshell*);

/** get charm-Quark mass */
double gm2calc_mssmnofv_get_MC(const MSSMNoFV_onshell*);

/** get top-Quark mass */
double gm2calc_mssmnofv_get_MT(const MSSMNoFV_onshell*);

/** get down-Quark mass */
double gm2calc_mssmnofv_get_MD(const MSSMNoFV_onshell*);

/** get strange-Quark mass */
double gm2calc_mssmnofv_get_MS(const MSSMNoFV_onshell*);

/** get bottom-Quark DR-bar mass mb(MZ) */
double gm2calc_mssmnofv_get_MB(const MSSMNoFV_onshell*);

/** get bottom-Quark MS-bar mass mb(mb) */
double gm2calc_mssmnofv_get_MBMB(const MSSMNoFV_onshell*);

/** get CP-odde Higgs pole mass */
double gm2calc_mssmnofv_get_MAh(const MSSMNoFV_onshell*);

/** get chargino pole mass */
double gm2calc_mssmnofv_get_MCha(const MSSMNoFV_onshell*, unsigned);

/** get neutralino pole mass */
double gm2calc_mssmnofv_get_MChi(const MSSMNoFV_onshell*, unsigned);

/** get smuon pole masses */
double gm2calc_mssmnofv_get_MSm(const MSSMNoFV_onshell*, unsigned);

/** get sneutrino pole mass */
double gm2calc_mssmnofv_get_MSvmL(const MSSMNoFV_onshell*);

/** get lepton Yukawa coupling */
double gm2calc_mssmnofv_get_Ye(const MSSMNoFV_onshell*, unsigned, unsigned);

/** get down-Quark Yukawa coupling */
double gm2calc_mssmnofv_get_Yd(const MSSMNoFV_onshell*, unsigned, unsigned);

/** get up-Quark Yukawa coupling */
double gm2calc_mssmnofv_get_Yu(const MSSMNoFV_onshell*, unsigned, unsigned);


/** convert parameters to mixed on-shell/DR-bar scheme */
int gm2calc_mssmnofv_convert_to_onshell(MSSMNoFV_onshell*);

/** convert parameters to mixed on-shell/DR-bar scheme */
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
