:Evaluate: BeginPackage["GM2Calc`"]

:Evaluate: loopOrder::usage =
    "Loop order of the calculation. Possible values: 0, 1, 2";

:Evaluate: tanBetaResummation::usage =
    "Enable/disable tan(beta) resummation. Possible values: True, False";

:Evaluate: forceOutput::usage =
    "Enforce output, even in case an error has occurred.";

:Evaluate: alphaMZ::usage =
    "Electromagnetic coupling in the MS-bar scheme at the scale MZ."

:Evaluate: alpha0::usage =
    "Electromagnetic coupling in the Thomson limit."

:Evaluate: alphaS::usage =
    "Strong coupling."

:Evaluate: MW::usage =
    "W-boson pole mass."

:Evaluate: MZ::usage =
    "Z-boson pole mass."

:Evaluate: MT::usage =
    "Top quark pole mass."

:Evaluate: mbmb::usage =
    "Bottom quark MS-bar mass mb at the scale mb."

:Evaluate: mbMZ::usage =
    "Bottom quark DR-bar mass at the scale MZ."

:Evaluate: ML::usage =
    "Tau lepton pole mass."

:Evaluate: MM::usage =
    "Muon pole mass."

:Evaluate: amu::usage =
    "Calculated value of the anomalous magnetic moment of the muon, amu = (g-2)/2.";

:Evaluate: Damu::usage =
    "Uncertainty of the calculated value of the anomalous magnetic moment of the muon.";

:Evaluate: UM::usage =
    "Mixing matrix of the negatively charged charginos.";

:Evaluate: UP::usage =
    "Mixing matrix of the positively charged charginos.";

:Evaluate: ZN::usage =
    "Mixing matrix of the neutralinos.";

:Evaluate: USt::usage =
    "Mixing matrix of the stops.";

:Evaluate: USb::usage =
    "Mixing matrix of the sbottoms.";

:Evaluate: UStau::usage =
    "Mixing matrix of the staus.";

:Evaluate: USm::usage =
    "Mixing matrix of the smuons.";

:Evaluate: MSveL::usage =
    "Electron-sneutrino mass.";

:Evaluate: MSvmL::usage =
    "Muon-sneutrino mass.";

:Evaluate: MSvtL::usage =
    "Tau-sneutrino mass.";

:Evaluate: MSe::usage =
    "Selectron masses.";

:Evaluate: MSm::usage =
    "Smuon masses.";

:Evaluate: MStau::usage =
    "Stau masses.";

:Evaluate: MSu::usage =
    "Sup masses.";

:Evaluate: MSd::usage =
    "Sdown masses.";

:Evaluate: MSc::usage =
    "Scharm masses.";

:Evaluate: MSs::usage =
    "Sstrange masses.";

:Evaluate: MSt::usage =
    "Stop masses.";

:Evaluate: MSb::usage =
    "Sbottom masses.";

:Evaluate: MChi::usage =
    "Neutralino masses.";

:Evaluate: MCha::usage =
    "Chargino masses.";

:Evaluate: MAh::usage =
    "CP-odd Higgs boson mass.";

:Evaluate: Mhh::usage =
    "CP-even Higgs bosons masses.";

:Evaluate: TB::usage =
    "tan(beta) = vu/vd.";

:Evaluate: Mu::usage =
    "Superpotential mu-parameter.";

:Evaluate: MassB::usage =
    "Soft-breaking bino mass parameter.";

:Evaluate: MassWB::usage =
    "Soft-breaking wino mass parameter.";

:Evaluate: MassG::usage =
    "Soft-breaking gluino mass parameter.";

:Evaluate: mq2::usage =
    "Soft-breaking squared left-handed squark mass parameters.";

:Evaluate: mu2::usage =
    "Soft-breaking squared right-handed up-type squark mass parameters.";

:Evaluate: md2::usage =
    "Soft-breaking squared right-handed down-type squark mass parameters.";

:Evaluate: ml2::usage =
    "Soft-breaking squared left-handed slepton mass parameters.";

:Evaluate: me2::usage =
    "Soft-breaking squared right-handed down-type slepton mass parameters.";

:Evaluate: Au::usage =
    "Soft-breaking trilinear couplings between Higgs bosons and up-type squarks.";

:Evaluate: Ad::usage =
    "Soft-breaking trilinear couplings between Higgs bosons and down-type squarks.";

:Evaluate: Ae::usage =
    "Soft-breaking trilinear couplings between Higgs bosons and down-type sleptons.";

:Evaluate: Q::usage =
    "Renormalization scale.";

:Evaluate: Yu::usage =
    "Yukawa couplings of up-type quarks.";

:Evaluate: Yd::usage =
    "Yukawa couplings of down-type quarks.";

:Evaluate: Ye::usage =
    "Yukawa couplings of down-type leptons.";

:Evaluate: GM2CalcSetFlags::usage =
    "GM2CalcSetFlags sets the configuration flags for GM2Calc.  Available flags are: {loopOrder, tanBetaResummation, forceOutput}.  Unset flags are set to their default values, see Options[GM2CalcSetFlags].  Use GM2CalcGetFlags[] to retrieve the flags currently set."

:Evaluate: GM2CalcGetFlags::usage =
    "GM2CalcGetFlags returns the current configuration flags for GM2Calc."

:Evaluate: GM2CalcSetSMParameters::usage =
    "GM2CalcSetSMParameters sets the Standard Model parameters input parameters.  Available Standard Model parameters are: {alphaMZ, alpha0, alphaS, MW, MZ, MT, mbmb, mbMZ, ML, MM}.  Unset parameters are set to their default values, see Options[GM2CalcSetSMParameters].  Use GM2CalcGetSMParameters[] to retrieve the current values of the Standard Model parameters."

:Evaluate: GM2CalcGetSMParameters::usage =
    "GM2CalcGetSMParameters returns the Standard Model parameters."

:Evaluate: GM2CalcAmuSLHAScheme::usage =
    "GM2CalcAmuSLHAScheme calculates amu and its uncertainty using the given SLHA parameters (SLHA interface).  Unset SLHA parameters are set to zero.  See Options[GM2CalcAmuSLHAScheme] for all SLHA parameters and their default values."

:Evaluate: GM2CalcAmuGM2CalcScheme::usage =
    "GM2CalcAmuGM2CalcScheme calculates amu and its uncertainty using the given parameters in the GM2Calc-specific renormalization scheme (GM2Calc interface).  Unset parameters are set to zero.  See Options[GM2CalcAmuGM2CalcScheme] for all input parameters in the GM2Calc scheme and their default values."

:Evaluate: GM2CalcAmuSLHAScheme::error = "`1`";
:Evaluate: GM2CalcAmuSLHAScheme::warning = "`1`";

:Evaluate: GM2CalcAmuGM2CalcScheme::error = "`1`";
:Evaluate: GM2CalcAmuGM2CalcScheme::warning = "`1`";

:Evaluate: Begin["`Private`"]

:Begin:
:Function: GM2CalcSetFlags
:Pattern: GM2CalcSetFlags[OptionsPattern[]]
:Arguments: {
   OptionValue[loopOrder],
   Boole[OptionValue[tanBetaResummation]],
   Boole[OptionValue[forceOutput]] }
:ArgumentTypes: {Integer, Integer, Integer}
:ReturnType: Integer
:End:

:Evaluate: Options[GM2CalcSetFlags] = {
   loopOrder -> 2,
   tanBetaResummation -> True,
   forceOutput -> False }

:Evaluate: GM2CalcSetFlags::wronglooporder = "Unsupported loop order: `1`";

:Begin:
:Function: GM2CalcGetFlags
:Pattern: GM2CalcGetFlags[]
:Arguments: {}
:ArgumentTypes: {}
:ReturnType: Manual
:End:

:Begin:
:Function: GM2CalcSetSMParameters
:Pattern: GM2CalcSetSMParameters[OptionsPattern[]]
:Arguments: {
   N @ OptionValue[alphaMZ],
   N @ OptionValue[alpha0],
   N @ OptionValue[alphaS],
   N @ OptionValue[MW],
   N @ OptionValue[MZ],
   N @ OptionValue[MT],
   N @ OptionValue[mbmb],
   N @ OptionValue[ML],
   N @ OptionValue[MM] }
:ArgumentTypes: { Real, Real, Real, Real, Real, Real, Real, Real, Real }
:ReturnType: Integer
:End:

:Evaluate: Options[GM2CalcSetSMParameters] = {
    alphaMZ -> 0.00775531,
    alpha0 -> 0.00729735,
    alphaS -> 0.1184,
    MW -> 80.385,
    MZ -> 91.1876,
    MT -> 173.34,
    mbmb -> 4.18,
    ML -> 1.777,
    MM -> 0.1056583715 }

:Begin:
:Function: GM2CalcGetSMParameters
:Pattern: GM2CalcGetSMParameters[]
:Arguments: {}
:ArgumentTypes: {}
:ReturnType: Manual
:End:

:Begin:
:Function: GM2CalcAmuSLHAScheme
:Pattern: GM2CalcAmuSLHAScheme[OptionsPattern[]]
:Arguments: {
    N @ OptionValue[MSvmL],
    N @ OptionValue[MSm][[1]],
    N @ OptionValue[MSm][[2]],
    N @ OptionValue[MChi][[1]],
    N @ OptionValue[MChi][[2]],
    N @ OptionValue[MChi][[3]],
    N @ OptionValue[MChi][[4]],
    N @ OptionValue[MCha][[1]],
    N @ OptionValue[MCha][[2]],
    N @ OptionValue[MAh],
    N @ OptionValue[TB],
    N @ OptionValue[Mu],
    N @ OptionValue[MassB],
    N @ OptionValue[MassWB],
    N @ OptionValue[MassG],
    N @ OptionValue[mq2][[1,1]],
    N @ OptionValue[mq2][[2,2]],
    N @ OptionValue[mq2][[3,3]],
    N @ OptionValue[ml2][[1,1]],
    N @ OptionValue[ml2][[2,2]],
    N @ OptionValue[ml2][[3,3]],
    N @ OptionValue[mu2][[1,1]],
    N @ OptionValue[mu2][[2,2]],
    N @ OptionValue[mu2][[3,3]],
    N @ OptionValue[md2][[1,1]],
    N @ OptionValue[md2][[2,2]],
    N @ OptionValue[md2][[3,3]],
    N @ OptionValue[me2][[1,1]],
    N @ OptionValue[me2][[2,2]],
    N @ OptionValue[me2][[3,3]],
    N @ OptionValue[Au][[3,3]],
    N @ OptionValue[Ad][[3,3]],
    N @ OptionValue[Ae][[2,2]],
    N @ OptionValue[Ae][[3,3]],
    N @ OptionValue[Q] }
:ArgumentTypes: {
   Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
   Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
   Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
   Real, Real, Real, Real, Real }
:ReturnType: Manual
:End:

:Evaluate: Options[GM2CalcAmuSLHAScheme] = {
    MSvmL     -> 0,
    MSm       -> {0, 0},
    MChi      -> {0, 0, 0, 0},
    MCha      -> {0, 0},
    MAh       -> 0,
    TB        -> 0,
    Mu        -> 0,
    MassB     -> 0,
    MassWB    -> 0,
    MassG     -> 0,
    mq2       -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    ml2       -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    mu2       -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    md2       -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    me2       -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    Au        -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    Ad        -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    Ae        -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    Q         -> 0 }

:Begin:
:Function: GM2CalcAmuGM2CalcScheme
:Pattern: GM2CalcAmuGM2CalcScheme[OptionsPattern[]]
:Arguments: {
   N @ OptionValue[MAh],
   N @ OptionValue[TB],
   N @ OptionValue[Mu],
   N @ OptionValue[MassB],
   N @ OptionValue[MassWB],
   N @ OptionValue[MassG],
   N @ OptionValue[mq2][[1,1]],
   N @ OptionValue[mq2][[2,2]],
   N @ OptionValue[mq2][[3,3]],
   N @ OptionValue[ml2][[1,1]],
   N @ OptionValue[ml2][[2,2]],
   N @ OptionValue[ml2][[3,3]],
   N @ OptionValue[mu2][[1,1]],
   N @ OptionValue[mu2][[2,2]],
   N @ OptionValue[mu2][[3,3]],
   N @ OptionValue[md2][[1,1]],
   N @ OptionValue[md2][[2,2]],
   N @ OptionValue[md2][[3,3]],
   N @ OptionValue[me2][[1,1]],
   N @ OptionValue[me2][[2,2]],
   N @ OptionValue[me2][[3,3]],
   N @ OptionValue[Au][[3,3]],
   N @ OptionValue[Ad][[3,3]],
   N @ OptionValue[Ae][[2,2]],
   N @ OptionValue[Ae][[3,3]],
   N @ OptionValue[Q] }
:ArgumentTypes: {
   Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
   Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
   Real, Real, Real, Real, Real, Real }
:ReturnType: Manual
:End:

:Evaluate: Options[GM2CalcAmuGM2CalcScheme] = {
    MAh    -> 0,
    TB     -> 0,
    Mu     -> 0,
    MassB  -> 0,
    MassWB -> 0,
    MassG  -> 0,
    mq2    -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    ml2    -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    mu2    -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    md2    -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    me2    -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    Au     -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    Ad     -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    Ae     -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    Q      -> 0 }

:Evaluate: End[]

:Evaluate: EndPackage[]


#include "gm2calc/gm2_1loop.h"
#include "gm2calc/gm2_2loop.h"
#include "gm2calc/gm2_uncertainty.h"
#include "gm2calc/gm2_version.h"
#include "gm2calc/MSSMNoFV_onshell.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mathlink.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NELEMS(x) (sizeof(x) / sizeof((x)[0]))

#define MLPutRule(link,s)                       \
 do {                                           \
    MLPutFunction(link, "Rule", 2);             \
    MLPutSymbol(link, s);                       \
 } while (0)

#define MLPutRuleToReal(link,v,s)               \
 do {                                           \
    MLPutRule(link, s);                         \
    MLPutReal(link, v);                         \
 } while (0)

#define MLPutRuleToInteger(link,v,s)            \
 do {                                           \
    MLPutRule(link, s);                         \
    MLPutInteger(link, v);                      \
 } while (0)

#define MLPutRuleToString(link,v,s)             \
 do {                                           \
    MLPutRule(link, s);                         \
    MLPutString(link, v);                       \
 } while (0)

#define MLPutRealMatrix(link,v,dim1,dim2)                  \
  do {                                                     \
    long dims[] = { dim1, dim2 };                          \
    MLPutDoubleArray(link, v, dims, NULL, NELEMS(dims));   \
  } while (0)

/* macros which access matrix/vector elements through GM2Calc interface
   functions */

#define MLPutRealVectorInterface(link,M,dim)            \
   do {                                                 \
      double M[dim];                                    \
      for (unsigned i = 0; i < dim; i++) {              \
         M[i] = gm2calc_mssmnofv_get_##M(model, i);     \
      }                                                 \
      MLPutRealList(link, M, dim);                      \
   } while (0)

#define MLPutRealMatrixInterface(link,M,dim1,dim2)           \
   do {                                                      \
      double M[dim1][dim2];                                  \
      for (unsigned i = 0; i < dim1; i++) {                  \
         for (unsigned k = 0; k < dim2; k++) {               \
            M[i][k] = gm2calc_mssmnofv_get_##M(model, i, k); \
         }                                                   \
      }                                                      \
      MLPutRealMatrix(link, (double*)M, dim1, dim2);         \
   } while (0)

#define MLPutComplexMatrixInterface(link,M,dim1,dim2)             \
   do {                                                           \
      MLPutFunction(link, "List", dim1);                          \
      for (unsigned i = 0; i < dim1; i++) {                       \
         MLPutFunction(link, "List", dim2);                       \
         for (unsigned k = 0; k < dim2; k++) {                    \
            double re = 0., im = 0.;                              \
            re = gm2calc_mssmnofv_get_##M(model, i, k, &im);      \
            MLPutComplex(link, re, im);                           \
         }                                                        \
      }                                                           \
   } while (0)

#define MLPutRuleToRealVectorInterface(link,M,name,dim) \
 do {                                                   \
    MLPutFunction(link, "Rule", 2);                     \
    MLPutSymbol(link, name);                            \
    MLPutRealVectorInterface(link,M,dim);               \
 } while (0)

#define MLPutRuleToRealMatrixInterface(link,v,name,dim1,dim2)   \
 do {                                                           \
    MLPutFunction(link, "Rule", 2);                             \
    MLPutSymbol(link, name);                                    \
    MLPutRealMatrixInterface(link,v,dim1,dim2);                 \
 } while (0)

#define MLPutRuleToComplexMatrixInterface(link,M,name,dim1,dim2)        \
 do {                                                                   \
    MLPutFunction(link, "Rule", 2);                                     \
    MLPutSymbol(link, name);                                            \
    MLPutComplexMatrixInterface(link,M,dim1,dim2);                      \
 } while (0)

/* global configuration flags */
struct Config_flags {
   int loopOrder;
   int tanBetaResummation;
   int forceOutput;
} config_flags = {
   .loopOrder = 2,
   .tanBetaResummation = 1,
   .forceOutput = 0
};

/* Standard Model parameters */
struct SM_parameters {
   double alphaMZ;
   double alpha0;
   double alphaS;
   double MW;
   double MZ;
   double MT;
   double mbmb;
   double ML;
   double MM;
} sm_parameters = {
   .alphaMZ = 0.00775531,
   .alpha0 = 0.00729735,
   .alphaS = 0.1184,
   .MW = 80.385,
   .MZ = 91.1876,
   .MT = 173.34,
   .mbmb = 4.18,
   .ML = 1.777,
   .MM = 0.1056583715
};

/* SUSY parameters for SLHA interface */
struct SLHA_parameters {
   double MSvmL;
   double MSm[2];
   double MChi[4];
   double MCha[2];
   double MAh;
   double TB;
   double Mu;
   double MassB;
   double MassWB;
   double MassG;
   double mq2[3][3];
   double ml2[3][3];
   double mu2[3][3];
   double md2[3][3];
   double me2[3][3];
   double Au[3][3];
   double Ad[3][3];
   double Ae[3][3];
   double Q;
};

/* SUSY parameters for GM2Calc interface */
struct GM2Calc_parameters {
   double MAh;
   double TB;
   double Mu;
   double MassB;
   double MassWB;
   double MassG;
   double mq2[3][3];
   double ml2[3][3];
   double mu2[3][3];
   double md2[3][3];
   double me2[3][3];
   double Au[3][3];
   double Ad[3][3];
   double Ae[3][3];
   double Q;
};

/******************************************************************/

void initialize_slha_parameters(struct SLHA_parameters* pars)
{
   pars->MSvmL     = 0.;
   memset(pars->MSm, 0, sizeof(pars->MSm[0]) * 2);
   memset(pars->MChi, 0, sizeof(pars->MChi[0]) * 4);
   memset(pars->MCha, 0, sizeof(pars->MCha[0]) * 2);
   pars->MAh       = 0.;
   pars->TB        = 0.;
   pars->Mu        = 0.;
   pars->MassB     = 0.;
   pars->MassWB    = 0.;
   pars->MassG     = 0.;
   memset(pars->mq2, 0, sizeof(pars->mq2[0][0]) * 3 * 3);
   memset(pars->ml2, 0, sizeof(pars->ml2[0][0]) * 3 * 3);
   memset(pars->mu2, 0, sizeof(pars->mu2[0][0]) * 3 * 3);
   memset(pars->md2, 0, sizeof(pars->md2[0][0]) * 3 * 3);
   memset(pars->me2, 0, sizeof(pars->me2[0][0]) * 3 * 3);
   memset(pars->Au, 0, sizeof(pars->Au[0][0]) * 3 * 3);
   memset(pars->Ad, 0, sizeof(pars->Ad[0][0]) * 3 * 3);
   memset(pars->Ae, 0, sizeof(pars->Ae[0][0]) * 3 * 3);
   pars->Q         = 0.;
}

/******************************************************************/

void initialize_gm2calc_parameters(struct GM2Calc_parameters* pars)
{
   pars->MAh     = 0.;
   pars->TB      = 0.;
   pars->Mu      = 0.;
   pars->MassB   = 0.;
   pars->MassWB  = 0.;
   pars->MassG   = 0.;
   memset(pars->mq2, 0, sizeof(pars->mq2[0][0]) * 3 * 3);
   memset(pars->ml2, 0, sizeof(pars->ml2[0][0]) * 3 * 3);
   memset(pars->mu2, 0, sizeof(pars->mu2[0][0]) * 3 * 3);
   memset(pars->md2, 0, sizeof(pars->md2[0][0]) * 3 * 3);
   memset(pars->me2, 0, sizeof(pars->me2[0][0]) * 3 * 3);
   memset(pars->Au, 0, sizeof(pars->Au[0][0]) * 3 * 3);
   memset(pars->Ad, 0, sizeof(pars->Ad[0][0]) * 3 * 3);
   memset(pars->Ae, 0, sizeof(pars->Ae[0][0]) * 3 * 3);
   pars->Q       = 0.;
}

/******************************************************************/

double sqr(double x) { return x*x; }

/******************************************************************/

void MLPutComplex(MLINK link, double re, double im)
{
   if (im == 0.) {
      MLPutReal(link, re);
   } else {
      MLPutFunction(link, "Complex", 2);
      MLPutReal(link, re);
      MLPutReal(link, im);
   }
}

/******************************************************************/

void put_error_message(const char* function_name,
                       const char* message_tag,
                       const char* message_str)
{
   MLPutFunction(stdlink, "CompoundExpression", 2);
   MLPutFunction(stdlink, "Message", 2);
   MLPutFunction(stdlink, "MessageName", 2);
   MLPutSymbol(stdlink, function_name);
   MLPutString(stdlink, message_tag);
   MLPutString(stdlink, message_str);
}

/******************************************************************/

gm2calc_error setup_model_slha_scheme(MSSMNoFV_onshell* model,
                                      const struct SLHA_parameters* pars)
{
   /* fill SM parameters */
   gm2calc_mssmnofv_set_alpha_MZ(model, sm_parameters.alphaMZ);
   gm2calc_mssmnofv_set_alpha_thompson(model, sm_parameters.alpha0);
   gm2calc_mssmnofv_set_g3(model, sqrt(4 * M_PI * sm_parameters.alphaS));
   gm2calc_mssmnofv_set_MT_pole(model, sm_parameters.MT);
   gm2calc_mssmnofv_set_MB_running(model, sm_parameters.mbmb);
   gm2calc_mssmnofv_set_MM_pole(model, sm_parameters.MM);
   gm2calc_mssmnofv_set_ML_pole(model, sm_parameters.ML);
   gm2calc_mssmnofv_set_MW_pole(model, sm_parameters.MW);
   gm2calc_mssmnofv_set_MZ_pole(model, sm_parameters.MZ);

   /* fill pole masses */
   gm2calc_mssmnofv_set_MSvmL_pole(model, pars->MSvmL);
   gm2calc_mssmnofv_set_MSm_pole(model, 0, pars->MSm[0]);
   gm2calc_mssmnofv_set_MSm_pole(model, 1, pars->MSm[1]);
   gm2calc_mssmnofv_set_MChi_pole(model, 0, pars->MChi[0]);
   gm2calc_mssmnofv_set_MChi_pole(model, 1, pars->MChi[1]);
   gm2calc_mssmnofv_set_MChi_pole(model, 2, pars->MChi[2]);
   gm2calc_mssmnofv_set_MChi_pole(model, 3, pars->MChi[3]);
   gm2calc_mssmnofv_set_MCha_pole(model, 0, pars->MCha[0]);
   gm2calc_mssmnofv_set_MCha_pole(model, 1, pars->MCha[1]);
   gm2calc_mssmnofv_set_MAh_pole(model, pars->MAh);

   /* fill DR-bar parameters */
   gm2calc_mssmnofv_set_TB(model, pars->TB);
   gm2calc_mssmnofv_set_Mu(model, pars->Mu);
   gm2calc_mssmnofv_set_MassB(model, pars->MassB);
   gm2calc_mssmnofv_set_MassWB(model, pars->MassWB);
   gm2calc_mssmnofv_set_MassG(model, pars->MassG);
   for (unsigned i = 0; i < 3; i++) {
      for (unsigned k = 0; k < 3; k++) {
         gm2calc_mssmnofv_set_mq2(model, i, k, pars->mq2[i][k]);
         gm2calc_mssmnofv_set_ml2(model, i, k, pars->ml2[i][k]);
         gm2calc_mssmnofv_set_mu2(model, i, k, pars->mu2[i][k]);
         gm2calc_mssmnofv_set_md2(model, i, k, pars->md2[i][k]);
         gm2calc_mssmnofv_set_me2(model, i, k, pars->me2[i][k]);
         gm2calc_mssmnofv_set_Au(model, i, k, pars->Au[i][k]);
         gm2calc_mssmnofv_set_Ad(model, i, k, pars->Ad[i][k]);
         gm2calc_mssmnofv_set_Ae(model, i, k, pars->Ae[i][k]);
      }
   }
   gm2calc_mssmnofv_set_scale(model, pars->Q);

   /* convert DR-bar parameters to on-shell */
   const gm2calc_error error = gm2calc_mssmnofv_convert_to_onshell(model);

   return error;
}

/******************************************************************/

gm2calc_error setup_model_gm2calc_scheme(MSSMNoFV_onshell* model,
                                         const struct GM2Calc_parameters* pars)
{
   /* fill SM parameters */
   gm2calc_mssmnofv_set_alpha_MZ(model, sm_parameters.alphaMZ);
   gm2calc_mssmnofv_set_alpha_thompson(model, sm_parameters.alpha0);
   gm2calc_mssmnofv_set_g3(model, sqrt(4 * M_PI * sm_parameters.alphaS));
   gm2calc_mssmnofv_set_MT_pole(model, sm_parameters.MT);
   gm2calc_mssmnofv_set_MB_running(model, sm_parameters.mbmb);
   gm2calc_mssmnofv_set_MM_pole(model, sm_parameters.MM);
   gm2calc_mssmnofv_set_ML_pole(model, sm_parameters.ML);
   gm2calc_mssmnofv_set_MW_pole(model, sm_parameters.MW);
   gm2calc_mssmnofv_set_MZ_pole(model, sm_parameters.MZ);

   /* fill DR-bar parameters */
   gm2calc_mssmnofv_set_TB(model, pars->TB);

   /* fill on-shell parameters */
   gm2calc_mssmnofv_set_Mu(model, pars->Mu);
   gm2calc_mssmnofv_set_MassB(model, pars->MassB);
   gm2calc_mssmnofv_set_MassWB(model, pars->MassWB);
   gm2calc_mssmnofv_set_MassG(model, pars->MassG);
   gm2calc_mssmnofv_set_MAh_pole(model, pars->MAh);
   gm2calc_mssmnofv_set_scale(model, pars->Q);

   /* fill DR-bar parameters */
   for (unsigned i = 0; i < 3; i++) {
      for (unsigned k = 0; k < 3; k++) {
         gm2calc_mssmnofv_set_mq2(model, i, k, pars->mq2[i][k]);
         gm2calc_mssmnofv_set_ml2(model, i, k, pars->ml2[i][k]);
         gm2calc_mssmnofv_set_mu2(model, i, k, pars->mu2[i][k]);
         gm2calc_mssmnofv_set_md2(model, i, k, pars->md2[i][k]);
         gm2calc_mssmnofv_set_me2(model, i, k, pars->me2[i][k]);
         gm2calc_mssmnofv_set_Au(model, i, k, pars->Au[i][k]);
         gm2calc_mssmnofv_set_Ad(model, i, k, pars->Ad[i][k]);
         gm2calc_mssmnofv_set_Ae(model, i, k, pars->Ae[i][k]);
      }
   }

   /* convert DR-bar parameters to on-shell */
   const gm2calc_error error = gm2calc_mssmnofv_calculate_masses(model);

   return error;
}

/******************************************************************/

double calculate_amu(MSSMNoFV_onshell* model)
{
   double amu = 0.;

   if (config_flags.loopOrder == 1) {
      if (config_flags.tanBetaResummation) {
         amu = gm2calc_mssmnofv_calculate_amu_1loop(model);
      } else {
         amu = gm2calc_mssmnofv_calculate_amu_1loop_non_tan_beta_resummed(model);
      }
   } else if (config_flags.loopOrder > 1) {
      if (config_flags.tanBetaResummation) {
         amu = gm2calc_mssmnofv_calculate_amu_1loop(model)
             + gm2calc_mssmnofv_calculate_amu_2loop(model);
      } else {
         amu = gm2calc_mssmnofv_calculate_amu_1loop_non_tan_beta_resummed(model)
             + gm2calc_mssmnofv_calculate_amu_2loop_non_tan_beta_resummed(model);
      }
   }

   return amu;
}

/******************************************************************/

double calculate_uncertainty(MSSMNoFV_onshell* model)
{
   double delta_amu = 0.;

   if (config_flags.loopOrder == 0) {
      delta_amu = gm2calc_mssmnofv_calculate_uncertainty_amu_0loop(model);
   } else if (config_flags.loopOrder == 1) {
      delta_amu = gm2calc_mssmnofv_calculate_uncertainty_amu_1loop(model);
   } else if (config_flags.loopOrder > 1) {
      delta_amu = gm2calc_mssmnofv_calculate_uncertainty_amu_2loop(model);
   }

   return delta_amu;
}

/******************************************************************/

static void print_package()
{
   static int do_print = 1;

   if (do_print) {
      printf("=======================================================\n");
      printf("GM2Calc " GM2CALC_VERSION "\n");
      printf("P. Athron, M. Bach, H. G. Fargnoli, C. Gnendiger,\n");
      printf("R. Greifenhagen, J.-h. Park, S. Paßehr, D. Stöckinger,\n");
      printf("H. Stöckinger-Kim, A. Voigt\n");
      printf("http://gm2calc.hepforge.org\n");
      printf("=======================================================\n");

      do_print = 0;
   }
}

/******************************************************************/

int GM2CalcSetFlags(int loopOrder_, int tanBetaResummation_, int forceOutput_)
{
   char loop_order_str[12];

   print_package();

   if (loopOrder_ < 0 || loopOrder_ > 2) {
      snprintf(loop_order_str, 12, "%d", loopOrder_);
      put_error_message("GM2CalcSetFlags", "wronglooporder", loop_order_str);
   }

   config_flags.loopOrder = loopOrder_;
   config_flags.tanBetaResummation = tanBetaResummation_;
   config_flags.forceOutput = forceOutput_;

   return 0;
}

/******************************************************************/

void GM2CalcGetFlags(void)
{
   MLPutFunction(stdlink, "List", 3);
   MLPutRuleToInteger(stdlink, config_flags.loopOrder, "loopOrder");
   MLPutRuleToString(stdlink, config_flags.tanBetaResummation ? "True" : "False", "tanBetaResummation");
   MLPutRuleToString(stdlink, config_flags.forceOutput ? "True" : "False", "forceOutput");
   MLEndPacket(stdlink);
}

/******************************************************************/

int GM2CalcSetSMParameters(
   double alphaMZ_,
   double alpha0_,
   double alphaS_,
   double MW_,
   double MZ_,
   double MT_,
   double mbmb_,
   double ML_,
   double MM_)
{
   sm_parameters.alphaMZ = alphaMZ_;
   sm_parameters.alpha0 = alpha0_;
   sm_parameters.alphaS = alphaS_;
   sm_parameters.MW = MW_;
   sm_parameters.MZ = MZ_;
   sm_parameters.MT = MT_;
   sm_parameters.mbmb = mbmb_;
   sm_parameters.ML = ML_;
   sm_parameters.MM = MM_;

   return 0;
}

/******************************************************************/

void GM2CalcGetSMParameters(void)
{
   MLPutFunction(stdlink, "List", 9);
   MLPutRuleToReal(stdlink, sm_parameters.alphaMZ, "alphaMZ");
   MLPutRuleToReal(stdlink, sm_parameters.alpha0, "alpha0");
   MLPutRuleToReal(stdlink, sm_parameters.alphaS, "alphaS");
   MLPutRuleToReal(stdlink, sm_parameters.MW, "MW");
   MLPutRuleToReal(stdlink, sm_parameters.MZ, "MZ");
   MLPutRuleToReal(stdlink, sm_parameters.MT, "MT");
   MLPutRuleToReal(stdlink, sm_parameters.mbmb, "mbmb");
   MLPutRuleToReal(stdlink, sm_parameters.ML, "ML");
   MLPutRuleToReal(stdlink, sm_parameters.MM, "MM");
   MLEndPacket(stdlink);
}

/******************************************************************/

void create_result_list(MSSMNoFV_onshell* model)
{
   double amu = calculate_amu(model);
   double Damu = calculate_uncertainty(model);

   MLPutFunction(stdlink, "List", 45);
   /* amu [2] */
   MLPutRuleToReal(stdlink, amu, "amu");
   MLPutRuleToReal(stdlink, Damu, "Damu");
   /* couplings [3] */
   MLPutRuleToReal(stdlink, sqr(gm2calc_mssmnofv_get_EL(model))/(4.*M_PI), "alphaMZ");
   MLPutRuleToReal(stdlink, sqr(gm2calc_mssmnofv_get_EL0(model))/(4.*M_PI), "alpha0");
   MLPutRuleToReal(stdlink, sqr(gm2calc_mssmnofv_get_g3(model))/(4.*M_PI), "alphaS");
   /* on-shell masses and parameters [40] */
   MLPutRuleToReal(stdlink, gm2calc_mssmnofv_get_MM(model), "MM");
   MLPutRuleToReal(stdlink, gm2calc_mssmnofv_get_MT(model), "MT");
   MLPutRuleToReal(stdlink, gm2calc_mssmnofv_get_MBMB(model), "mbmb");
   MLPutRuleToReal(stdlink, gm2calc_mssmnofv_get_MB(model), "mbMZ");
   MLPutRuleToReal(stdlink, gm2calc_mssmnofv_get_ML(model), "ML");
   MLPutRuleToReal(stdlink, gm2calc_mssmnofv_get_MW(model), "MW");
   MLPutRuleToReal(stdlink, gm2calc_mssmnofv_get_MZ(model), "MZ");
   MLPutRuleToRealVectorInterface(stdlink, MSm, "MSm", 2);
   MLPutRuleToRealMatrixInterface(stdlink, USm, "USm", 2, 2);
   MLPutRuleToReal(stdlink, gm2calc_mssmnofv_get_MSvmL(model), "MSvmL");
   MLPutRuleToRealVectorInterface(stdlink, MStau, "MStau", 2);
   MLPutRuleToRealMatrixInterface(stdlink, UStau, "UStau", 2, 2);
   MLPutRuleToRealVectorInterface(stdlink, MSt, "MSt", 2);
   MLPutRuleToRealMatrixInterface(stdlink, USt, "USt", 2, 2);
   MLPutRuleToRealVectorInterface(stdlink, MSb, "MSb", 2);
   MLPutRuleToRealMatrixInterface(stdlink, USb, "USb", 2, 2);
   MLPutRuleToRealVectorInterface(stdlink, MCha, "MCha", 2);
   MLPutRuleToComplexMatrixInterface(stdlink, UM, "UM", 2, 2);
   MLPutRuleToComplexMatrixInterface(stdlink, UP, "UP", 2, 2);
   MLPutRuleToRealVectorInterface(stdlink, MChi, "MChi", 4);
   MLPutRuleToComplexMatrixInterface(stdlink, ZN, "ZN", 4, 4);
   MLPutRuleToReal(stdlink, gm2calc_mssmnofv_get_MAh(model), "MAh");
   MLPutRuleToRealVectorInterface(stdlink, Mhh, "Mhh", 2);
   MLPutRuleToReal(stdlink, gm2calc_mssmnofv_get_TB(model), "TB");
   MLPutRuleToRealMatrixInterface(stdlink, Yu, "Yu", 3, 3);
   MLPutRuleToRealMatrixInterface(stdlink, Yd, "Yd", 3, 3);
   MLPutRuleToRealMatrixInterface(stdlink, Ye, "Ye", 3, 3);
   MLPutRuleToReal(stdlink, gm2calc_mssmnofv_get_Mu(model), "Mu");
   MLPutRuleToReal(stdlink, gm2calc_mssmnofv_get_MassB(model), "MassB");
   MLPutRuleToReal(stdlink, gm2calc_mssmnofv_get_MassWB(model), "MassWB");
   MLPutRuleToReal(stdlink, gm2calc_mssmnofv_get_MassG(model), "MassG");
   MLPutRuleToRealMatrixInterface(stdlink, mq2, "mq2", 3, 3);
   MLPutRuleToRealMatrixInterface(stdlink, ml2, "ml2", 3, 3);
   MLPutRuleToRealMatrixInterface(stdlink, mu2, "mu2", 3, 3);
   MLPutRuleToRealMatrixInterface(stdlink, md2, "md2", 3, 3);
   MLPutRuleToRealMatrixInterface(stdlink, me2, "me2", 3, 3);
   MLPutRuleToRealMatrixInterface(stdlink, Au, "Au", 3, 3);
   MLPutRuleToRealMatrixInterface(stdlink, Ad, "Ad", 3, 3);
   MLPutRuleToRealMatrixInterface(stdlink, Ae, "Ae", 3, 3);
   MLPutRuleToReal(stdlink, gm2calc_mssmnofv_get_scale(model), "Q");
   MLEndPacket(stdlink);
}

/******************************************************************/

void create_error_output(void)
{
   MLPutFunction(stdlink, "List", 0);
   MLEndPacket(stdlink);
}

/******************************************************************/

void GM2CalcAmuSLHAScheme(
   double MSvmL_,
   double MSm_1_,
   double MSm_2_,
   double MChi_1_,
   double MChi_2_,
   double MChi_3_,
   double MChi_4_,
   double MCha_1_,
   double MCha_2_,
   double MAh_,
   double TB_,
   double Mu_,
   double MassB_,
   double MassWB_,
   double MassG_,
   double mq2_11_,
   double mq2_22_,
   double mq2_33_,
   double ml2_11_,
   double ml2_22_,
   double ml2_33_,
   double mu2_11_,
   double mu2_22_,
   double mu2_33_,
   double md2_11_,
   double md2_22_,
   double md2_33_,
   double me2_11_,
   double me2_22_,
   double me2_33_,
   double Au_33_,
   double Ad_33_,
   double Ae_22_,
   double Ae_33_,
   double Q_
)
{
   struct SLHA_parameters pars;
   initialize_slha_parameters(&pars);

   pars.MSvmL     = MSvmL_;
   pars.MSm[0]    = MSm_1_;
   pars.MSm[1]    = MSm_2_;
   pars.MChi[0]   = MChi_1_;
   pars.MChi[1]   = MChi_2_;
   pars.MChi[2]   = MChi_3_;
   pars.MChi[3]   = MChi_4_;
   pars.MCha[0]   = MCha_1_;
   pars.MCha[1]   = MCha_2_;
   pars.MAh       = MAh_;
   pars.TB        = TB_;
   pars.Mu        = Mu_;
   pars.MassB     = MassB_;
   pars.MassWB    = MassWB_;
   pars.MassG     = MassG_;
   pars.mq2[0][0] = mq2_11_;
   pars.mq2[1][1] = mq2_22_;
   pars.mq2[2][2] = mq2_33_;
   pars.ml2[0][0] = ml2_11_;
   pars.ml2[1][1] = ml2_22_;
   pars.ml2[2][2] = ml2_33_;
   pars.mu2[0][0] = mu2_11_;
   pars.mu2[1][1] = mu2_22_;
   pars.mu2[2][2] = mu2_33_;
   pars.md2[0][0] = md2_11_;
   pars.md2[1][1] = md2_22_;
   pars.md2[2][2] = md2_33_;
   pars.me2[0][0] = me2_11_;
   pars.me2[1][1] = me2_22_;
   pars.me2[2][2] = me2_33_;
   pars.Au[2][2]  = Au_33_;
   pars.Ad[2][2]  = Ad_33_;
   pars.Ae[1][1]  = Ae_22_;
   pars.Ae[2][2]  = Ae_33_;
   pars.Q         = Q_;

   MSSMNoFV_onshell* model = gm2calc_mssmnofv_new();
   const gm2calc_error error = setup_model_slha_scheme(model, &pars);

   if (gm2calc_mssmnofv_have_warning(model)) {
      char msg[400];
      gm2calc_mssmnofv_get_warnings(model, msg, sizeof(msg));
      put_error_message("GM2CalcAmuSLHAScheme", "warning", msg);
   }

   if (error != gm2calc_NoError) {
      put_error_message("GM2CalcAmuSLHAScheme", "error",
                        gm2calc_error_str(error));
   }

   if (gm2calc_mssmnofv_have_problem(model)) {
      char msg[400];
      gm2calc_mssmnofv_get_problems(model, msg, sizeof(msg));
      put_error_message("GM2CalcAmuSLHAScheme", "error", msg);
   }

   if (error == gm2calc_NoError || config_flags.forceOutput) {
      create_result_list(model);
   } else {
      create_error_output();
   }

   gm2calc_mssmnofv_free(model);
}

/******************************************************************/

void GM2CalcAmuGM2CalcScheme(
   double MAh_,
   double TB_,
   double Mu_,
   double MassB_,
   double MassWB_,
   double MassG_,
   double mq2_11_,
   double mq2_22_,
   double mq2_33_,
   double ml2_11_,
   double ml2_22_,
   double ml2_33_,
   double mu2_11_,
   double mu2_22_,
   double mu2_33_,
   double md2_11_,
   double md2_22_,
   double md2_33_,
   double me2_11_,
   double me2_22_,
   double me2_33_,
   double Au_33_,
   double Ad_33_,
   double Ae_22_,
   double Ae_33_,
   double Q_)
{
   struct GM2Calc_parameters pars;
   initialize_gm2calc_parameters(&pars);

   pars.MAh       = MAh_;
   pars.TB        = TB_;
   pars.Mu        = Mu_;
   pars.MassB     = MassB_;
   pars.MassWB    = MassWB_;
   pars.MassG     = MassG_;
   pars.mq2[0][0] = mq2_11_;
   pars.mq2[1][1] = mq2_22_;
   pars.mq2[2][2] = mq2_33_;
   pars.ml2[0][0] = ml2_11_;
   pars.ml2[1][1] = ml2_22_;
   pars.ml2[2][2] = ml2_33_;
   pars.mu2[0][0] = mu2_11_;
   pars.mu2[1][1] = mu2_22_;
   pars.mu2[2][2] = mu2_33_;
   pars.md2[0][0] = md2_11_;
   pars.md2[1][1] = md2_22_;
   pars.md2[2][2] = md2_33_;
   pars.me2[0][0] = me2_11_;
   pars.me2[1][1] = me2_22_;
   pars.me2[2][2] = me2_33_;
   pars.Au[2][2]  = Au_33_;
   pars.Ad[2][2]  = Ad_33_;
   pars.Ae[1][1]  = Ae_22_;
   pars.Ae[2][2]  = Ae_33_;
   pars.Q         = Q_;

   MSSMNoFV_onshell* model = gm2calc_mssmnofv_new();
   const gm2calc_error error = setup_model_gm2calc_scheme(model, &pars);

   if (gm2calc_mssmnofv_have_warning(model)) {
      char msg[400];
      gm2calc_mssmnofv_get_warnings(model, msg, sizeof(msg));
      put_error_message("GM2CalcAmuGM2CalcScheme", "warning", msg);
   }

   if (error != gm2calc_NoError) {
      put_error_message("GM2CalcAmuGM2CalcScheme", "error",
                        gm2calc_error_str(error));
   }

   if (gm2calc_mssmnofv_have_problem(model)) {
      char msg[400];
      gm2calc_mssmnofv_get_problems(model, msg, sizeof(msg));
      put_error_message("GM2CalcAmuGM2CalcScheme", "error", msg);
   }

   if (error == gm2calc_NoError || config_flags.forceOutput) {
      create_result_list(model);
   } else {
      create_error_output();
   }

   gm2calc_mssmnofv_free(model);
}

/******************************************************************/

int main(int argc, char *argv[])
{
   return MLMain(argc, argv);
}
