#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"

#include "gm2calc/MSSMNoFV_onshell.hpp"

#include "gm2_config_options.hpp"
#include "gm2_numerics.hpp"
#include "gm2_slha_io.hpp"

#include <sstream>

#define CHECK_CLOSE(a,b,eps)                            \
   do {                                                 \
      CHECK((a) == doctest::Approx(b).epsilon(eps));    \
   } while (0)


TEST_CASE("fill_valid_config")
{
   gm2calc::Config_options config;

   char const * const slha_input = R"(
Block GM2CalcConfig
     0     1     # output format (0, 1, 2, 3, 4)
     1     1     # loop order (0, 1 or 2)
     2     0     # disable/enable tan(beta) resummation (0 or 1)
     3     1     # force output (0 or 1)
     4     1     # verbose output (0 or 1)
     5     1     # calculate uncertainty
)";

   std::istringstream stream(slha_input);

   gm2calc::GM2_slha_io slha;
   slha.read_from_stream(stream);
   slha.fill(config);

   CHECK(config.output_format == 1);
   CHECK(config.loop_order == 1);
   CHECK(config.tanb_resummation == false);
   CHECK(config.force_output == true);
   CHECK(config.verbose_output == true);
   CHECK(config.calculate_uncertainty == true);
}


// checks that default configuration values are used if invalid
// configuration values are passed
TEST_CASE("fill_invalid_config")
{
   const gm2calc::Config_options default_config;
   gm2calc::Config_options config;

   char const * const slha_input = R"(
Block GM2CalcConfig
     0     9     # output format (0, 1, 2, 3, 4)
     1     9     # loop order (0, 1 or 2)
     2     2     # disable/enable tan(beta) resummation (0 or 1)
     3     2     # force output (0 or 1)
     4     2     # verbose output (0 or 1)
     5     2     # calculate uncertainty (0 or 1)
)";

   std::istringstream stream(slha_input);

   gm2calc::GM2_slha_io slha;
   slha.read_from_stream(stream);
   slha.fill(config);

   CHECK(config.output_format         == default_config.output_format);
   CHECK(config.loop_order            == default_config.loop_order);
   CHECK(config.tanb_resummation      == default_config.tanb_resummation);
   CHECK(config.force_output          == default_config.force_output);
   CHECK(config.verbose_output        == default_config.verbose_output);
   CHECK(config.calculate_uncertainty == default_config.calculate_uncertainty);
}


TEST_CASE("fill_gm2calc")
{
   gm2calc::MSSMNoFV_onshell model;

   char const * const slha_input = R"(
Block GM2CalcInput
     0     8.66360379E+02   # ren. scale Q          [2L]
     1     0.00775531       # alpha(MZ)             [1L]
     2     0.00729735       # alpha(0)              [2L]
     3     10               # tan(beta) DR-bar at Q [1L]
     4     619.858          # Mu parameter on-shell [1L]
     5     211.722          # M1 on-shell           [1L]
     6     401.057          # M2 on-shell           [1L]
     7     1.10300877E+03   # M3                    [2L]
     8     707.025          # MA(pole)              [2L]
     9     3.51653258E+02   # msl(1,1)              [2L]
    10     356.09           # msl(2,2) on-shell     [1L]
    11     3.50674223E+02   # msl(3,3)              [2L]
    12     2.21215037E+02   # mse(1,1)              [2L]
    13     225.076          # mse(2,2) on-shell     [1L]
    14     2.18022142E+02   # mse(3,3)              [2L]
    15     1.00711403E+03   # msq(1,1)              [2L]
    16     1.00711149E+03   # msq(2,2)              [2L]
    17     9.29083096E+02   # msq(3,3)              [2L]
    18     9.69369660E+02   # msu(1,1)              [2L]
    19     9.69366965E+02   # msu(2,2)              [2L]
    20     7.99712943E+02   # msu(3,3)              [2L]
    21     9.64756473E+02   # msd(1,1)              [2L]
    22     9.64753818E+02   # msd(2,2)              [2L]
    23     9.60016201E+02   # msd(3,3)              [2L]
    24     0                # Ae(1,1)               [irrelevant]
    25     -2.93720212E+02  # Ae(2,2) DR-bar        [1L]
    26     -2.92154796E+02  # Ae(3,3)               [2L]
    27     0                # Ad(1,1)               [irrelevant]
    28     0                # Ad(2,2)               [irrelevant]
    29     -1.28330100E+03  # Ad(3,3)               [2L]
    30     0                # Au(1,1)               [irrelevant]
    31     0                # Au(2,2)               [irrelevant]
    32     -8.70714986E+02  # Au(3,3)               [2L]
Block SMINPUTS
     3     0.1185           # alpha_s(MZ) SM MS-bar [2L]
     4     91.1876          # M_Z(pole)             [1L]
     5     4.18             # m_b(m_b) SM MS-bar    [2L]
     6     173.34           # M_top(pole)           [2L]
     7     1.777            # M_tau(pole)           [2L]
     8     0.00000000E+00   # mnu3(pole)            [irrelevant]
     9     80.385           # M_W(pole)             [1L]
    11     0.000510998928   # melectron(pole)       [irrelevant]
    12     0.00000000E+00   # mnu1(pole)            [irrelevant]
    13     0.1056583715     # M_muon(pole)          [1L]
    14     0.00000000E+00   # mnu2(pole)            [irrelevant]
    21     4.76052706E-03   # md                    [irrelevant]
    22     2.40534062E-03   # mu                    [irrelevant]
    23     1.04230487E-01   # ms                    [irrelevant]
    24     1.27183378E+00   # mc                    [irrelevant]
)";

   std::istringstream stream(slha_input);

   gm2calc::GM2_slha_io slha;
   slha.read_from_stream(stream);
   slha.fill_gm2calc(model);

   const double eps = std::numeric_limits<double>::epsilon();
   const double Pi = 3.14159265358979323846;

   // Block GM2CalcInput
   CHECK_CLOSE(model.get_scale(), 8.66360379E+02, eps);
   CHECK_CLOSE(model.get_EL(), std::sqrt(4*Pi*0.00775531), eps);
   CHECK_CLOSE(model.get_EL0(), std::sqrt(4*Pi*0.00729735), eps);
   CHECK_CLOSE(model.get_TB(), 10.0, eps);
   CHECK_CLOSE(model.get_Mu(), 619.858, eps);
   CHECK_CLOSE(model.get_MassB(), 211.722, eps);
   CHECK_CLOSE(model.get_MassWB(), 401.057, eps);
   CHECK_CLOSE(model.get_MassG(), 1.10300877E+03, eps);
   CHECK_CLOSE(model.get_MA0(), 707.025, eps);
   CHECK_CLOSE(model.get_ml2(0,0), gm2calc::sqr(3.51653258E+02), eps);
   CHECK_CLOSE(model.get_ml2(1,1), gm2calc::sqr(356.09        ), eps);
   CHECK_CLOSE(model.get_ml2(2,2), gm2calc::sqr(3.50674223E+02), eps);
   CHECK_CLOSE(model.get_me2(0,0), gm2calc::sqr(2.21215037E+02), eps);
   CHECK_CLOSE(model.get_me2(1,1), gm2calc::sqr(225.076       ), eps);
   CHECK_CLOSE(model.get_me2(2,2), gm2calc::sqr(2.18022142E+02), eps);
   CHECK_CLOSE(model.get_mq2(0,0), gm2calc::sqr(1.00711403E+03), eps);
   CHECK_CLOSE(model.get_mq2(1,1), gm2calc::sqr(1.00711149E+03), eps);
   CHECK_CLOSE(model.get_mq2(2,2), gm2calc::sqr(9.29083096E+02), eps);
   CHECK_CLOSE(model.get_mu2(0,0), gm2calc::sqr(9.69369660E+02), eps);
   CHECK_CLOSE(model.get_mu2(1,1), gm2calc::sqr(9.69366965E+02), eps);
   CHECK_CLOSE(model.get_mu2(2,2), gm2calc::sqr(7.99712943E+02), eps);
   CHECK_CLOSE(model.get_md2(0,0), gm2calc::sqr(9.64756473E+02), eps);
   CHECK_CLOSE(model.get_md2(1,1), gm2calc::sqr(9.64753818E+02), eps);
   CHECK_CLOSE(model.get_md2(2,2), gm2calc::sqr(9.60016201E+02), eps);
   CHECK_CLOSE(model.get_Ae(0,0), 0              , eps);
   CHECK_CLOSE(model.get_Ae(1,1), -2.93720212E+02, eps);
   CHECK_CLOSE(model.get_Ae(2,2), -2.92154796E+02, eps);
   CHECK_CLOSE(model.get_Ad(0,0), 0              , eps);
   CHECK_CLOSE(model.get_Ad(1,1), 0              , eps);
   CHECK_CLOSE(model.get_Ad(2,2), -1.28330100E+03, eps);
   CHECK_CLOSE(model.get_Au(0,0), 0              , eps);
   CHECK_CLOSE(model.get_Au(1,1), 0              , eps);
   CHECK_CLOSE(model.get_Au(2,2), -8.70714986E+02, eps);

   // Block SMINPUTS
   CHECK_CLOSE(model.get_g3(), std::sqrt(4*Pi*0.1185), eps);
   CHECK_CLOSE(model.get_MZ(), 91.1876, eps);
   CHECK_CLOSE(model.get_MBMB(), 4.18, eps);
   CHECK_CLOSE(model.get_MT(), 173.34, eps);
   CHECK_CLOSE(model.get_ML(), 1.777, eps);
   CHECK_CLOSE(model.get_MW(), 80.385, eps);
   CHECK_CLOSE(model.get_MM(), 0.1056583715, eps);
   CHECK_CLOSE(model.get_MD(), 4.76052706E-03, eps);
   CHECK_CLOSE(model.get_MU(), 2.40534062E-03, eps);
   CHECK_CLOSE(model.get_MS(), 1.04230487E-01, eps);
   CHECK_CLOSE(model.get_MC(), 1.27183378E+00, eps);
}


TEST_CASE("fill_slha")
{
   gm2calc::MSSMNoFV_onshell model;

   char const * const slha_input = R"(
Block GM2CalcInput
     1     0.00775531       # alpha(MZ)             [1L]
     2     0.00729735       # alpha(0)              [2L]
Block SMINPUTS
     3     0.1185           # alpha_s(MZ) SM MSbar  [2L]
     4     9.11876000E+01   # MZ(pole)              [1L]
     5     4.18000000E+00   # mb(mb) SM MSbar       [2L]
     6     1.73340000E+02   # mtop(pole)            [2L]
     7     1.77700000E+00   # mtau(pole)            [2L]
     8     0.00000000E+00   # mnu3(pole)            [irrelevant]
     9     80.1             # MW(pole)              [1L]
    11     0.000510998928   # melectron(pole)       [irrelevant]
    12     0.00000000E+00   # mnu1(pole)            [irrelevant]
    13     0.1056583715     # mmuon(pole)           [1L]
    14     0.00000000E+00   # mnu2(pole)            [irrelevant]
    21     4.76052706E-03   # md                    [irrelevant]
    22     2.40534062E-03   # mu                    [irrelevant]
    23     1.04230487E-01   # ms                    [irrelevant]
    24     1.27183378E+00   # mc                    [irrelevant]
Block MASS                      # Mass spectrum
        24     8.03773317e+01   # W                 [1L, preferred value]
        25     1.25712136e+02   # h0                [irrelevant]
        35     1.50002058e+03   # H0                [irrelevant]
        36     1.50000000e+03   # A0                [2L]
        37     1.50241108e+03   # H+                [irrelevant]
   1000021     2.27707784e+03   # ~g                [irrelevant]
   1000022     2.01611468e+02   # neutralino(1)     [1L]
   1000023     4.10040273e+02   # neutralino(2)     [1L]
   1000024     4.09989890e+02   # chargino(1)       [1L]
   1000025    -5.16529941e+02   # neutralino(3)     [1L]
   1000035     5.45628749e+02   # neutralino(4)     [1L]
   1000037     5.46057190e+02   # chargino(2)       [1L]
   1000001     7.06303219e+03   # ~d_L              [irrelevant]
   1000002     7.06271372e+03   # ~u_L              [irrelevant]
   1000003     7.06305914e+03   # ~s_L              [irrelevant]
   1000004     7.06274067e+03   # ~c_L              [irrelevant]
   1000005     7.10826650e+03   # ~b_1              [irrelevant]
   1000006     7.16237728e+03   # ~t_1              [irrelevant]
   1000011     5.25167746e+02   # ~e_L              [irrelevant]
   1000012     5.18841083e+02   # ~nu_e_L           [irrelevant]
   1000013     5.25187016e+02   # ~mu_L             [1L]
   1000014     5.18860573e+02   # ~nu_mu_L          [1L]
   1000015     3.00151415e+03   # ~tau_1            [irrelevant]
   1000016     3.00880751e+03   # ~nu_tau_L         [irrelevant]
   2000001     7.04938393e+03   # ~d_R              [irrelevant]
   2000002     7.05136356e+03   # ~u_R              [irrelevant]
   2000003     7.04943468e+03   # ~s_R              [irrelevant]
   2000004     7.05136702e+03   # ~c_R              [irrelevant]
   2000005     7.16151312e+03   # ~b_2              [irrelevant]
   2000006     7.19205727e+03   # ~t_2              [irrelevant]
   2000011     5.05054724e+02   # ~e_R              [irrelevant]
   2000013     5.05095249e+02   # ~mu_R             [1L]
   2000015     3.01426808e+03   # ~tau_2            [irrelevant]
BLOCK NMIX
    1   1   9.999854734645232e-01   # neutralino mixing matrix (1,1)
    1   2  -6.016557041166651e-05   # neutralino mixing matrix (1,2)
    1   3   5.116898001169601e-03   # neutralino mixing matrix (1,3)
    1   4  -1.693102147969638e-03   # neutralino mixing matrix (1,4)
    2   1   1.298196396152619e-04   # neutralino mixing matrix (2,1)
    2   2   9.999130251329618e-01   # neutralino mixing matrix (2,2)
    2   3  -1.144509441463206e-02   # neutralino mixing matrix (2,3)
    2   4   6.552490377711157e-03   # neutralino mixing matrix (2,4)
    3   1   2.420724208822000e-03   # neutralino mixing matrix (3,1)
    3   2  -3.459991614555609e-03   # neutralino mixing matrix (3,2)
    3   3  -7.070901296490089e-01   # neutralino mixing matrix (3,3)
    3   4  -7.071108237789040e-01   # neutralino mixing matrix (3,4)
    4   1   4.814156322730807e-03   # neutralino mixing matrix (4,1)
    4   2  -1.272662593070186e-02   # neutralino mixing matrix (4,2)
    4   3  -7.070122882394858e-01   # neutralino mixing matrix (4,3)
    4   4   7.070703509338281e-01   # neutralino mixing matrix (4,4)
Block HMIX Q= 1.00000000e+03
     1     4.89499929e+02      # mu(Q) MSSM DRbar   [initial guess]
     2     3.93371545e+01      # tan(beta)(Q) MSSM DRbar Feynman gauge [1L]
     4     2.5e+05             # mA^2(Q) MSSM DRbar [irrelevant]
Block MSOFT Q= 1.00000000e+03  # MSSM DRbar SUSY breaking parameters
     1     2.00000000e+02      # M_1(Q)             [initial guess]
     2     4.00000000e+02      # M_2(Q)             [initial guess]
     3     2.00000000e+03      # M_3(Q)             [2L]
    31     5.00000000e+02      # meL(Q)             [2L]
    32     5.00000000e+02      # mmuL(Q)            [irrelevant]
    33     3.00000000e+03      # mtauL(Q)           [2L]
    34     4.99999999e+02      # meR(Q)             [2L]
    35     4.99999999e+02      # mmuR(Q)            [initial guess]
    36     3.00000000e+03      # mtauR(Q)           [2L]
    41     7.00000000e+03      # mqL1(Q)            [2L]
    42     7.00000000e+03      # mqL2(Q)            [2L]
    43     6.99999999e+03      # mqL3(Q)            [2L]
    44     7.00000000e+03      # muR(Q)             [2L]
    45     7.00000000e+03      # mcR(Q)             [2L]
    46     6.99999999e+03      # mtR(Q)             [2L]
    47     7.00000000e+03      # mdR(Q)             [2L]
    48     7.00000000e+03      # msR(Q)             [2L]
    49     7.00000000e+03      # mbR(Q)             [2L]
Block AU Q= 1.00000000e+03
  3  3     1.57871614e-05      # At(Q)              [2L]
Block AD Q= 1.00000000e+03
  3  3     8.99561673e-06      # Ab(Q)              [2L]
Block AE Q= 1.00000000e+03
  2  2     2.84230475e-06      # Amu(Q) MSSM DR-bar [1L]
  3  3     3.02719242e-06      # Atau(Q)            [2L]
)";

   std::istringstream stream(slha_input);

   gm2calc::GM2_slha_io slha;
   slha.read_from_stream(stream);
   slha.fill_slha(model);

   const double eps = std::numeric_limits<double>::epsilon();
   const double Pi = 3.14159265358979323846;

   // Block GM2CalcInput
   CHECK_CLOSE(model.get_EL(), std::sqrt(4*Pi*0.00775531), eps);
   CHECK_CLOSE(model.get_EL0(), std::sqrt(4*Pi*0.00729735), eps);

   // Block SMINPUTS
   CHECK_CLOSE(model.get_g3(), std::sqrt(4*Pi*0.1185), eps);
   CHECK_CLOSE(model.get_MZ(), 91.1876, eps);
   CHECK_CLOSE(model.get_MBMB(), 4.18, eps);
   CHECK_CLOSE(model.get_MT(), 173.34, eps);
   CHECK_CLOSE(model.get_ML(), 1.777, eps);
   CHECK_CLOSE(model.get_MM(), 0.1056583715, eps);
   CHECK_CLOSE(model.get_MD(), 4.76052706E-03, eps);
   CHECK_CLOSE(model.get_MU(), 2.40534062E-03, eps);
   CHECK_CLOSE(model.get_MS(), 1.04230487E-01, eps);
   CHECK_CLOSE(model.get_MC(), 1.27183378E+00, eps);

   // Block MASS
   CHECK_CLOSE(model.get_MW()                , 8.03773317e+01, eps);
   CHECK_CLOSE(model.get_physical().Mhh(0)   , 1.25712136e+02, eps);
   CHECK_CLOSE(model.get_physical().Mhh(1)   , 1.50002058e+03, eps);
   CHECK_CLOSE(model.get_physical().MAh(1)   , 1.50000000e+03, eps);
   CHECK_CLOSE(model.get_physical().MHpm(1)  , 1.50241108e+03, eps);
   CHECK_CLOSE(model.get_physical().MGlu     , 2.27707784e+03, eps);
   CHECK_CLOSE(model.get_physical().MChi(0)  , 2.01611468e+02, eps);
   CHECK_CLOSE(model.get_physical().MChi(1)  , 4.10040273e+02, eps);
   CHECK_CLOSE(model.get_physical().MCha(0)  , 4.09989890e+02, eps);
   CHECK_CLOSE(model.get_physical().MChi(2)  , 5.16529941e+02, eps);
   CHECK_CLOSE(model.get_physical().MChi(3)  , 5.45628749e+02, eps); // positive
   CHECK_CLOSE(model.get_physical().MCha(1)  , 5.46057190e+02, eps);
   CHECK_CLOSE(model.get_physical().MSd(0)   , 7.06303219e+03, eps);
   CHECK_CLOSE(model.get_physical().MSu(0)   , 7.06271372e+03, eps);
   CHECK_CLOSE(model.get_physical().MSs(0)   , 7.06305914e+03, eps);
   CHECK_CLOSE(model.get_physical().MSc(0)   , 7.06274067e+03, eps);
   CHECK_CLOSE(model.get_physical().MSb(0)   , 7.10826650e+03, eps);
   CHECK_CLOSE(model.get_physical().MSt(0)   , 7.16237728e+03, eps);
   CHECK_CLOSE(model.get_physical().MSe(0)   , 5.25167746e+02, eps);
   CHECK_CLOSE(model.get_physical().MSveL    , 5.18841083e+02, eps);
   CHECK_CLOSE(model.get_physical().MSm(0)   , 5.25187016e+02, eps);
   CHECK_CLOSE(model.get_physical().MSvmL    , 5.18860573e+02, eps);
   CHECK_CLOSE(model.get_physical().MStau(0) , 3.00151415e+03, eps);
   CHECK_CLOSE(model.get_physical().MSvtL    , 3.00880751e+03, eps);
   CHECK_CLOSE(model.get_physical().MSd(1)   , 7.04938393e+03, eps);
   CHECK_CLOSE(model.get_physical().MSu(1)   , 7.05136356e+03, eps);
   CHECK_CLOSE(model.get_physical().MSs(1)   , 7.04943468e+03, eps);
   CHECK_CLOSE(model.get_physical().MSc(1)   , 7.05136702e+03, eps);
   CHECK_CLOSE(model.get_physical().MSb(1)   , 7.16151312e+03, eps);
   CHECK_CLOSE(model.get_physical().MSt(1)   , 7.19205727e+03, eps);
   CHECK_CLOSE(model.get_physical().MSe(1)   , 5.05054724e+02, eps);
   CHECK_CLOSE(model.get_physical().MSm(1)   , 5.05095249e+02, eps);
   CHECK_CLOSE(model.get_physical().MStau(1) , 3.01426808e+03, eps);

   CHECK_CLOSE(model.get_MA0(), model.get_physical().MAh(1), eps);

   // Block NMIX (scale is optional)
   for (int i = 0; i < 4; i++) {
      CHECK(model.get_physical().ZN.row(i).cwiseAbs().maxCoeff() > 0.0);
   }

   // Block HMIX
   const double tanb = model.get_TB();
   const double scb = tanb / (1 + tanb*tanb); // sin(beta)*cos(beta)

   CHECK_CLOSE(model.get_Mu()   , 4.89499929e+02, eps);
   CHECK_CLOSE(model.get_TB()   , 3.93371545e+01, eps);
   CHECK_CLOSE(model.get_BMu()  , 2.5e+05*scb   , eps);
   CHECK_CLOSE(model.get_scale(), 1.00000000e+03, eps);

   CHECK_CLOSE(model.get_MUDIM(), model.get_scale(), eps);

   // Block MSOFT
   CHECK_CLOSE(model.get_MassB() , 2.00000000e+02, eps);
   CHECK_CLOSE(model.get_MassWB(), 4.00000000e+02, eps);
   CHECK_CLOSE(model.get_MassG() , 2.00000000e+03, eps);
   CHECK_CLOSE(model.get_ml2(0,0), gm2calc::sqr(5.00000000e+02), eps);
   CHECK_CLOSE(model.get_ml2(1,1), gm2calc::sqr(5.00000000e+02), eps);
   CHECK_CLOSE(model.get_ml2(2,2), gm2calc::sqr(3.00000000e+03), eps);
   CHECK_CLOSE(model.get_me2(0,0), gm2calc::sqr(4.99999999e+02), eps);
   CHECK_CLOSE(model.get_me2(1,1), gm2calc::sqr(4.99999999e+02), eps);
   CHECK_CLOSE(model.get_me2(2,2), gm2calc::sqr(3.00000000e+03), eps);
   CHECK_CLOSE(model.get_mq2(0,0), gm2calc::sqr(7.00000000e+03), eps);
   CHECK_CLOSE(model.get_mq2(1,1), gm2calc::sqr(7.00000000e+03), eps);
   CHECK_CLOSE(model.get_mq2(2,2), gm2calc::sqr(6.99999999e+03), eps);
   CHECK_CLOSE(model.get_mu2(0,0), gm2calc::sqr(7.00000000e+03), eps);
   CHECK_CLOSE(model.get_mu2(1,1), gm2calc::sqr(7.00000000e+03), eps);
   CHECK_CLOSE(model.get_mu2(2,2), gm2calc::sqr(6.99999999e+03), eps);
   CHECK_CLOSE(model.get_md2(0,0), gm2calc::sqr(7.00000000e+03), eps);
   CHECK_CLOSE(model.get_md2(1,1), gm2calc::sqr(7.00000000e+03), eps);
   CHECK_CLOSE(model.get_md2(2,2), gm2calc::sqr(7.00000000e+03), eps);

   // Block AU
   CHECK_CLOSE(model.get_Au(0,0), 0             , eps);
   CHECK_CLOSE(model.get_Au(1,1), 0             , eps);
   CHECK_CLOSE(model.get_Au(2,2), 1.57871614e-05, eps);

   // Block AD
   CHECK_CLOSE(model.get_Ad(0,0), 0             , eps);
   CHECK_CLOSE(model.get_Ad(1,1), 0             , eps);
   CHECK_CLOSE(model.get_Ad(2,2), 8.99561673e-06, eps);

   // Block AE
   CHECK_CLOSE(model.get_Ae(0,0), 0             , eps);
   CHECK_CLOSE(model.get_Ae(1,1), 2.84230475e-06, eps);
   CHECK_CLOSE(model.get_Ae(2,2), 3.02719242e-06, eps);

   const auto TYu = (model.get_Au() * model.get_Yu()).eval();
   const auto TYd = (model.get_Ad() * model.get_Yd()).eval();
   const auto TYe = (model.get_Ae() * model.get_Ye()).eval();

   for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
         CHECK_CLOSE(model.get_TYu(i,k), TYu(i,k), eps);
         CHECK_CLOSE(model.get_TYd(i,k), TYd(i,k), eps);
         CHECK_CLOSE(model.get_TYe(i,k), TYe(i,k), eps);
      }
   }
}


// tests that EL and EL0 have non-zero default values
TEST_CASE("fill_slha_defaults")
{
   char const * const slha_input = "Block HMIX Q= 1.00000000e+03";

   std::istringstream stream(slha_input);
   gm2calc::GM2_slha_io slha;
   slha.read_from_stream(stream);

   gm2calc::MSSMNoFV_onshell model;
   slha.fill_slha(model);

   // Block GM2CalcInput
   CHECK(model.get_EL() > 0.0);
   CHECK(model.get_EL0() > 0.0);
}


// Tests SLHA input with mutliple blocks at different scales.
// The last blocks should be preferred.
TEST_CASE("fill_slha_multiple_blocks")
{
   gm2calc::MSSMNoFV_onshell model;

   char const * const slha_input = R"(
# should be ignored:
Block HMIX Q= 500
     1     1                   # mu
Block MSOFT Q= 500
     1     4                   # M_1
Block AU Q= 500
  3  3     5                   # At
Block AD Q= 500
  3  3     6                   # Ab
Block AE Q= 500
  3  3     7                   # Atau

# should be used:
Block HMIX Q= 1000
     1     10                  # mu
Block MSOFT Q= 1000
     1     40                  # M_1
Block AU Q= 1000
  3  3     50                  # At
Block AD Q= 1000
  3  3     60                  # Ab
Block AE Q= 1000
  3  3     70                  # Atau

# should be ignored:
Block MSOFT Q= 2000
     1     400                 # M_1
Block AU Q= 2000
  3  3     500                 # At
Block AD Q= 2000
  3  3     600                 # Ab
Block AE Q= 2000
  3  3     700                 # Atau
)";

   std::istringstream stream(slha_input);
   gm2calc::GM2_slha_io slha;
   slha.read_from_stream(stream);
   slha.fill_slha(model);

   const double eps = std::numeric_limits<double>::epsilon();

   CHECK_CLOSE(model.get_scale(), 1000, eps);
   CHECK_CLOSE(model.get_Mu()   , 10  , eps);
   CHECK_CLOSE(model.get_MassB(), 40  , eps);
   CHECK_CLOSE(model.get_Au(2,2), 50  , eps);
   CHECK_CLOSE(model.get_Ad(2,2), 60  , eps);
   CHECK_CLOSE(model.get_Ae(2,2), 70  , eps);
}
