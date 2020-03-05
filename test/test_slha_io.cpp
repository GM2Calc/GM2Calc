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
     3     0.1184           # alpha_s(MZ) SM MS-bar [2L]
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
   CHECK_CLOSE(model.get_g3(), std::sqrt(4*Pi*0.1184), eps);
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
