#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"

#include "gm2calc/SM.hpp"

#include "gm2_numerics.hpp"
#include "gm2_slha_io.hpp"

#include <sstream>

#define CHECK_CLOSE(a,b,eps)                            \
   do {                                                 \
      CHECK((a) == doctest::Approx(b).epsilon(eps));    \
   } while (0)


TEST_CASE("fill_SM")
{
   char const * const slha_input = R"(
Block SMINPUTS
     1     1.27934000E+02   # alpha_em(MZ)^(-1) SM MS-bar
     3     0.1184           # alpha_s(MZ) SM MS-bar
     4     9.11876000E+01   # MZ(pole)
     5     4.18000000E+00   # mb(mb) SM MS-bar
     6     1.73340000E+02   # mtop(pole)
     7     1.77700000E+00   # mtau(pole)
     8     0.00000000E+00   # mnu3(pole)
     9     8.03773317E+01   # mW(pole)
    11     0.000510998928   # melectron(pole)
    12     0.00000000E+00   # mnu1(pole)
    13     0.1056583715     # mmuon(pole)
    14     0.00000000E+00   # mnu2(pole)
    21     4.76052706E-03   # md(2 GeV)
    22     2.40534062E-03   # mu(2 GeV)
    23     1.04230487E-01   # ms(2 GeV)
    24     1.27183378E+00   # mc(2 GeV)
)";

   gm2calc::SM sm;
   std::istringstream stream(slha_input);

   gm2calc::GM2_slha_io slha;
   slha.read_from_stream(stream);
   slha.fill(sm);

   const double eps = std::numeric_limits<double>::epsilon();

   CHECK_CLOSE(sm.get_alpha_em_mz(), 1.27934000E+02, eps);
   CHECK_CLOSE(sm.get_alpha_s_mz() , 0.1184        , eps);
   CHECK_CLOSE(sm.get_mz()         , 9.11876000E+01, eps);
   CHECK_CLOSE(sm.get_md(2)        , 4.18000000E+00, eps);
   CHECK_CLOSE(sm.get_mu(2)        , 1.73340000E+02, eps);
   CHECK_CLOSE(sm.get_ml(2)        , 1.77700000E+00, eps);
   CHECK_CLOSE(sm.get_mw()         , 8.03773317E+01, eps);
   CHECK_CLOSE(sm.get_ml(0)        , 0.000510998928, eps);
   CHECK_CLOSE(sm.get_ml(1)        , 0.1056583715  , eps);
   CHECK_CLOSE(sm.get_md(0)        , 4.76052706E-03, eps);
   CHECK_CLOSE(sm.get_mu(0)        , 2.40534062E-03, eps);
   CHECK_CLOSE(sm.get_md(1)        , 1.04230487E-01, eps);
   CHECK_CLOSE(sm.get_mu(1)        , 1.27183378E+00, eps);
}
