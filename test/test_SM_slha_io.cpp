#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "gm2calc/SM.hpp"
#include "gm2_slha_io.hpp"
#include <sstream>

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
     8     1.82000000E-02   # mnu3(pole)
     9     8.03773317E+01   # mW(pole)
    11     0.000510998928   # melectron(pole)
    12     1.00000000E-09   # mnu1(pole)
    13     0.1056583715     # mmuon(pole)
    14     1.70000000E-04   # mnu2(pole)
    21     4.76052706E-03   # md(2 GeV)
    22     2.40534062E-03   # mu(2 GeV)
    23     1.04230487E-01   # ms(2 GeV)
    24     1.27183378E+00   # mc(2 GeV)
Block GM2CalcInput
    33     125.0            # SM Higgs boson mass
Block VCKMIN                # CKM matrix in Wolfenstein parametrization
     1     0.2257           # lambda
     2     0.814            # A
     3     0.135            # rho
     4     0.349            # eta
)";

   gm2calc::SM sm;
   std::istringstream stream(slha_input);

   gm2calc::GM2_slha_io slha;
   slha.read_from_stream(stream);
   slha.fill(sm);

   CHECK(sm.get_alpha_em_mz() == 1/1.27934000E+02);
   CHECK(sm.get_alpha_s_mz()  == 0.1184        );
   CHECK(sm.get_mz()          == 9.11876000E+01);
   CHECK(sm.get_md(2)         == 4.18000000E+00);
   CHECK(sm.get_mu(2)         == 1.73340000E+02);
   CHECK(sm.get_ml(2)         == 1.77700000E+00);
   CHECK(sm.get_mv(2)         == 1.82000000E-02);
   CHECK(sm.get_mw()          == 8.03773317E+01);
   CHECK(sm.get_ml(0)         == 0.000510998928);
   CHECK(sm.get_mv(0)         == 1.00000000E-09);
   CHECK(sm.get_ml(1)         == 0.1056583715  );
   CHECK(sm.get_mv(1)         == 1.70000000E-04);
   CHECK(sm.get_md(0)         == 4.76052706E-03);
   CHECK(sm.get_mu(0)         == 2.40534062E-03);
   CHECK(sm.get_md(1)         == 1.04230487E-01);
   CHECK(sm.get_mu(1)         == 1.27183378E+00);
   CHECK(sm.get_mh()          == 125.0         );

   gm2calc::SM sm_ckm;
   sm_ckm.set_ckm_from_wolfenstein(0.2257, 0.814, 0.135, 0.349);
   const auto ckm_expected = sm_ckm.get_ckm();
   const auto ckm = sm.get_ckm();

   for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
         CHECK(ckm(i,k) == ckm_expected(i,k));
      }
   }
}


TEST_CASE("prefer_mw_from_MASS_block")
{
   char const * const slha_input = R"(
Block SMINPUTS
     9     8.03773317E+01   # mW(pole)
Block MASS
    24     8.00000000E+01   # W boson mass
)";

   gm2calc::SM sm;
   std::istringstream stream(slha_input);

   gm2calc::GM2_slha_io slha;
   slha.read_from_stream(stream);
   slha.fill(sm);

   CHECK(sm.get_mw() == 8.00000000E+01);
}
