#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "gm2_mf.hpp"


#define CHECK_CLOSE(a,b,eps)                            \
   do {                                                 \
      CHECK((a) == doctest::Approx(b).epsilon(eps));    \
   } while (0)


// This test checks the implementation of the calculation of
// mb(MZ,DR-bar) against Fig.1 of hep-ph/0207126
TEST_CASE("mb_DRbar_at_mz")
{
   const double alpha_s = 0.1172;
   const double mz = 91.1876;
   const double mb_mb = 4.2;

   const auto mb_mz = gm2calc::calculate_mb_SM5_DRbar(mb_mb, alpha_s, mz);

   CHECK_CLOSE(mb_mz, 2.84, 0.01);
}


// This test checks the implementation of the calculation of
// mb(MZ,MS-bar) against Fig.1 of hep-ph/0207126
TEST_CASE("mb_MSbar_at_mz")
{
   const double alpha_s_mz = 0.1172;
   const double mz = 91.1876;
   const double mt_pole = mz; // run only to mz
   const double mb_mb = 4.18;
   const double scale = mz; // run only to mz

   const auto mb_mz = gm2calc::calculate_mb_SM6_MSbar(mb_mb, mt_pole, alpha_s_mz, mz, scale);

   CHECK_CLOSE(mb_mz, 2.88, 0.01);
}
