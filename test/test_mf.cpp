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


// This test checks the renormalization group running of mb(MZ,MS-bar)
TEST_CASE("mb_beta_function")
{
   const double pi = 3.14159265358979323846;
   const double alpha_s_mz = 0.1184;
   const double mz = 91.1876;
   const double mt_pole = mz; // run only to mz, to not change alpha_s
   const double mb_mb = 4.18;
   const double Q1 = mz;
   const double Q2 = mz*(1 + 1e-4);

   const auto mb_1 = gm2calc::calculate_mb_SM6_MSbar(mb_mb, mt_pole, alpha_s_mz, mz, Q1);
   const auto mb_2 = gm2calc::calculate_mb_SM6_MSbar(mb_mb, mt_pole, alpha_s_mz, mz, Q2);

   const auto beta = (mb_1/mb_2 - 1)/std::log(Q1/Q2);

   // Note: The QCD beta function of mb is
   // dmb/d(log(Q)) = -2*alpha_s*mb/pi

   CHECK_CLOSE(beta, -2*alpha_s_mz/pi, 1e-4);
}


// This test checks the calculation of mt(MS-bar,Q=mt)
TEST_CASE("mt_MSbar_at_mt")
{
   const double alpha_s_mz = 0.1184;
   const double mz = 91.1876;
   const double mt_pole = 173.34;
   const double scale = mt_pole; // run only to mt_pole

   const auto mt_mt = gm2calc::calculate_mt_SM6_MSbar(mt_pole, alpha_s_mz, mz, scale);

   // obtained with FlexibleSUSY 2.6, using 1-loop running
   const auto mt_mt_expected = 1.64464650E+02;

   CHECK_CLOSE(mt_mt, mt_mt_expected, 0.01);
}


// This test checks the calculation of mt(MS-bar,Q=1000GeV)
TEST_CASE("mt_MSbar_at_1TeV")
{
   const double alpha_s_mz = 0.1184;
   const double mz = 91.1876;
   const double mt_pole = 173.34;
   const double scale = 1000;

   const auto mt = gm2calc::calculate_mt_SM6_MSbar(mt_pole, alpha_s_mz, mz, scale);

   // obtained with FlexibleSUSY 2.6, using 1-loop running
   const auto mt_expected = 1.49828567E+02;

   CHECK_CLOSE(mt, mt_expected, 0.02);
}


// This test checks the calculation of mt(MS-bar,1000GeV)
TEST_CASE("mtau_MSbar_at_1TeV")
{
   const double alpha_em_mz = 1/127.934;
   const double mtau_pole = 1.777;
   const double scale = 1000;

   const auto mtau = gm2calc::calculate_mtau_SM6_MSbar(mtau_pole, alpha_em_mz, scale);

   // obtained with FlexibleSUSY 2.6, using 1-loop running
   const auto mtau_expected = 1.75137065E+00;

   CHECK_CLOSE(mtau, mtau_expected, 0.01);
}


// This test checks for non-convergence of the iteration to determine
// Lambda_QCD due to a too small renormalization scale
TEST_CASE("non-convergence")
{
   const double alpha_s = 0.1184;
   const double scale = 1e-15;
   const double mb_mb = 4.2;

   const auto mb_mz = gm2calc::calculate_mb_SM5_DRbar(mb_mb, alpha_s, scale);

   CHECK_LT(mb_mz, 5.0);
   CHECK_GT(mb_mz, 1.0);
}
