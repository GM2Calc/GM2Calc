#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "gm2calc/GeneralTHDM.hpp"
#include "gm2calc/gm2_1loop.hpp"
#include "gm2calc/gm2_2loop.hpp"
#include "GeneralTHDM/gm2_2loop_helpers.hpp"
#include "read_data.hpp"

#include <cmath>

#define CHECK_CLOSE(a,b,eps)                            \
   do {                                                 \
      CHECK((a) == doctest::Approx(b).epsilon(eps));    \
   } while (0)


namespace {

double sqr(double x) noexcept { return x*x; }

const double pi = 3.1415926535897932;

} // anonymous namespace


TEST_CASE("tree-level-spectrum")
{
   const double eps = 1e-14;

   // parameter point where choice of range
   // -pi/2 <= beta - alpha_h <= pi/2
   // matters
   gm2calc::GeneralTHDM::General_basis basis;
   basis.lambda1 = 0.26249;
   basis.lambda2 = 0.23993;
   basis.lambda3 = 2.09923;
   basis.lambda4 = -1.27781;
   basis.lambda5 = -0.71038;
   basis.lambda6 = 0.0;
   basis.lambda7 = 0.0;
   basis.tan_beta = 3.0;
   basis.M122 = sqr(200.0);

   gm2calc::GeneralTHDM model(basis);

   CHECK(!model.get_problems().have_problem());
   CHECK(model.get_MVG() == 0.0);
   CHECK(model.get_MVP() == 0.0);
   CHECK_CLOSE(model.get_MVZ(), model.get_sm().get_mz(), eps);
   CHECK_CLOSE(model.get_MVWm(), model.get_sm().get_mw(), eps);

   const double m122 = model.get_M122();
   const double tb = model.get_tan_beta();
   const double ctb = 1./tb;
   const double v_sqr = model.get_v_sqr();
   const double sb = tb/std::sqrt(1.0 + sqr(tb));
   const double cb = 1./std::sqrt(1.0 + sqr(tb));
   const double s2b = 2*sb*cb;
   const double sb2 = sqr(sb);
   const double cb2 = sqr(cb);
   const double s3b = 3*sb - 4*sb*sb2;
   const double c3b = 4*cb*cb2 - 3*cb;
   const double c2b = cb2 - sb2;
   const double l1 = model.get_Lambda1();
   const double l2 = model.get_Lambda2();
   const double l3 = model.get_Lambda3();
   const double l4 = model.get_Lambda4();
   const double l5 = model.get_Lambda5();
   const double l6 = model.get_Lambda6();
   const double l7 = model.get_Lambda7();

   // CP-odd Higgs boson
   const double mA2 = m122/sb/cb - 0.5*v_sqr*(2*l5 + l6*ctb + l7*tb);

   CHECK_CLOSE(model.get_MAh(0), model.get_MVZ(), eps);
   CHECK_CLOSE(model.get_MAh(1), std::sqrt(mA2), eps);

   // charged Higgs boson
   const double mHp2 = mA2 + 0.5*v_sqr*(l5 - l4);

   CHECK_CLOSE(model.get_MHm(0), model.get_MVWm(), eps);
   CHECK_CLOSE(model.get_MHm(1), std::sqrt(mHp2), eps);

   // CP-even Higgs bosons
   const double M112 =  mA2*sb2 + v_sqr*(l1*cb2 + 2*l6*sb*cb + l5*sb2);
   const double M122 = -mA2*sb*cb + v_sqr*((l3 + l4)*sb*cb+l6*cb2 + l7*sb2);
   const double M222 =  mA2*cb2 + v_sqr*(l2*sb2 + 2*l7*sb*cb + l5*cb2);
   const double mh2 = 0.5*(M112 + M222 - std::sqrt(sqr(M112 - M222) + 4*sqr(M122)));
   const double mH2 = 0.5*(M112 + M222 + std::sqrt(sqr(M112 - M222) + 4*sqr(M122)));

   CHECK_CLOSE(model.get_Mhh(0), std::sqrt(mh2), eps);
   CHECK_CLOSE(model.get_Mhh(1), std::sqrt(mH2), eps);

   // CP-even Higgs mixing angle alpha_h
   const double l345 = l3 + l4 + l5;
   const double lhat = 0.5*s2b*(l1*cb2 - l2*sb2 - l345*c2b) - l6*cb*c3b - l7*sb*s3b;
   const double lA = c2b*(l1*cb2 - l2*sb2) + l345*s2b*s2b - l5 + 2*l6*cb*s3b - 2*l7*sb*c3b;
   const double s2ba = 2.*lhat*v_sqr;
   const double c2ba = -(mA2 - lA*v_sqr);
   const double bma = 0.5*std::atan2(s2ba, c2ba);
   const double alpha_h = model.get_beta() - bma;

   CHECK_CLOSE(model.get_alpha_h(), alpha_h, eps);
   CHECK_CLOSE(model.get_eta(), pi/2 - bma, eps);

   // fermions
   CHECK_CLOSE(model.get_MFu(0), model.get_sm().get_mu(0), eps);
   CHECK_CLOSE(model.get_MFu(1), model.get_sm().get_mu(1), eps);
   CHECK_CLOSE(model.get_MFu(2), model.get_sm().get_mu(2), eps);
   CHECK_CLOSE(model.get_MFd(0), model.get_sm().get_md(0), eps);
   CHECK_CLOSE(model.get_MFd(1), model.get_sm().get_md(1), eps);
   CHECK_CLOSE(model.get_MFd(2), model.get_sm().get_md(2), eps);
   CHECK_CLOSE(model.get_MFe(0), model.get_sm().get_ml(0), eps);
   CHECK_CLOSE(model.get_MFe(1), model.get_sm().get_ml(1), eps);
   CHECK_CLOSE(model.get_MFe(2), model.get_sm().get_ml(2), eps);
   CHECK_CLOSE(model.get_MFv(0), 0.0, eps);
   CHECK_CLOSE(model.get_MFv(1), 0.0, eps);
   CHECK_CLOSE(model.get_MFv(2), 0.0, eps);
}


// This test ensures that the Goldstone bosons are always at the first
// position (index 0) in the MAh and MHm multiplets, even if they are
// heavier than the corresponding Higgs bosons.
TEST_CASE("test_light_Higgs_spectrum")
{
   const double eps = 1e-14;
   const double tb = 2.0;
   const double v = 245.0;

   gm2calc::GeneralTHDM_mass_eigenstates model;
   model.set_tan_beta_and_v(tb, v);
   model.set_alpha_em_and_cw(1.0/137.0, 80.0/91.0);
   model.set_Lambda1(0.0);
   model.set_Lambda2(0.0);
   model.set_Lambda3(0.0);
   model.set_Lambda4(0.0);
   model.set_Lambda5(0.0);
   model.set_Lambda6(0.0);
   model.set_Lambda7(0.0);
   model.set_M122(0.1);
   model.calculate_MSbar_masses();

   CHECK(!model.get_problems().have_problem());

   // Higgs bosons are lighter than W and Z
   CHECK_LT(model.get_MAh(1), model.get_MVZ());
   CHECK_LT(model.get_MHm(1), model.get_MVWm());
   // Goldstone bosons are at position 0
   CHECK_CLOSE(model.get_MAh(0), model.get_MVZ(), eps);
   CHECK_CLOSE(model.get_MHm(0), model.get_MVWm(), eps);
}


TEST_CASE("general_basis")
{
   const double eps = 1e-14;

   gm2calc::GeneralTHDM::General_basis basis;
   basis.lambda1 = 0.7;
   basis.lambda2 = 0.6;
   basis.lambda3 = 0.5;
   basis.lambda4 = 0.4;
   basis.lambda5 = 0.3;
   basis.lambda6 = 0.2;
   basis.lambda7 = 0.1;
   basis.tan_beta = 20;
   basis.M122 = sqr(200);

   // initialize using the basis
   gm2calc::GeneralTHDM model2(basis);

   // initialize by hand
   gm2calc::GeneralTHDM_mass_eigenstates model1;
   model1.set_tan_beta_and_v(basis.tan_beta, model2.get_sm().get_v());
   model1.set_alpha_em_and_cw(model2.get_sm().get_alpha_em_mz(), model2.get_sm().get_mw()/model2.get_sm().get_mz());
   model1.set_Lambda1(basis.lambda1);
   model1.set_Lambda2(basis.lambda2);
   model1.set_Lambda3(basis.lambda3);
   model1.set_Lambda4(basis.lambda4);
   model1.set_Lambda5(basis.lambda5);
   model1.set_Lambda6(basis.lambda6);
   model1.set_Lambda7(basis.lambda7);
   model1.set_M122(basis.M122);
   model1.calculate_MSbar_masses();

   CHECK_CLOSE(model1.get_Mhh(0), model2.get_Mhh(0), eps);
   CHECK_CLOSE(model1.get_Mhh(1), model2.get_Mhh(1), eps);
   CHECK_CLOSE(model1.get_MAh(1), model2.get_MAh(1), eps);
   CHECK_CLOSE(model1.get_MHm(1), model2.get_MHm(1), eps);
}


TEST_CASE("physical_basis")
{
   const double eps = 1e-14;

   gm2calc::GeneralTHDM::Physical_basis basis;
   basis.mh = 125;
   basis.mH = 400;
   basis.mA = 420;
   basis.mHp = 440;
   basis.sin_beta_minus_alpha = 0.999;
   basis.lambda6 = 0.1;
   basis.lambda7 = 0.2;
   basis.tan_beta = 3;
   basis.M122 = 4000;

   // initialize using the basis
   gm2calc::GeneralTHDM model2(basis);
   CHECK(!model2.get_problems().have_problem());

   // initialize by hand
   gm2calc::GeneralTHDM_mass_eigenstates model1;
   model1.set_g1(model2.get_g1());
   model1.set_g2(model2.get_g2());
   model1.set_v1(model2.get_v1());
   model1.set_v2(model2.get_v2());
   model1.set_Lambda1(model2.get_Lambda1());
   model1.set_Lambda2(model2.get_Lambda2());
   model1.set_Lambda3(model2.get_Lambda3());
   model1.set_Lambda4(model2.get_Lambda4());
   model1.set_Lambda5(model2.get_Lambda5());
   model1.set_Lambda6(model2.get_Lambda6());
   model1.set_Lambda7(model2.get_Lambda7());
   model1.set_M122(basis.M122);
   // recalculate mass spectrum from Lagrangian parameters
   model1.calculate_MSbar_masses();
   CHECK(!model1.get_problems().have_problem());

   CHECK_CLOSE(model1.get_Mhh(0), basis.mh, eps);
   CHECK_CLOSE(model1.get_Mhh(1), basis.mH, eps);
   CHECK_CLOSE(model1.get_MAh(1), basis.mA, eps);
   CHECK_CLOSE(model1.get_MHm(1), basis.mHp, eps);

   CHECK_CLOSE(model1.get_Mhh(0), model2.get_Mhh(0), eps);
   CHECK_CLOSE(model1.get_Mhh(1), model2.get_Mhh(1), eps);
   CHECK_CLOSE(model1.get_MAh(1), model2.get_MAh(1), eps);
   CHECK_CLOSE(model1.get_MHm(1), model2.get_MHm(1), eps);
}


TEST_CASE("2HDMC-demo-point")
{
   gm2calc::GeneralTHDM::General_basis basis;
   basis.lambda1 = 4.81665;
   basis.lambda2 = 0.23993;
   basis.lambda3 = 2.09923;
   basis.lambda4 = -1.27781;
   basis.lambda5 = -0.71038;
   basis.lambda6 = 0.0;
   basis.lambda7 = 0.0;
   basis.tan_beta = 3.0;
   basis.M122 = sqr(200.0);

   gm2calc::GeneralTHDM model(basis);
   CHECK(!model.get_problems().have_problem());

   const auto amu1L = gm2calc::calculate_amu_1loop(model);
   const auto amu2L = gm2calc::calculate_amu_2loop(model);

   // Notes on the 2HDMC result:
   // * the 1- and 2-loop SM Higgs contributions are not subtracted
   // * at 2-loop only the ferimonic Barr-Zee contributions from
   //   neutral Higgs bosons are implemented

   const auto amu1LSM = 2.08436e-14;

   CHECK_CLOSE((amu1L + amu1LSM)*1e14, 1.95524, 0.05);
   // CHECK_CLOSE(amu2L*1e14, -6.80278, 0.05);
}


TEST_CASE("test-point-GAMBIT")
{
   gm2calc::GeneralTHDM::General_basis basis;
   basis.lambda1 =  2.02924518279587396;
   basis.lambda2 =  0.25812066515822629;
   basis.lambda3 =  0.81575007334344507;
   basis.lambda4 =  0.43433870128700558;
   basis.lambda5 = -0.55866546170766029;
   basis.lambda6 =  0.0;
   basis.lambda7 =  0.0;
   basis.tan_beta = 20.0;
   basis.M122 = 1428;
   basis.Xe(1,1) = 0.1;

   gm2calc::SM sm;
   sm.set_alpha_em_mz(1.0/132.23323);

   gm2calc::GeneralTHDM model(basis, sm);

   CHECK(!model.get_problems().have_problem());

   const auto amu1L = gm2calc::calculate_amu_1loop(model);
   const auto amu2L = gm2calc::calculate_amu_2loop_fermionic(model);

   CHECK_CLOSE(amu1L*1e8, 6.9952544, 1e-7);
   // CHECK_CLOSE(amu2L*1e8, -5.6596159, 1e-7);
}


// This test compares the 2-loop fermionic contributions to 2HDMC.
// 2HDMC includes the following contributions: 2-loop fermionic
// Barr-Zee diagrams with only a photon and neutral Higgs bosons
// connecting to the muon line.  I.e. the following contributions are
// missing in 2HDMC:
//
// * 2-loop bosonic contributions
// * 2-loop fermionic contributions with a Z connecting to the muon line
// * 2-loop fermionic contributions with a charged Higgs connecting to the muon line
//
// Furhermore, 2HDMC does not subtract the SM Higgs contribution.
//
// To generate the 2HDMC, the running of the Yukawa couplings has been
// disabled in the function THBM::get_coupling_huu() in 2HDMC.
TEST_CASE("2HDMC-mA-scan")
{
   const auto data = gm2calc::test::read_from_file<double>(
      std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "2HDMC" +
      PATH_SEPARATOR + "2HDMC-2L-scan-mA.txt");

   gm2calc::SM sm;
   sm.set_alpha_em_mz(1.0/127.934);
   sm.set_mu(2, 172.5);
   sm.set_mu(1, 1.42);
   sm.set_md(2, 4.75);
   sm.set_ml(2, 1.77684);

   const double amu2LSM = -1.6644614586036608e-11;
   const double amu2LCharged = 8.1484174579402084e-12;

   for (const auto& p: data) {
      const auto mA = p.at(0);

      gm2calc::GeneralTHDM::Physical_basis basis;
      basis.mh = 125;
      basis.mH = 400;
      basis.mA = mA;
      basis.mHp = 440;
      basis.sin_beta_minus_alpha = 0.999;
      basis.lambda6 = 0;
      basis.lambda7 = 0;
      basis.tan_beta = 3;
      basis.M122 = 40000;

      gm2calc::GeneralTHDM model(basis, sm);

      INFO("mA = " << mA);

      const auto amu2LF = gm2calc::calculate_amu_2loop_fermionic(model);
      const auto amu2L2HDMC = p.at(1);

      if (mA > 200 && mA < 450) {
         CHECK_CLOSE(amu2LF*1e12, (amu2L2HDMC + amu2LCharged - amu2LSM)*1e12, 0.05);
      }
   }
}
