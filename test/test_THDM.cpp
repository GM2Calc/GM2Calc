#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "gm2calc/THDM.hpp"
#include "gm2calc/gm2_1loop.hpp"
#include "gm2calc/gm2_2loop.hpp"
#include "THDM/gm2_2loop_helpers.hpp"
#include "read_data.hpp"

#include <cmath>
#include <limits>

#define CHECK_CLOSE(a,b,eps)                            \
   do {                                                 \
      CHECK((a) == doctest::Approx(b).epsilon(eps));    \
   } while (0)


namespace {

double sqr(double x) noexcept { return x*x; }

const double pi = 3.1415926535897932;

} // anonymous namespace


TEST_CASE("int_to_yukawa_type")
{
   CHECK_THROWS(gm2calc::thdm::int_to_cpp_yukawa_type(0));

   CHECK(gm2calc::thdm::int_to_cpp_yukawa_type(1) == gm2calc::thdm::Yukawa_type::type_1);
   CHECK(gm2calc::thdm::int_to_cpp_yukawa_type(2) == gm2calc::thdm::Yukawa_type::type_2);
   CHECK(gm2calc::thdm::int_to_cpp_yukawa_type(3) == gm2calc::thdm::Yukawa_type::type_X);
   CHECK(gm2calc::thdm::int_to_cpp_yukawa_type(4) == gm2calc::thdm::Yukawa_type::type_Y);
   CHECK(gm2calc::thdm::int_to_cpp_yukawa_type(5) == gm2calc::thdm::Yukawa_type::aligned);
   CHECK(gm2calc::thdm::int_to_cpp_yukawa_type(6) == gm2calc::thdm::Yukawa_type::general);

   CHECK_THROWS(gm2calc::thdm::int_to_cpp_yukawa_type(7));
}


void test_tree_level_spectrum(gm2calc::thdm::Yukawa_type yukawa_type)
{
   const double eps = 1e-14;

   // parameter point where choice of range
   // -pi/2 <= beta - alpha_h <= pi/2
   // matters
   gm2calc::thdm::Gauge_basis basis;
   basis.yukawa_type = yukawa_type;
   basis.lambda << 0.26249, 0.23993, 2.09923, -1.27781, -0.71038, 0.0, 0.0;
   basis.tan_beta = 3.0;
   basis.m122 = sqr(200.0);
   basis.zeta_u = 5;
   basis.zeta_d = 6;
   basis.zeta_l = 7;
   basis.Pi_u = 0.2*decltype(basis.Pi_u)::Identity();
   basis.Pi_d = 0.3*decltype(basis.Pi_d)::Identity();
   basis.Pi_l = 0.4*decltype(basis.Pi_l)::Identity();

   gm2calc::THDM model(basis);

   CHECK(model.get_MVG() == 0.0);
   CHECK(model.get_MVP() == 0.0);
   CHECK_CLOSE(model.get_MVZ(), model.get_sm().get_mz(), eps);
   CHECK_CLOSE(model.get_MVWm(), model.get_sm().get_mw(), eps);

   const double m122 = model.get_m122();
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
   const double l1 = model.get_lambda1();
   const double l2 = model.get_lambda2();
   const double l3 = model.get_lambda3();
   const double l4 = model.get_lambda4();
   const double l5 = model.get_lambda5();
   const double l6 = model.get_lambda6();
   const double l7 = model.get_lambda7();

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

   // CKM matrix
   const Eigen::Matrix<std::complex<double>,3,3> ckm_input = model.get_sm().get_ckm();
   const Eigen::Matrix<std::complex<double>,3,3> ckm_output = model.get_Vu() * model.get_Vd().adjoint();

   CHECK( (ckm_input - ckm_output).cwiseAbs().maxCoeff() < eps );
}


TEST_CASE("tree-level-spectrum")
{
   test_tree_level_spectrum(gm2calc::thdm::Yukawa_type::type_1);
   test_tree_level_spectrum(gm2calc::thdm::Yukawa_type::type_2);
   test_tree_level_spectrum(gm2calc::thdm::Yukawa_type::type_X);
   test_tree_level_spectrum(gm2calc::thdm::Yukawa_type::type_Y);
   test_tree_level_spectrum(gm2calc::thdm::Yukawa_type::aligned);
   test_tree_level_spectrum(gm2calc::thdm::Yukawa_type::general);
}


TEST_CASE("tree-level-spectrum-with-tachyons")
{
   gm2calc::thdm::Gauge_basis basis;
   basis.lambda << -0.1, -0.2, 0, 0, 0, 0, 0;
   basis.tan_beta = 20;
   basis.m122 = sqr(200);

   REQUIRE_THROWS(gm2calc::THDM(basis));
}


// This test ensures that the Goldstone bosons are always at the first
// position (index 0) in the MAh and MHm multiplets, even if they are
// heavier than the corresponding Higgs bosons.
TEST_CASE("test_light_Higgs_spectrum")
{
   const double eps = 1e-14;
   const double tb = 2.0;
   const double v = 245.0;

   gm2calc::THDM_mass_eigenstates model;
   model.set_tan_beta_and_v(tb, v);
   model.set_alpha_em_and_cw(1.0/137.0, 80.0/91.0);
   model.set_lambda1(0.0);
   model.set_lambda2(0.0);
   model.set_lambda3(0.0);
   model.set_lambda4(0.0);
   model.set_lambda5(0.0);
   model.set_lambda6(0.0);
   model.set_lambda7(0.0);
   model.set_m122(0.1);
   model.calculate_MSbar_masses();

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

   gm2calc::thdm::Gauge_basis basis;
   basis.lambda << 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1;
   basis.tan_beta = 20;
   basis.m122 = sqr(200);

   // initialize using the basis
   gm2calc::THDM model2(basis);

   // initialize by hand
   gm2calc::THDM_mass_eigenstates model1;
   model1.set_tan_beta_and_v(basis.tan_beta, model2.get_sm().get_v());
   model1.set_alpha_em_and_cw(model2.get_sm().get_alpha_em_mz(), model2.get_sm().get_mw()/model2.get_sm().get_mz());
   model1.set_lambda1(basis.lambda(0));
   model1.set_lambda2(basis.lambda(1));
   model1.set_lambda3(basis.lambda(2));
   model1.set_lambda4(basis.lambda(3));
   model1.set_lambda5(basis.lambda(4));
   model1.set_lambda6(basis.lambda(5));
   model1.set_lambda7(basis.lambda(6));
   model1.set_m122(basis.m122);
   model1.calculate_MSbar_masses();

   CHECK_CLOSE(model1.get_Mhh(0), model2.get_Mhh(0), eps);
   CHECK_CLOSE(model1.get_Mhh(1), model2.get_Mhh(1), eps);
   CHECK_CLOSE(model1.get_MAh(1), model2.get_MAh(1), eps);
   CHECK_CLOSE(model1.get_MHm(1), model2.get_MHm(1), eps);
}


TEST_CASE("physical_basis")
{
   const double eps = 1e-14;

   gm2calc::thdm::Mass_basis basis;
   basis.mh = 125;
   basis.mH = 400;
   basis.mA = 420;
   basis.mHp = 440;
   basis.sin_beta_minus_alpha = 0.999;
   basis.lambda_6 = 0.1;
   basis.lambda_7 = 0.2;
   basis.tan_beta = 3;
   basis.m122 = 4000;

   // initialize using the basis
   gm2calc::THDM model2(basis);

   // initialize by hand
   gm2calc::THDM_mass_eigenstates model1;
   model1.set_g1(model2.get_g1());
   model1.set_g2(model2.get_g2());
   model1.set_v1(model2.get_v1());
   model1.set_v2(model2.get_v2());
   model1.set_lambda1(model2.get_lambda1());
   model1.set_lambda2(model2.get_lambda2());
   model1.set_lambda3(model2.get_lambda3());
   model1.set_lambda4(model2.get_lambda4());
   model1.set_lambda5(model2.get_lambda5());
   model1.set_lambda6(model2.get_lambda6());
   model1.set_lambda7(model2.get_lambda7());
   model1.set_m122(basis.m122);
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
   gm2calc::thdm::Gauge_basis basis;
   basis.lambda << 4.81665, 0.23993, 2.09923, -1.27781, -0.71038, 0.0, 0.0;
   basis.tan_beta = 3.0;
   basis.m122 = sqr(200.0);

   gm2calc::THDM model(basis);

   const auto amu1L = gm2calc::calculate_amu_1loop(model);
   // const auto amu2L = gm2calc::calculate_amu_2loop(model);

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
   gm2calc::thdm::Gauge_basis basis;
   basis.yukawa_type = gm2calc::thdm::Yukawa_type::general;
   basis.lambda <<  2.02924518279587396,
                    0.25812066515822629,
                    0.81575007334344507,
                    0.43433870128700558,
                   -0.55866546170766029,
                    0.0,
                    0.0;
   basis.tan_beta = 20.0;
   basis.m122 = 1428;
   basis.Pi_l(1,1) = 0.1;

   gm2calc::thdm::Config config;
   config.running_couplings = false;

   gm2calc::SM sm;
   sm.set_alpha_em_mz(1.0/132.23323);
   sm.set_ckm_from_wolfenstein(0, 0, 0, 0);

   gm2calc::THDM model(basis, sm, config);

   const auto amu1L = gm2calc::calculate_amu_1loop(model);
   const auto amu2L = gm2calc::calculate_amu_2loop_fermionic(model);

   CHECK_CLOSE(amu1L*1e8, 6.9952544, 1e-7);
   CHECK_CLOSE(amu2L*1e8, 265.56157, 1e-7);
}


TEST_CASE("test-point-GAMBIT-real-CKM")
{
   gm2calc::thdm::Gauge_basis basis;
   basis.yukawa_type = gm2calc::thdm::Yukawa_type::general;
   basis.lambda <<  2.02924518279587396,
                    0.25812066515822629,
                    0.81575007334344507,
                    0.43433870128700558,
                   -0.55866546170766029,
                    0.0,
                    0.0;
   basis.tan_beta = 20.0;
   basis.m122 = 1428;
   basis.Pi_l(1,1) = 0.1;

   Eigen::Matrix<std::complex<double>,3,3> ckm;
   ckm << 0.97383946649250375, 0.22720257917454678, 0.003959989650152115,
         -0.22716464426353369, 0.97294117488820131, 0.042210124422804286,
          0.0057374121533755448, -0.04200545468865001, 0.99910090775565907;

   // check unitarity of CKM matrix
   {
      const double eps = 10*std::numeric_limits<double>::epsilon();
      const Eigen::Matrix<std::complex<double>,3,3> unit = Eigen::Matrix<std::complex<double>,3,3>::Identity();
      const Eigen::Matrix<std::complex<double>,3,3> cca = ckm * ckm.adjoint();
      const double max_diff = (cca - unit).cwiseAbs().maxCoeff();
      CHECK(max_diff <= eps);
   }

   gm2calc::thdm::Config config;
   config.running_couplings = false;

   gm2calc::SM sm;
   sm.set_alpha_em_mz(1.0/132.23323);
   sm.set_ckm(ckm);

   gm2calc::THDM model(basis, sm, config);

   const auto amu1L = gm2calc::calculate_amu_1loop(model);
   const auto amu2L = gm2calc::calculate_amu_2loop_fermionic(model);

   // effect of CKM mixing in 2-loop fermionic Barr-Zee charged Higgs contribution
   const auto ckm_mixing_charged_higgs = 0.084792223;

   CHECK_CLOSE(amu1L*1e8, 6.9952544, 1e-7);
   CHECK_CLOSE(amu2L*1e8, 265.47676 + ckm_mixing_charged_higgs, 1e-8);
}


TEST_CASE("test-point-GAMBIT-complex-CKM")
{
   gm2calc::thdm::Gauge_basis basis;
   basis.yukawa_type = gm2calc::thdm::Yukawa_type::general;
   basis.lambda <<  2.02924518279587396,
                    0.25812066515822629,
                    0.81575007334344507,
                    0.43433870128700558,
                   -0.55866546170766029,
                    0.0,
                    0.0;
   basis.tan_beta = 20.0;
   basis.m122 = 1428;
   basis.Pi_u << 0.0, 0.0, 0.0, 0.0, 0.3, 0.05, 0.0, 0.05, 0.3;
   basis.Pi_d << 0.0, 0.0, 0.0, 0.0, 0.2, 0.03, 0.0, 0.03, 0.2;
   basis.Pi_l << 0.0, 0.0, 0.0, 0.0, 0.1, 0.01, 0.0, 0.01, 0.1;

   Eigen::Matrix<std::complex<double>,3,3> ckm;
   ckm << 1.0, 0.22537, std::complex<double>(0.0010901809,-0.0032891784),
         -0.22537, 1.0, 0.041344392,
         std::complex<double>(0.0082276048,-0.0032891784), -0.041344392, 1.0;

   gm2calc::thdm::Config config;
   config.running_couplings = false;

   gm2calc::SM sm;
   sm.set_alpha_em_mz(1.0/132.23323);
   sm.set_ckm(ckm);

   gm2calc::THDM model(basis, sm, config);

   // const auto amu2L = gm2calc::calculate_amu_2loop_fermionic(model);

   // CHECK_CLOSE(amu2L*1e8, 173.724, 1e-6);
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
TEST_CASE("2HDMC-mA-scan-no-running")
{
   const auto data = gm2calc::test::read_from_file<double>(
      std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "2HDMC" +
      PATH_SEPARATOR + "2HDMC-scan-mA-no-running.txt");

   gm2calc::thdm::Config config;
   config.running_couplings = false;

   gm2calc::SM sm;
   sm.set_alpha_em_mz(1.0/127.934);
   sm.set_mu(2, 172.5);
   sm.set_mu(1, 1.42);
   sm.set_md(2, 4.75);
   sm.set_ml(2, 1.77684);

   const double amu1LSM = 2.2334814296705501e-14;
   const double amu2LSM = -1.6644614586036608e-11;
   const double amu2LCharged = 8.1484174579402084e-12;

   for (const auto& p: data) {
      const auto mA = p.at(0);

      gm2calc::thdm::Mass_basis basis;
      basis.mh = 125;
      basis.mH = 400;
      basis.mA = mA;
      basis.mHp = 440;
      basis.sin_beta_minus_alpha = 0.999;
      basis.lambda_6 = 0;
      basis.lambda_7 = 0;
      basis.tan_beta = 3;
      basis.m122 = 40000;

      gm2calc::THDM model(basis, sm, config);

      INFO("mA = " << mA);

      const auto amu1L = gm2calc::calculate_amu_1loop(model);
      const auto amu2LF = gm2calc::calculate_amu_2loop_fermionic(model);
      const auto amu1L2HDMC = p.at(1);
      const auto amu2L2HDMC = p.at(2);

      CHECK_CLOSE(amu1L*1e13, (amu1L2HDMC - amu1LSM)*1e13, 0.03);

      // For mA < 200 the Z contribution in GM2Calc becomes larger, so
      // the test must be restricted to mA > 200.
      //
      // For mA > 450 the fact that running fermion masses are used in
      // 2HDMC becomes more relevant.  So the test must be restricted
      // to mA < 450.
      if (mA > 200 && mA < 450) {
         CHECK_CLOSE(amu2LF*1e12, (amu2L2HDMC + amu2LCharged - amu2LSM)*1e12, 0.05);
      }
   }
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
TEST_CASE("2HDMC-mA-scan-with-running")
{
   const auto data = gm2calc::test::read_from_file<double>(
      std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "2HDMC" +
      PATH_SEPARATOR + "2HDMC-scan-mA-with-running.txt");

   gm2calc::thdm::Config config;
   config.running_couplings = true;

   gm2calc::SM sm;
   sm.set_alpha_em_mz(1.0/127.934);
   sm.set_mu(2, 172.5);
   sm.set_mu(1, 1.42);
   sm.set_md(2, 4.75);
   sm.set_ml(2, 1.77684);

   const double amu1LSM = 2.2334814296705501e-14;
   const double amu2LSM = -1.6644614586036608e-11;
   const double amu2LCharged = 8.1484174579402084e-12;

   for (const auto& p: data) {
      const auto mA = p.at(0);

      gm2calc::thdm::Mass_basis basis;
      basis.mh = 125;
      basis.mH = 400;
      basis.mA = mA;
      basis.mHp = 440;
      basis.sin_beta_minus_alpha = 0.999;
      basis.lambda_6 = 0;
      basis.lambda_7 = 0;
      basis.tan_beta = 3;
      basis.m122 = 40000;

      gm2calc::THDM model(basis, sm, config);

      INFO("mA = " << mA);

      const auto amuGM2Calc = gm2calc::calculate_amu_1loop(model)
                            + gm2calc::calculate_amu_2loop_fermionic(model);
      const auto amu2HDMC = p.at(1) + amu2LCharged - amu1LSM - amu2LSM;

      // For larger mA the fact that proper RG running for the fermion
      // masses is used in 2HDMC becomes more relevant.  So the test
      // is less precise for larger mA
      if (mA < 200) {
         CHECK_CLOSE(amuGM2Calc*1e13, amu2HDMC*1e13, 0.02);
      } else if (mA < 430) {
         CHECK_CLOSE(amuGM2Calc*1e13, amu2HDMC*1e13, 0.03);
      } else {
         CHECK_CLOSE(amuGM2Calc*1e13, amu2HDMC*1e13, 0.04);
      }
   }
}


// This test checks that the input parametrization with Pi_f is
// equivalent the input parametrization with Delta_f.
TEST_CASE("parametrizations")
{
   gm2calc::SM sm;
   sm.set_ckm(Eigen::Matrix<double,3,3>::Identity());

   const double v = sm.get_v();
   const Eigen::Matrix<double,3,3> mu = sm.get_mu().asDiagonal();
   const Eigen::Matrix<double,3,3> md = sm.get_md().asDiagonal();
   const Eigen::Matrix<double,3,3> ml = sm.get_ml().asDiagonal();
   const double mh = 125;
   const double mH = 400;
   const double mA = 420;
   const double mHp = 440;
   const double sin_beta_minus_alpha = 0.995;
   const double lambda_6 = 0.1;
   const double lambda_7 = 0.2;
   const double tan_beta = 3;
   const double cos_beta = 1/std::sqrt(1 + tan_beta*tan_beta);
   const double m122 = 40000;
   const double zeta_u = 10;
   const double zeta_d = 20;
   const double zeta_l = 30;
   const Eigen::Matrix<double,3,3> Delta_u((Eigen::Matrix<double,3,3>() << 1, 2, 3, 4, 5, 6, 7, 8, 9).finished());
   const Eigen::Matrix<double,3,3> Delta_d(2.0*Delta_u);
   const Eigen::Matrix<double,3,3> Delta_l(3.0*Delta_u);

   gm2calc::thdm::Config config;
   config.running_couplings = false;

   gm2calc::thdm::Mass_basis aligned_basis;
   aligned_basis.yukawa_type = gm2calc::thdm::Yukawa_type::aligned;
   aligned_basis.mh = mh;
   aligned_basis.mH = mH;
   aligned_basis.mA = mA;
   aligned_basis.mHp = mHp;
   aligned_basis.sin_beta_minus_alpha = sin_beta_minus_alpha;
   aligned_basis.lambda_6 = lambda_6;
   aligned_basis.lambda_7 = lambda_7;
   aligned_basis.tan_beta = tan_beta;
   aligned_basis.m122 = m122;
   aligned_basis.zeta_u = zeta_u;
   aligned_basis.zeta_d = zeta_d;
   aligned_basis.zeta_l = zeta_l;
   aligned_basis.Delta_u = Delta_u;
   aligned_basis.Delta_d = Delta_d;
   aligned_basis.Delta_l = Delta_l;

   gm2calc::thdm::Mass_basis general_basis;
   general_basis.yukawa_type = gm2calc::thdm::Yukawa_type::general;
   general_basis.mh = mh;
   general_basis.mH = mH;
   general_basis.mA = mA;
   general_basis.mHp = mHp;
   general_basis.sin_beta_minus_alpha = sin_beta_minus_alpha;
   general_basis.lambda_6 = lambda_6;
   general_basis.lambda_7 = lambda_7;
   general_basis.tan_beta = tan_beta;
   general_basis.m122 = m122;
   general_basis.zeta_u = zeta_u;
   general_basis.zeta_d = zeta_d;
   general_basis.zeta_l = zeta_l;
   // change parametrization: Input Pi_f instead of Delta_f
   general_basis.Pi_u = cos_beta*(std::sqrt(2.0)*mu/v*(zeta_u + tan_beta) + Delta_u);
   general_basis.Pi_d = cos_beta*(std::sqrt(2.0)*md/v*(zeta_d + tan_beta) + Delta_d);
   general_basis.Pi_l = cos_beta*(std::sqrt(2.0)*ml/v*(zeta_l + tan_beta) + Delta_l);

   const gm2calc::THDM aligned_thdm(aligned_basis, sm, config);
   const gm2calc::THDM general_thdm(general_basis, sm, config);

   const double aligned_amu = gm2calc::calculate_amu_1loop(aligned_thdm)
                            + gm2calc::calculate_amu_2loop_fermionic(aligned_thdm);
   const double general_amu = gm2calc::calculate_amu_1loop(general_thdm)
                            + gm2calc::calculate_amu_2loop_fermionic(general_thdm);

   CHECK_CLOSE(aligned_amu*1e10, general_amu*1e10, 1e-10);
}


double calc_amu_for_type(gm2calc::thdm::Yukawa_type yukawa_type, double tan_beta,
                         double zeta_u, double zeta_d, double zeta_l)
{
   gm2calc::thdm::Mass_basis basis;
   basis.yukawa_type = yukawa_type;
   basis.mh = 125;
   basis.mH = 400;
   basis.mA = 420;
   basis.mHp = 440;
   basis.sin_beta_minus_alpha = 0.995;
   basis.tan_beta = tan_beta;
   basis.m122 = 200*200;
   basis.zeta_u = zeta_u;
   basis.zeta_d = zeta_d;
   basis.zeta_l = zeta_l;

   const gm2calc::THDM thdm(basis);

   return gm2calc::calculate_amu_1loop(thdm) + gm2calc::calculate_amu_2loop(thdm);
}


// This test checks that the type I, II, X, Y models are special cases
// of the aligned THDM for specific choices of the zeta_f.
TEST_CASE("alignment_limits")
{
   const double eps = 1e-14;
   const double tan_beta = 3;
   const double cot_beta = 1/tan_beta;

   {
      const auto type_1  = calc_amu_for_type(gm2calc::thdm::Yukawa_type::type_1 , tan_beta, 0, 0, 0);
      const auto aligned = calc_amu_for_type(gm2calc::thdm::Yukawa_type::aligned, tan_beta, cot_beta, cot_beta, cot_beta);
      CHECK_CLOSE(type_1*1e10, aligned*1e10, eps);
   }

   {
      const auto type_2  = calc_amu_for_type(gm2calc::thdm::Yukawa_type::type_2 , tan_beta, 0, 0, 0);
      const auto aligned = calc_amu_for_type(gm2calc::thdm::Yukawa_type::aligned, tan_beta, cot_beta, -tan_beta, -tan_beta);
      CHECK_CLOSE(type_2*1e10, aligned*1e10, eps);
   }

   {
      const auto type_X  = calc_amu_for_type(gm2calc::thdm::Yukawa_type::type_X , tan_beta, 0, 0, 0);
      const auto aligned = calc_amu_for_type(gm2calc::thdm::Yukawa_type::aligned, tan_beta, cot_beta, cot_beta, -tan_beta);
      CHECK_CLOSE(type_X*1e10, aligned*1e10, eps);
   }

   {
      const auto type_Y  = calc_amu_for_type(gm2calc::thdm::Yukawa_type::type_Y , tan_beta, 0, 0, 0);
      const auto aligned = calc_amu_for_type(gm2calc::thdm::Yukawa_type::aligned, tan_beta, cot_beta, -tan_beta, cot_beta);
      CHECK_CLOSE(type_Y*1e10, aligned*1e10, eps);
   }
}


gm2calc::THDM calc_point(double sin_beta_minus_alpha)
{
   gm2calc::thdm::Mass_basis basis;
   basis.yukawa_type = gm2calc::thdm::Yukawa_type::general;
   basis.mh = 125.09;
   basis.mH = 3485;
   basis.mA = 3485;
   basis.mHp = 3429;
   basis.sin_beta_minus_alpha = sin_beta_minus_alpha;
   basis.tan_beta = 0.98;
   basis.m122 = 5885476;
   basis.Pi_u << 0, 0, 0,
                 0, 1.15, -0.64,
                 0, -0.64, 0.60;
   basis.Pi_d << 0, 0, 0,
                 0, 0.10, 0.004,
                 0, 0.004, 0.017;
   basis.Pi_l << 0, 0, 0,
                 0, -0.04, 0.75,
                 0, 0.75, -0.36;

   return gm2calc::THDM(basis);
}


// tests that the value for sin(beta - alpha_h) returned by the model
// matches the input
void test_sin_beta_minus_alpha(double sin_beta_minus_alpha)
{
   const auto model = calc_point(sin_beta_minus_alpha);
   CHECK_CLOSE(model.get_sin_beta_minus_alpha(), sin_beta_minus_alpha, 1e-14);
}


// test stability of alpha_h when h and H are swapped
TEST_CASE("alpha_h_swapped")
{
   test_sin_beta_minus_alpha(1.0);    // order of h and H switched: (H, h)
   test_sin_beta_minus_alpha(0.9999); // normal ordering: (h, H)
}
