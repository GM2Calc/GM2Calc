#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "gm2calc/GeneralTHDM.hpp"
#include "gm2calc/gm2_1loop.hpp"
#include "gm2calc/gm2_2loop.hpp"

#include <cmath>

#define CHECK_CLOSE(a,b,eps)                            \
   do {                                                 \
      CHECK((a) == doctest::Approx(b).epsilon(eps));    \
   } while (0)


namespace {

double sqr(double x) noexcept { return x*x; }

const double pi = 3.1415926535897932;
const double alpha_em = 1./127.95000; // 1./137.0359997;
const double e = std::sqrt(alpha_em*4*pi);
const double alpha_s = 0.1184;
const double g3 = std::sqrt(alpha_s*4*pi);
const double mhSM = 125.09;
const double mw = 80.385;
const double mz = 91.1876;
const double cw = mw/mz;

const Eigen::Matrix<double,3,3> mu{
   (Eigen::Matrix<double,3,3>()
    << 2.2e-3, 0.0, 0.0,
    0.0, 1.28, 0.0,
    0.0, 0.0, 173.34).finished()
};

const Eigen::Matrix<double,3,3> md{
   (Eigen::Matrix<double,3,3>()
    << 4.7e-3, 0.0, 0.0,
    0.0, 0.096, 0.0,
    0.0, 0.0, 4.18).finished()
};

const Eigen::Matrix<double,3,3> mv{Eigen::Matrix<double,3,3>::Zero()};

const Eigen::Matrix<double,3,3> ml{
   (Eigen::Matrix<double,3,3>()
    << 510.998946e-6, 0.0, 0.0,
    0.0, 0.105658375, 0.0,
    0.0, 0.0, 1.77686).finished()
};

struct THDM_pars {
   double lambda1{0.0};
   double lambda2{0.0};
   double lambda3{0.0};
   double lambda4{0.0};
   double lambda5{0.0};
   double lambda6{0.0};
   double lambda7{0.0};
   double tan_beta{0.0};
   double M122{0.0};
   Eigen::Matrix<double,3,3> Xu{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> Xd{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> Xe{Eigen::Matrix<double,3,3>::Zero()};
   gm2calc::SM sm;
};

gm2calc::GeneralTHDM setup(const THDM_pars& pars)
{
   // parameter point from 2HDMC Demo.cpp
   gm2calc::GeneralTHDM model;
   model.set_sm(pars.sm);
   model.set_tan_beta(pars.tan_beta);
   model.set_Lambda1(pars.lambda1);
   model.set_Lambda2(pars.lambda2);
   model.set_Lambda3(pars.lambda3);
   model.set_Lambda4(pars.lambda4);
   model.set_Lambda5(pars.lambda5);
   model.set_Lambda6(pars.lambda6);
   model.set_Lambda7(pars.lambda7);
   model.set_M122(pars.M122);
   model.set_Xu(pars.Xu);
   model.set_Xd(pars.Xd);
   model.set_Xe(pars.Xe);

   const double vd = model.get_v1();
   const double vu = model.get_v2();

   model.set_Yu(std::sqrt(2.0)*mu/vu - vd/vu*pars.Xu);
   model.set_Yd(std::sqrt(2.0)*md/vd - vu/vd*pars.Xd);
   model.set_Ye(std::sqrt(2.0)*ml/vd - vu/vd*pars.Xe);

   model.solve_ewsb();
   model.calculate_MSbar_masses();

   model.get_sm().set_mh(mhSM);

   return model;
}

} // anonymous namespace


TEST_CASE("tree-level-spectrum")
{
   const double eps = 1e-14;

   // parameter point where choice of range
   // -pi/2 <= beta - alpha_h <= pi/2
   // matters
   THDM_pars pars;
   pars.lambda1 = 0.26249;
   pars.lambda2 = 0.23993;
   pars.lambda3 = 2.09923;
   pars.lambda4 = -1.27781;
   pars.lambda5 = -0.71038;
   pars.lambda6 = 0.0;
   pars.lambda7 = 0.0;
   pars.tan_beta = 3.0;
   pars.M122 = sqr(200.0);

   auto model = setup(pars);
   const bool have_problem = model.get_problems().have_problem();

   CHECK(!have_problem);
   CHECK(model.get_MVG() == 0.0);
   CHECK(model.get_MVP() == 0.0);
   CHECK_CLOSE(model.get_MVZ(), mz, eps);
   CHECK_CLOSE(model.get_MVWm(), mw, eps);

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
   CHECK_CLOSE(model.get_MFu(0), mu(0,0), eps);
   CHECK_CLOSE(model.get_MFu(1), mu(1,1), eps);
   CHECK_CLOSE(model.get_MFu(2), mu(2,2), eps);
   CHECK_CLOSE(model.get_MFd(0), md(0,0), eps);
   CHECK_CLOSE(model.get_MFd(1), md(1,1), eps);
   CHECK_CLOSE(model.get_MFd(2), md(2,2), eps);
   CHECK_CLOSE(model.get_MFv(0), mv(0,0), eps);
   CHECK_CLOSE(model.get_MFv(1), mv(1,1), eps);
   CHECK_CLOSE(model.get_MFv(2), mv(2,2), eps);
   CHECK_CLOSE(model.get_MFe(0), ml(0,0), eps);
   CHECK_CLOSE(model.get_MFe(1), ml(1,1), eps);
   CHECK_CLOSE(model.get_MFe(2), ml(2,2), eps);
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

   // initialize by hand
   gm2calc::GeneralTHDM model1;
   model1.set_alpha_em_and_cw(alpha_em, cw);
   model1.set_tan_beta_and_v(basis.tan_beta, model1.get_v_from_mW(mw));
   model1.set_Lambda1(basis.lambda1);
   model1.set_Lambda2(basis.lambda2);
   model1.set_Lambda3(basis.lambda3);
   model1.set_Lambda4(basis.lambda4);
   model1.set_Lambda5(basis.lambda5);
   model1.set_Lambda6(basis.lambda6);
   model1.set_Lambda7(basis.lambda7);
   model1.set_M122(basis.M122);
   model1.solve_ewsb();
   model1.calculate_MSbar_masses();

   // initialize using set_basis
   gm2calc::GeneralTHDM model2;
   model2.set_alpha_em_and_cw(alpha_em, cw);
   model2.get_sm().set_mw(mw);
   CHECK_NOTHROW(model2.set_basis(basis));

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

   // initialize using set_basis
   gm2calc::GeneralTHDM model2;
   model2.set_alpha_em_and_cw(alpha_em, cw);
   model2.get_sm().set_mw(mw);
   CHECK_NOTHROW(model2.set_basis(basis));
   CHECK(!model2.get_problems().have_problem());

   // initialize by hand
   gm2calc::GeneralTHDM model1(model2);
   // recalculate mass spectrum from Lagrangian parameters
   model1.solve_ewsb();
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
   THDM_pars pars;
   pars.lambda1 = 4.81665;
   pars.lambda2 = 0.23993;
   pars.lambda3 = 2.09923;
   pars.lambda4 = -1.27781;
   pars.lambda5 = -0.71038;
   pars.lambda6 = 0.0;
   pars.lambda7 = 0.0;
   pars.tan_beta = 3.0;
   pars.M122 = sqr(200.0);

   auto model = setup(pars);
   const bool have_problem = model.get_problems().have_problem();

   CHECK(!have_problem);

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
   THDM_pars pars;
   pars.lambda1 =  2.0292452;
   pars.lambda2 =  0.25812067;
   pars.lambda3 =  0.81575007;
   pars.lambda4 =  0.4343387;
   pars.lambda5 = -0.55866546;
   pars.lambda6 =  0.0;
   pars.lambda7 =  0.0;
   pars.tan_beta = 20.0;
   pars.M122 = 1428;
   pars.Xe(1,1) = 0.1;

   auto model = setup(pars);
   const bool have_problem = model.get_problems().have_problem();

   CHECK(!have_problem);

   const auto amu = gm2calc::calculate_amu_1loop(model);

   CHECK_CLOSE(amu*1e8, 6.9952544, 0.005);
}
