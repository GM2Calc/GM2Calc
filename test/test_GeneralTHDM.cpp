#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "gm2calc/GeneralTHDM.hpp"

#include <cmath>

#define CHECK_CLOSE(a,b,eps)                            \
   do {                                                 \
      CHECK((a) == doctest::Approx(b).epsilon(eps));    \
   } while (0)


namespace {

double sqr(double x) noexcept { return x*x; }

const double pi = 3.1415926535897932;
const double alpha = 1./137.036;
const double e = std::sqrt(alpha*4*pi);
const double alpha_s = 0.1184;
const double g3 = std::sqrt(alpha_s*4*pi);
const double mw = 80.379;
const double mz = 91.1876;
const double cw = mw/mz;
const double sw = std::sqrt(1 - sqr(cw));
const double gY = e/cw;
const double g1 = gY*std::sqrt(5.0/3.0);
const double g2 = e/sw;
const double v = 2*mw/g2;

const Eigen::Matrix<double,3,3> ml{
   (Eigen::Matrix<double,3,3>()
    << 510.999e-6, 0.0, 0.0,
    0.0, 0.10565837, 0.0,
    0.0, 0.0, 1.7768).finished()
};

gm2calc::GeneralTHDM setup()
{
   const double tb = 20.0;
   const double sb = tb/std::sqrt(1.0 + sqr(tb));
   const double cb = 1./std::sqrt(1.0 + sqr(tb));
   const double vu = v*sb;
   const double vd = v*cb;

   gm2calc::GeneralTHDM model;
   model.set_g1(g1);
   model.set_g2(g2);
   model.set_g3(g3);
   model.set_v1(vd);
   model.set_v2(vu);
   model.set_Lambda1(0.7);
   model.set_Lambda2(0.6);
   model.set_Lambda3(0.5);
   model.set_Lambda4(0.4);
   model.set_Lambda5(0.3);
   model.set_Lambda6(0.2);
   model.set_Lambda7(0.1);
   model.set_M122(sqr(100.0));
   model.set_Ye(std::sqrt(2.0)*ml/vd);

   model.solve_ewsb();
   model.calculate_MSbar_masses();

   return model;
}

} // anonymous namespace


TEST_CASE("tree-level-spectrum")
{
   const double eps = 1e-14;
   auto model = setup();
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

   // CP-even Higgs mixing angle alpha_H
   const double l345 = l3 + l4 + l5;
   const double lhat = 0.5*s2b*(l1*cb2 - l2*sb2 - l345*c2b) - l6*cb*c3b - l7*sb*s3b;
   const double lA = c2b*(l1*cb2 - l2*sb2) + l345*s2b*s2b - l5 + 2*l6*cb*s3b - 2*l7*sb*c3b;
   const double s2ba = 2.*lhat*v_sqr;
   const double c2ba = -(mA2 - lA*v_sqr);
   const double bma = 0.5*std::atan2(s2ba, c2ba);
   const double alpha_H = model.get_beta() - bma;

   CHECK_CLOSE(model.get_alpha(), alpha_H, eps);
}
